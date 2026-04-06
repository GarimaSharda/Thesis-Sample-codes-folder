################################################################################
# ADVANCED SURVIVAL ANALYSIS: LA LIGA CAREER LENGTH
# Description: Evaluates career span using both static baseline models and 
#              dynamic time-dependent covariates (counting process format).
# Datasets:    la_liga_stats...csv, liga_injuries...csv, liga_history...csv
################################################################################

# ==============================================================================
# SECTION 1: Environment Setup & Package Loading
# ==============================================================================
# Load necessary libraries
library(dplyr)      # Data manipulation
library(tidyr)      # Reshaping data (replace_na, etc.)
library(stringr)    # String parsing and cleaning
library(lubridate)  # Date parsing for injury tracking
library(readr)      # Fast file reading
library(ggplot2)    # Core visualization
library(patchwork)  # Combining ggplot panels (p1 | p2 / p3)
library(survival)   # Core survival models (Cox, AFT, cluster, strata, tt, pspline)
library(survminer)  # Enhanced survival model visualizations

# Set seed for reproducibility
set.seed(123)

# Define working directory (Update this path to where your CSVs are located)
# setwd("C:/Users/garim/OneDrive/Documents") 


# ==============================================================================
# SECTION 2: Data Ingestion
# ==============================================================================
cat("\n--- Loading Datasets ---\n")
df_static    <- read_csv("la_liga_stats_2010_2015_depured2.csv", show_col_types = FALSE)
df_injuries  <- read_csv("liga_injuries_detailed.csv", show_col_types = FALSE)
df_history   <- read_csv("liga_history_detailed.csv", show_col_types = FALSE)


# ==============================================================================
# SECTION 3: Static Data Engineering (Baseline Profile)
# ==============================================================================
cat("\n--- Processing Static Data ---\n")
# Creates the dataset for standard Kaplan-Meier and Baseline Models.
# IMPORTANT: Keeps continuous variables intact and converts aggregates to rates.
surv_df <- df_static %>%
  select(
    player_name = player_name,
    player_url = player_url_fbref, # Unique ID
    nationality = citizenship_tra, 
    career_length,
    retired = retired_tra, 
    debut_age = debut_age_transfermarket,
    debut_year = debut_year_fbref,
    injury_count_total = injury_count_total_tra,
    seasons_club = total_seasons_national,
    foot = foot_tra,
    position_simple,
    height = height_tra
  ) %>%
  mutate(
    # Clean continuous numeric data (No arbitrary binning)
    debut_age_num  = as.numeric(gsub("[^0-9.]", "", debut_age)),
    height_clean   = as.numeric(height),
    debut_year_num = as.numeric(debut_year),
    
    # [ANTI-LEAKAGE] Calculate Intensity Rates for static models
    injuries_per_season = ifelse(!is.na(seasons_club) & seasons_club > 0, 
                                 injury_count_total / seasons_club, 0),
    
    # Clean categorical covariates
    foot_clean = as.factor(case_when(
      tolower(foot) %in% c("left", "right", "both") ~ tolower(foot),
      TRUE ~ "unknown"
    )),
    
    # Binarize Event
    event = case_when(
      retired == "Yes" ~ 1L,
      retired == "No"  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # Group Nationalities for Frailty/Strata Models
    nationality_clean = str_squish(str_replace_all(nationality, "\u00A0", " ")), 
    region = as.factor(case_when(
      str_detect(nationality_clean, "England|Scotland|Wales|Northern Ireland|Ireland") ~ "UK_Ireland",
      str_detect(nationality_clean, "France|Spain|Germany|Italy|Netherlands|Portugal|Belgium") ~ "EU_Major",
      str_detect(nationality_clean, "Brazil|Argentina|Uruguay") ~ "South_America",
      TRUE ~ "Rest_of_World"
    ))
  ) %>%
  # Filter only complete records for accurate model comparisons
  filter(complete.cases(career_length, event, debut_age_num, debut_year_num, 
                        position_simple, height_clean, region, injuries_per_season))


# ==============================================================================
# SECTION 4: Exploratory Data Analysis (EDA)
# ==============================================================================
cat("\n--- Generating Exploratory Plots ---\n")

# 4.1 Static Variables Dashboard
p_pos <- ggplot(surv_df, aes(x = position_simple)) +
  geom_bar(fill = "firebrick", alpha = 0.8) +
  theme_minimal() + labs(title = "Positions", x = "", y = "Count")

p_age <- ggplot(surv_df, aes(x = debut_age_num, y = career_length)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "loess", se = FALSE, colour = "red") +
  theme_minimal() + labs(title = "Debut Age vs Career Length", x = "Debut Age", y = "Years")

print(p_pos | p_age)

# 4.2 Injury Dashboard
eda_injuries <- df_injuries %>%
  mutate(
    severity_cat = case_when(
      days <= 7 ~ "Minor (≤1w)", days <= 28 ~ "Moderate (1-4w)",
      days <= 90 ~ "Serious (1-3m)", TRUE ~ "Severe (>3m)"
    ),
    injury_group = case_when(
      grepl("muscle|hamstring|thigh|calf|groin", tolower(injury_type)) ~ "Muscular",
      grepl("knee|ankle|ligament|acl|meniscus", tolower(injury_type)) ~ "Joint",
      grepl("fracture|bone|break", tolower(injury_type)) ~ "Bone",
      TRUE ~ "Other"
    )
  )

p_inj1 <- ggplot(eda_injuries %>% count(injury_group), aes(x = reorder(injury_group, n), y = n)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() + labs(title = "Type", x = "", y = "Count")
p_inj2 <- ggplot(eda_injuries %>% count(severity_cat), aes(x = severity_cat, y = n, fill = severity_cat)) +
  geom_col() + theme_minimal() + theme(legend.position = "none") + labs(title = "Severity", x = "", y = "")

print(p_inj1 | p_inj2)

# ==============================================================================
# SECTION 4: Exploratory Data Analysis (EDA) - Univariate & Bivariate
# ==============================================================================
cat("\n--- Generating Univariate & Bivariate Exploratory Plots ---\n")

# ------------------------------------------------------------------------------
# 4.1 Univariate Analysis (Distributions of Individual Variables)
# ------------------------------------------------------------------------------

# A. Continuous Numeric Variables (Histograms)
num_vars <- c("career_length", "debut_age_num", "height_clean", "injuries_per_season")

for (v in num_vars) {
  if (v %in% names(surv_df)) {
    p <- ggplot(surv_df, aes(x = .data[[v]])) +
      geom_histogram(bins = 30, colour = "white", fill = "steelblue") +
      theme_minimal() +
      labs(title = paste("Univariate Distribution of", v), x = v, y = "Frequency")
    print(p) # Prints each histogram individually
  }
}

# B. Categorical Variables (Bar Charts)
p_pos <- ggplot(surv_df, aes(x = position_simple)) +
  geom_bar(fill = "firebrick", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Univariate: Player Count by Position", x = "Position", y = "Count")
print(p_pos)

p_foot <- ggplot(surv_df, aes(x = foot_clean)) +
  geom_bar(fill = "forestgreen", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Univariate: Player Count by Preferred Foot", x = "Foot", y = "Count")
print(p_foot)

p_region <- ggplot(surv_df, aes(x = region)) +
  geom_bar(fill = "purple", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Univariate: Player Count by Region", x = "Region", y = "Count")
print(p_region)


# ------------------------------------------------------------------------------
# 4.2 Bivariate Analysis (Predictors vs. Target Variable: Career Length)
# ------------------------------------------------------------------------------

# A. Categorical Predictors vs Career Length (Boxplots)
p_biv_pos <- ggplot(surv_df, aes(x = position_simple, y = career_length, fill = position_simple)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Bivariate: Career Length by Position", x = "Position", y = "Career Length (Years)")
print(p_biv_pos)

p_biv_foot <- ggplot(surv_df, aes(x = foot_clean, y = career_length, fill = foot_clean)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Bivariate: Career Length by Preferred Foot", x = "Foot", y = "Career Length (Years)")
print(p_biv_foot)

# B. Continuous Predictors vs Career Length (Scatterplots with Trend Lines)
p_biv_age <- ggplot(surv_df, aes(x = debut_age_num, y = career_length)) +
  geom_point(alpha = 0.3, colour = "darkgray") +
  geom_smooth(method = "loess", se = TRUE, colour = "red") +
  theme_minimal() +
  labs(title = "Bivariate: Career Length vs Debut Age", x = "Debut Age", y = "Career Length (Years)")
print(p_biv_age)

p_biv_height <- ggplot(surv_df, aes(x = height_clean, y = career_length)) +
  geom_point(alpha = 0.3, colour = "darkgray") +
  geom_smooth(method = "loess", se = TRUE, colour = "blue") +
  theme_minimal() +
  labs(title = "Bivariate: Career Length vs Height", x = "Height (cm)", y = "Career Length (Years)")
print(p_biv_height)

p_biv_inj <- ggplot(surv_df, aes(x = injuries_per_season, y = career_length)) +
  geom_point(alpha = 0.3, colour = "darkgray") +
  geom_smooth(method = "loess", se = TRUE, colour = "darkgreen") +
  theme_minimal() +
  labs(title = "Bivariate: Career Length vs Injury Rate", x = "Injuries per Season", y = "Career Length (Years)")
print(p_biv_inj)


# ------------------------------------------------------------------------------
# 4.3 Injury-Specific EDA
# ------------------------------------------------------------------------------
# Create a temporary dataframe to visualize the raw injury types and severity
eda_injuries <- df_injuries %>%
  mutate(
    severity_cat = case_when(
      days <= 7 ~ "Minor (≤1 week)",
      days <= 28 ~ "Moderate (1-4 weeks)",
      days <= 90 ~ "Serious (1-3 months)",
      TRUE ~ "Severe (>3 months)"
    ),
    injury_group = case_when(
      grepl("muscle|hamstring|thigh|calf|groin", tolower(injury_type)) ~ "Muscular",
      grepl("knee|ankle|ligament|acl|meniscus", tolower(injury_type)) ~ "Joint/Ligament",
      grepl("fracture|bone|break", tolower(injury_type)) ~ "Bone",
      TRUE ~ "Other"
    )
  )

# Plot Distribution of Injury Types
p_inj_type <- ggplot(eda_injuries %>% count(injury_group), aes(x = reorder(injury_group, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Total Count of Injury Types", x = "Injury Category", y = "Frequency")
print(p_inj_type)

# Plot Distribution of Injury Severity
p_inj_sev <- ggplot(eda_injuries %>% count(severity_cat), aes(x = severity_cat, y = n, fill = severity_cat)) +
  geom_col() +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total Count of Injuries by Severity", x = "Severity", y = "Frequency")
print(p_inj_sev)



# ==============================================================================
# SECTION 5: Kaplan-Meier Survival Estimates
# ==============================================================================
surv_static <- Surv(time = surv_df$career_length, event = surv_df$event)

km_pos <- survfit(surv_static ~ position_simple, data = surv_df)
p_km <- ggsurvplot(km_pos, data = surv_df, pval = TRUE, risk.table = TRUE, 
                   title = "Survival Probability by Position")
print(p_km)

# ==============================================================================
# SECTION 5.2: Additional Kaplan-Meier Curves
# ==============================================================================
# Temporarily create categorical bands just for the KM plots
km_df <- surv_df %>%
  mutate(
    debut_band = cut(debut_age_num, breaks = c(15, 18, 21, 24, 30), right = FALSE),
    debut_decade = as.factor(paste0(floor(debut_year_num / 10) * 10, "s"))
  )

# Survival by Debut Age
km_debut <- survfit(surv_static ~ debut_band, data = km_df)
p_km_debut <- ggsurvplot(km_debut, data = km_df, pval = TRUE, 
                         title = "Survival by Debut Age Band")
print(p_km_debut)

# Survival by Era/Decade (Are modern players lasting longer?)
km_decade <- survfit(surv_static ~ debut_decade, data = km_df)
p_km_decade <- ggsurvplot(km_decade, data = km_df, pval = TRUE, xlim = c(0, 15),
                          title = "Career Survival by Debut Decade", palette = "jco")
print(p_km_decade)



# ==============================================================================
# SECTION 6: Static Baseline Models (Cox & AFT)
# ==============================================================================
cat("\n--- Running Static Models ---\n")

# 6.1 Cox Proportional Hazards (Using Splines for continuous height and age)
cox_static <- coxph(Surv(career_length, event) ~ position_simple + pspline(debut_age_num) + 
                      pspline(height_clean) + region + foot_clean + injuries_per_season, 
                    data = surv_df)
print(summary(cox_static))

# 6.2 Accelerated Failure Time (Weibull) - Robustness Check for PH violations
aft_weibull <- survreg(Surv(career_length, event) ~ position_simple + debut_age_num + 
                         height_clean + region + injuries_per_season, 
                       data = surv_df, dist = "weibull")
print(summary(aft_weibull))


# ==============================================================================
# SECTION 6.1: Static Cox Model & Stepwise AIC Selection
# ==============================================================================
cat("\n--- Running Static Cox Model & AIC Selection ---\n")

# Define the "Full" Model with all potential covariates
cox_static_full <- coxph(Surv(career_length, event) ~ position_simple + pspline(debut_age_num) + 
                           pspline(height_clean) + region + foot_clean + injuries_per_season, 
                         data = surv_df)

# Perform Stepwise Covariate Selection based on AIC (as requested by committee)
cox_static_aic <- step(cox_static_full, direction = "both", trace = 0)

cat("\n--- Final Selected Static Cox Model (via AIC) ---\n")
print(summary(cox_static_aic))

# ==============================================================================
# SECTION 7: Dynamic Time-Dependent Engineering (Counting Process)
# ==============================================================================
cat("\n--- Engineering Time-Dependent Covariates ---\n")

# 7.1 Process Dynamic Injury Data by Season
injuries_td <- eda_injuries %>%
  mutate(season_year = as.numeric(season_simple)) %>%
  group_by(player_url, season_year) %>%
  summarise(
    n_injuries = n(),
    total_days_missed = sum(days, na.rm = TRUE),
    serious_injuries = sum(severity_cat %in% c("Serious (1-3m)", "Severe (>3m)")),
    .groups = "drop"
  )

# 7.2 Process Performance Workload by Season
performance_td <- df_history %>%
  mutate(season_year = as.numeric(str_extract(season_simple, "^[0-9]{4}"))) %>%
  group_by(player_url, season_year) %>%
  summarise(matches_played = sum(matches_played, na.rm = TRUE), .groups = "drop")

# 7.3 Merge Base History, Apply Age-Based Time Scale (Fixes Left Truncation)
td_df_full <- df_history %>%
  mutate(season_year = as.numeric(str_extract(season_simple, "^[0-9]{4}"))) %>%
  distinct(player_url, season_year) %>%
  
  # Join metrics
  left_join(injuries_td, by = c("player_url", "season_year")) %>%
  left_join(performance_td, by = c("player_url", "season_year")) %>%
  mutate(across(c(n_injuries, total_days_missed, serious_injuries, matches_played), ~replace_na(., 0))) %>%
  
  # Join static data to calculate Age
  left_join(surv_df %>% select(player_url, debut_year_num, debut_age_num, event_static = event, 
                               position_simple, height_clean, region, foot_clean), by = "player_url") %>%
  filter(!is.na(debut_age_num), !is.na(position_simple)) %>% # Remove unmatched players
  
  arrange(player_url, season_year) %>%
  group_by(player_url) %>%
  mutate(
    # Rolling Workload and Injury Metrics
    cum_injuries = lag(cumsum(n_injuries), default = 0),
    had_serious_last_2y = (lag(serious_injuries, 1, default = 0) + lag(serious_injuries, 2, default = 0)) > 0,
    cum_matches = lag(cumsum(matches_played), default = 0),
    high_load = lag(matches_played, 1, default = 0) > 30,
    
    # METHODOLOGY FIX: Age as the Time Scale (Solves Left Truncation)
    age_start = debut_age_num + (season_year - debut_year_num),
    age_stop  = age_start + 1,
    
    # Event only triggers on the final row if the player officially retired
    is_last_row = row_number() == n(),
    event_dynamic = ifelse(is_last_row & event_static == 1, 1, 0)
  ) %>%
  ungroup() %>%
  # Ensure no logical errors in time computation
  filter(age_stop > age_start)


# ==============================================================================
# SECTION 8: Final Time-Dependent Cox Model (Methodology Corrected)
# ==============================================================================
cat("\n--- Running Final Time-Dependent Cox Model ---\n")

# This final model uses:
# 1. Age as the underlying time scale (Surv(age_start, age_stop)).
# 2. Continuous height with penalized splines to handle U-shaped physical relationships.
# 3. Stratification by region to correct spatial PH violations.
# 4. Time-Transform tt() functions to correct temporal PH violations for workload metrics.

cox_td_final <- coxph(
  Surv(time = age_start, time2 = age_stop, event = event_dynamic) ~ 
    
    # Baseline Characteristics
    position_simple + 
    foot_clean + 
    pspline(height_clean, df = 3) + 
    strata(region) + # Fixes PH violation for Region
    
    # Dynamic Injury States
    n_injuries + 
    cum_injuries + 
    had_serious_last_2y +
    
    # Dynamic Workload (Base Effects)
    matches_played + 
    cum_matches + 
    high_load +
    
    # Time-Varying Coefficients (Fixes PH violations for continuous workload)
    # Allows the protective effect of playing matches to diminish as age increases
    tt(matches_played) + 
    tt(cum_matches) + 
    tt(high_load) +
    
    # Adjust standard errors for multiple rows per player
    cluster(player_url), 
  
  data = td_df_full, 
  id = player_url,
  tt = function(x, t, ...) x * log(t + 1)
)

print(summary(cox_td_final))

# Note: cox.zph(cox_td_final) is intentionally omitted because global Schoenfeld 
# residuals are mathematically undefined when using explicit tt() time-transform terms.

