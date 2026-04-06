################################################################################
# ADVANCED SURVIVAL ANALYSIS: LA LIGA CAREER LENGTH
# Description: Evaluates career span using both static baseline hazard models 
#              and dynamic time-dependent covariates (counting process format).
# Datasets:    la_liga_stats...csv, liga_injuries...csv, liga_history...csv
################################################################################

# ==============================================================================
# SECTION 1: Environment Setup & Package Loading
# ==============================================================================
library(dplyr)      # Data manipulation
library(tidyr)      # Reshaping data (replace_na, etc.)
library(stringr)    # String parsing and cleaning
library(lubridate)  # Date parsing for injury tracking
library(readr)      # Fast file reading
library(ggplot2)    # Core visualization
library(patchwork)  # Combining ggplot panels (p1 | p2) / p3
library(survival)   # Core survival models (Cox, AFT, cluster, strata)
library(survminer)  # Enhanced survival model visualizations

set.seed(123)


# ==============================================================================
# SECTION 2: Data Ingestion
# ==============================================================================
# Assuming files are in the working directory or a "data/" folder
setwd("D:/thesis data") # Update this to the folder containing your CSVs
df_static    <- read_csv("la_liga_stats_2010_2015_depured2.csv", show_col_types = FALSE)
df_injuries  <- read_csv("liga_injuries_detailed.csv", show_col_types = FALSE)
df_history   <- read_csv("liga_history_detailed.csv", show_col_types = FALSE)


# ==============================================================================
# SECTION 3: Static Data Engineering (Baseline Profile)
# ==============================================================================
# Creates the dataset for standard Kaplan-Meier and Baseline Cox Models
# Important: Converts lifetime aggregates into rates to prevent data leakage.
surv_df <- df_static %>%
  select(
    player_name = player_name,
    player_url = player_url_fbref, # Used as unique ID across tables
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
    # Clean numeric data
    debut_age_num = as.numeric(gsub("[^0-9.]", "", debut_age)),
    height_clean  = as.numeric(height),
    debut_year_num = as.numeric(debut_year),
    
    # [ANTI-LEAKAGE] Calculate Intensity Rates
    injuries_per_season = ifelse(!is.na(seasons_club) & seasons_club > 0, 
                                 injury_count_total / seasons_club, 0),
    
    # Clean categorical covariates
    foot_clean = case_when(
      tolower(foot) %in% c("left", "right", "both") ~ as.factor(tolower(foot)),
      TRUE ~ as.factor("unknown")
    ),
    event = case_when(
      retired == "Yes" ~ 1L,
      retired == "No"  ~ 0L,
      TRUE ~ NA_integer_
    ),
    debut_decade = as.factor(paste0(floor(debut_year_num / 10) * 10, "s")),
    
    # Categorical banding
    debut_band = cut(debut_age_num, breaks = c(15, 18, 21, 24, 30), include.lowest = TRUE, right = FALSE),
    height_group = cut(height_clean, breaks = quantile(height_clean, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE), include.lowest = TRUE),
    
    # Group Nationalities for Frailty Models
    nationality_clean = str_squish(str_replace_all(nationality, "\u00A0", " ")), 
    region = as.factor(case_when(
      str_detect(nationality_clean, "England|Scotland|Wales|Northern Ireland|Ireland") ~ "UK_Ireland",
      str_detect(nationality_clean, "France|Spain|Germany|Italy|Netherlands|Portugal|Belgium") ~ "EU_Major",
      str_detect(nationality_clean, "Brazil|Argentina|Uruguay") ~ "South_America",
      TRUE ~ "Rest_of_World"
    ))
  ) %>%
  # Keep only complete records required for fair Stepwise AIC comparison
  filter(complete.cases(career_length, event, debut_age_num, position_simple, height_group, region, injuries_per_season))
# ==============================================================================
# SECTION 3.5: Exploratory Data Analysis (EDA) - Static Variables
# ==============================================================================
cat("\n--- Generating Static EDA Plots ---\n")

# 1. Univariate Histograms for Numeric Variables
num_vars <- c("career_length", "debut_age_num", "injury_count_total", 
              "injuries_per_season", "height_clean")

for (v in num_vars) {
  if (v %in% names(surv_df)) {
    p <- ggplot(surv_df, aes(x = .data[[v]])) +
      geom_histogram(bins = 30, colour = "white", fill = "steelblue") +
      theme_minimal() +
      labs(title = paste("Distribution of", v), x = v, y = "Count")
    print(p)
  }
}

# 2. Categorical Bar Charts (Position and Foot)
p_pos <- ggplot(surv_df, aes(x = position_simple)) +
  geom_bar(fill = "firebrick", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Count by Position", x = "Position", y = "Count")
print(p_pos)

p_foot <- ggplot(surv_df, aes(x = foot_clean)) +
  geom_bar(fill = "forestgreen", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Count by Preferred Foot", x = "Foot", y = "Count")
print(p_foot)

# 3. Bivariate Relationships
p_biv1 <- ggplot(surv_df, aes(x = position_simple, y = career_length)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(title = "Career Length by Position", x = "Position", y = "Career Length")
print(p_biv1)

p_biv2 <- ggplot(surv_df, aes(x = debut_age_num, y = career_length)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  theme_minimal() +
  labs(title = "Career Length vs Debut Age", x = "Debut Age", y = "Career Length")
print(p_biv2)


# ==============================================================================
# SECTION 3.6: Exploratory Data Analysis (EDA) - Injuries 
# ==============================================================================
cat("\n--- Generating Injury EDA Plots ---\n")

# Process a temporary dataset for Injury EDA
eda_injuries <- df_injuries %>%
  mutate(
    severity_cat = case_when(
      days <= 7 ~ "Minor (≤1 week)",
      days <= 28 ~ "Moderate (1-4 weeks)",
      days <= 90 ~ "Serious (1-3 months)",
      days > 90 ~ "Severe (>3 months)"
    ),
    injury_group = case_when(
      grepl("muscle|hamstring|thigh|calf|groin", tolower(injury_type)) ~ "Muscular",
      grepl("knee|ankle|ligament|acl|meniscus", tolower(injury_type)) ~ "Joint/Ligament",
      grepl("fracture|bone|break", tolower(injury_type)) ~ "Bone",
      grepl("head|concussion", tolower(injury_type)) ~ "Head",
      grepl("back|spine", tolower(injury_type)) ~ "Back/Spine",
      TRUE ~ "Other"
    )
  )

# Plot 1: Distribution by Injury Group
injury_type_summary <- eda_injuries %>% count(injury_group, sort = TRUE)
p1 <- ggplot(injury_type_summary, aes(x = reorder(injury_group, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution by Type", x = "", y = "Number of injuries")

# Plot 2: Distribution by Severity
severity_summary <- eda_injuries %>% count(severity_cat)
p2 <- ggplot(severity_summary, aes(x = severity_cat, y = n, fill = severity_cat)) +
  geom_col() +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Injuries per Severity", x = "", y = "Number")

# Plot 3: Top 10 Specific Injuries
top_injuries <- eda_injuries %>% count(injury_type, sort = TRUE) %>% head(10)
p3 <- ggplot(top_injuries, aes(x = reorder(injury_type, n), y = n)) +
  geom_col(fill = "coral") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Specific Injuries", x = "", y = "Frequency")

# Combine plots using the 'patchwork' library
injury_dashboard <- (p1 | p2) / p3
print(injury_dashboard)

# ==============================================================================
# SECTION 4: Exploratory Data Analysis & Kaplan-Meier Estimates
# ==============================================================================
# 4.1 Survival Object Setup
surv_static <- Surv(time = surv_df$career_length, event = surv_df$event)

# 4.2 Kaplan-Meier by Position
km_pos <- survfit(surv_static ~ position_simple, data = surv_df)
ggsurvplot(km_pos, data = surv_df, pval = TRUE, risk.table = TRUE, 
           title = "Survival by Position", legend.title = "Position")

# 4.3 Kaplan-Meier by Era (Decade)
km_decade <- survfit(surv_static ~ debut_decade, data = surv_df)
ggsurvplot(km_decade, data = surv_df, pval = TRUE, risk.table = TRUE,
           xlim = c(0, 15), # Limit axis for comparability across cohorts
           title = "Career Survival by Debut Decade", palette = "jco")


# ==============================================================================
# SECTION 5: Static Baseline Models (Cox & Frailty)
# ==============================================================================
# 5.1 Full Cox Model
cox_static_full <- coxph(Surv(career_length, event) ~ position_simple + debut_band + 
                           height_group + region + foot_clean + debut_decade + 
                           injuries_per_season, data = surv_df)

# 5.2 Stepwise Feature Selection based on AIC
cat("\n--- Running Stepwise AIC Selection for Static Model ---\n")
cox_static_aic <- step(cox_static_full, direction = "both", trace = 0)
summary(cox_static_aic)

# 5.3 Frailty Model (Random Effects by Region)
# Adds a penalty/variance term for unobserved regional heterogeneity
cox_frailty <- coxph(Surv(career_length, event) ~ position_simple + debut_band + 
                       injuries_per_season + frailty(region), data = surv_df)
print(cox_frailty)


# ==============================================================================
# SECTION 6: Time-Dependent Dynamic Engineering (Counting Process)
# ==============================================================================
# 6.1 Process Dynamic Injury Data by Season
injuries_by_season <- df_injuries %>%
  mutate(
    from_date = dmy(from_date),
    season_year = as.numeric(season_simple),
    severity_cat = case_when(
      days <= 7  ~ "Minor (≤1 week)",
      days <= 28 ~ "Moderate (1-4 weeks)",
      days <= 90 ~ "Serious (1-3 months)",
      days > 90  ~ "Severe (>3 months)"
    ),
    injury_group = case_when(
      grepl("muscle|hamstring|thigh|calf|groin", tolower(injury_type)) ~ "Muscular",
      grepl("knee|ankle|ligament|acl|meniscus", tolower(injury_type))  ~ "Joint/Ligament",
      grepl("fracture|bone|break", tolower(injury_type)) ~ "Bone",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(player_url, season_year) %>%
  summarise(
    n_injuries = n(),
    total_days_missed = sum(days, na.rm = TRUE),
    serious_injuries = sum(severity_cat %in% c("Serious (1-3 months)", "Severe (>3 months)")),
    .groups = "drop"
  )

# 6.2 Merge Base History with Injuries
# Establish all active seasons per player, joining injury aggregates
all_seasons <- df_history %>%
  mutate(season_year = as.numeric(str_extract(season_simple, "^[0-9]{4}"))) %>%
  distinct(player_url, season_year) %>%
  left_join(injuries_by_season, by = c("player_url", "season_year")) %>%
  mutate(across(c(n_injuries, total_days_missed, serious_injuries), ~replace_na(., 0))) %>%
  arrange(player_url, season_year) %>%
  group_by(player_url) %>%
  mutate(
    # Create rolling and lagged windows for time-dependent covariates
    cum_injuries = lag(cumsum(n_injuries), default = 0),
    injuries_last_2y = lag(n_injuries, 1, default = 0) + lag(n_injuries, 2, default = 0),
    had_serious_last_2y = (lag(serious_injuries, 1, default = 0) + lag(serious_injuries, 2, default = 0)) > 0,
    
    # Establish tstart and tstop for interval censoring
    # Each row represents a season: from t=0 to t=1, t=1 to t=2, etc.
    tstart = row_number() - 1,
    tstop  = row_number(),
    
    # Event only triggers on the final row of a player's career if they actually retired
    is_last_row = row_number() == n()
  ) %>%
  ungroup()

# 6.3 Final Join to Static Baseline
td_df <- all_seasons %>%
  left_join(surv_df %>% select(player_url, event_static = event, position_simple, debut_band, height_group, region, foot_clean), by = "player_url") %>%
  mutate(event_dynamic = ifelse(is_last_row & event_static == 1, 1, 0)) %>%
  filter(!is.na(position_simple))


# ==============================================================================
# SECTION 7: Dynamic Time-Dependent Cox Models
# ==============================================================================
# 7.1 Cox Model with Rolling Injury Metrics
# Tracks the hazard of retirement updated season-to-season based on recent physical trauma
cox_td_inj <- coxph(Surv(tstart, tstop, event_dynamic) ~ 
                      position_simple + debut_band + height_group + region + foot_clean + 
                      n_injuries + cum_injuries + injuries_last_2y + had_serious_last_2y +
                      cluster(player_url), # Corrects standard errors for repeated player observations
                    data = td_df, id = player_url)

cat("\n--- Time-Dependent Cox Model (Injuries) ---\n")
summary(cox_td_inj)

# 7.2 Proportional Hazards Test for Time-Dependent Model
ph_test_td <- cox.zph(cox_td_inj)
print(ph_test_td)


# ==============================================================================
# SECTION 8: Accelerated Failure Time (AFT) Models
# ==============================================================================
# An alternative parametric approach to Cox proportional hazards.
# Predicts actual time to event rather than risk/hazard ratio.
aft_weibull_full <- survreg(Surv(career_length, event) ~ position_simple + debut_band + 
                              region + injuries_per_season, 
                            data = surv_df, dist = "weibull")

cat("\n--- AFT Weibull Model ---\n")
summary(aft_weibull_full)


# ==============================================================================
# SECTION 9: Model Architecture Comparison (AIC)
# ==============================================================================
# Compare the predictive validity of different architectural approaches.
# Note: Static AICs cannot be directly compared to Time-Dependent AICs 
# due to different dataset sizes (N=players vs N=seasons), but can be 
# compared within their own structural classes.
model_comparison <- data.frame(
  Model_Type = c("Static Baseline", "Static AIC Selected", "Frailty (Region)", "Dynamic Time-Dependent"),
  AIC = c(
    AIC(cox_static_full),
    AIC(cox_static_aic),
    AIC(cox_frailty),
    AIC(cox_td_inj)
  )
)

print(model_comparison)

# ==============================================================================
# SECTION 10: Adding Workload/Performance to Dynamic Model
# ==============================================================================
cat("\n--- Time-Dependent Cox Model (Injuries + Workload) ---\n")

# Process performance metrics from the history dataset
performance_td <- df_history %>%
  mutate(season_year = as.numeric(str_extract(season_simple, "^[0-9]{4}"))) %>%
  group_by(player_url, season_year) %>%
  summarise(matches_played = sum(matches_played, na.rm = TRUE), .groups = "drop")

# Merge into our existing Time-Dependent dataset and create rolling workload metrics
td_df_full <- td_df %>%
  left_join(performance_td, by = c("player_url", "season_year")) %>%
  mutate(matches_played = replace_na(matches_played, 0)) %>%
  group_by(player_url) %>%
  arrange(season_year) %>%
  mutate(
    cum_matches = lag(cumsum(matches_played), default = 0),
    matches_last_2y = lag(matches_played, 1, default = 0) + lag(matches_played, 2, default = 0),
    
    # Define high workload (e.g., playing more than 30 matches a season)
    high_load = lag(matches_played, 1, default = 0) > 30,
    
    # Define sharp workload increase (e.g., playing 15 more matches than the previous season)
    sharp_increase = (lag(matches_played, 1, default = 0) - lag(matches_played, 2, default = 0)) > 15
  ) %>%
  ungroup()

# Run the comprehensive model
cox_td_full <- coxph(Surv(tstart, tstop, event_dynamic) ~ 
                       position_simple + debut_band + height_group + region + foot_clean + 
                       n_injuries + cum_injuries + had_serious_last_2y +
                       matches_played + cum_matches + high_load +
                       cluster(player_url), 
                     data = td_df_full, id = player_url)

summary(cox_td_full)


# ==============================================================================
# SECTION 11: Correcting Proportional Hazards Violations
# ==============================================================================
cat("\n--- Stratified Time-Dependent Model (Fixing PH Violations) ---\n")

# Since `region` violated the PH assumption in Section 7.2 (p = 0.0038),
# we must use strata(region) instead of including it as a standard covariate.
# This allows each region to have its own underlying baseline hazard over time.

cox_td_stratified <- coxph(Surv(tstart, tstop, event_dynamic) ~ 
                             position_simple + debut_band + height_group + foot_clean + 
                             n_injuries + cum_injuries + had_serious_last_2y +
                             matches_played + cum_matches + high_load +
                             strata(region) +  # <--- Fixes the PH violation
                             cluster(player_url), 
                           data = td_df_full, id = player_url)

summary(cox_td_stratified)

# Re-test the assumptions on the stratified model
ph_test_stratified <- cox.zph(cox_td_stratified)
print(ph_test_stratified)
# ==============================================================================
# SECTION 12: Correcting Continuous PH Violations (Time-Varying Effects)
# ==============================================================================
cat("\n--- Time-Dependent Model with Time-Varying Coefficients (tt) ---\n")

# To fix the PH violation for matches_played, cum_matches, and high_load, 
# we interact them with time using the tt() function. 
# This allows the effect of workload to change as the player gets deeper into their career.

cox_td_tt <- coxph(Surv(tstart, tstop, event_dynamic) ~ 
                     position_simple + debut_band + height_group + foot_clean + 
                     n_injuries + cum_injuries + had_serious_last_2y +
                     
                     # Base variables
                     matches_played + cum_matches + high_load +
                     
                     # Time-varying interactions (allows effect to change over time)
                     tt(matches_played) + tt(cum_matches) + tt(high_load) +
                     
                     strata(region) + 
                     cluster(player_url), 
                   data = td_df_full, id = player_url,
                   
                   # Define how the interaction with time works (using log(time) is standard)
                   tt = function(x, t, ...) x * log(t + 1))

summary(cox_td_tt)


