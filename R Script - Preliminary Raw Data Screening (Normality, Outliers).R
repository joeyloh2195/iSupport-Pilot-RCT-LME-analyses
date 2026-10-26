# R Script for Preliminary Raw Data Screening (Outlier detection, Normality deviation)
# ==============================================================================
# Script for long format data
# ==============================================================================

#load packages to read .sav file (SPSS file)
library(tidyverse) # Loads ggplot2, dplyr, etc.
library(psych)
library(rstatix)

# Step 1: Load Your SPSS Data (ensure the data is already in long data format)
LME_longdata <- read_sav(file.choose())
View(LME_longdata)

# Step 2: Pre-process the Data 
# Assuming 'timepoint' is coded 0 (baseline), 1 (1-month), 2 (3-month)
LME_longdata <- LME_longdata %>%
  mutate(
    # 1. *** NEW: Create a true numeric time variable for unequal intervals ***
    # This reflects the actual time in months: 0, 1, and 3
    time_in_months = case_when(
      timepoint == 0 ~ 0,
      timepoint == 1 ~ 1,
      timepoint == 2 ~ 3
    ),
    # 2. Convert categorical variables to factors (remove/ add variable as appropriate)
    id_demografik = factor(id_demografik),
    group = factor(group, levels = c(1, 2), labels = c("Web", "eBook")),
    study_site = factor(study_site, levels = c(1, 2, 3, 4), labels = c("UMMC", "SGH", "SHC", "Social")),
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    relationship = factor(relationship, levels = c(1, 2), labels = c("Child", "Others")),
    # This 'timepoint' factor is for the FIXED effects (to get estimates for t0, t1, t2)
    timepoint = factor(timepoint, levels = c(0, 1, 2), labels = c("t0", "t1", "t2"))
  )


## ---------------------------------------------------------
# 1. OUTLIER DETECTION (Z-SCORES & BOXPLOTS)
# ---------------------------------------------------------
print("--- Flagging Extreme Univariate Outliers (Z > |3.29|) ---")
# We calculate Z-scores WITHIN each group and timepoint
LME_longdata <- LME_longdata %>%
  group_by(group, timepoint) %>% 
  mutate(
    mean_score = mean(hadsA_scores, na.rm = TRUE),
    sd_score = sd(hadsA_scores, na.rm = TRUE)
  ) %>%
  mutate(
    z_score = (hadsA_scores - mean_score) / sd_score
  ) %>%
  ungroup()

# Now simply filter the updated dataset to find the outliers
z_score_check <- LME_longdata %>%
  dplyr::filter(abs(z_score) > 3) %>%  # Tabachnick & Fidell cutoff > 3.29 
  select(id_demografik, group, timepoint, hadsA_scores, z_score)

print("--- OUTLIERS DETECTED VIA Z-SCORE (> 3 SD) ---")
print(z_score_check)


# Generate the Boxplot to visually inspect (using jitters command & z-score outlier detection)
boxplot_raw <- ggplot(LME_longdata, 
                      aes(x = timepoint, y = hadsA_scores, fill = group)) +
  # Hide default boxplot outliers to prevent double-plotting with jitter
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  # Add jittered points for all participants
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             alpha = 0.4, size = 1.5, color = "black") +
  # Highlight and label the extreme outliers identified by Z-scores
  geom_text(data = z_score_check,
            aes(label = id_demografik, y = hadsA_scores),
            position = position_dodge(width = 0.75), # Align with the correct box
            vjust = -1, color = "magenta", fontface = "bold", size = 3) +
  labs(title = "Raw Scores Boxplot: Identifying Univariate Outliers",
       subtitle = "Red labels indicate extreme outliers (Z > 3.29)",
       #"; Pink labels indicate potential outliers (Z > 3.00)",
       x = "Timepoint", y = "HADS-A Scores") +
  theme_minimal()
print(boxplot_raw)


# ---------------------------------------------------------
# 2. NORMALITY TESTING (SKEWNESS, KURTOSIS, SHAPIRO, LILLIEFORS)
# ---------------------------------------------------------
print("--- Descriptive Stats: Skewness & Kurtosis ---")
# psych::describeBy automatically calculates skewness and excess kurtosis by group!
# We create a combined grouping variable first:
LME_longdata$group_time <- paste(LME_longdata$group, LME_longdata$timepoint, sep="_")

desc_stats <- describeBy(LME_longdata$hadsA_scores, group = LME_longdata$group_time, mat = TRUE) %>%
  select(group1, n, mean, sd, skew, kurtosis) %>%
  mutate(
    # 1. Calculate the Standard Errors based on subgroup sample size (n)
    se_skew = sqrt(6 / n),
    se_kurt = sqrt(24 / n),
    # 2. Calculate the Z-scores for skewness and kurtosis
    z_skew = skew / se_skew,
    z_kurt = kurtosis / se_kurt,
    # 3. Flag violations: If absolute Z-score > 1.96, it violates normality
    skew_violates = abs(z_skew) > 1.96,
    kurt_violates = abs(z_kurt) > 1.96
  ) %>%
  # Clean up numeric formatting to 3 decimal places for readability
  mutate(across(where(is.numeric), ~ round(., 3)))

print(desc_stats)

print("--- Shapiro-Wilk Test by Group & Timepoint ---")
# Using rstatix to easily run Shapiro across all subgroups
shapiro_results <- LME_longdata %>%
  group_by(group, timepoint) %>%
  shapiro_test(hadsA_scores)
print(shapiro_results)

print("--- Lilliefors Test by Group & Timepoint ---")
# Custom dplyr summarize to run Lilliefors on each subgroup
lillie_results <- LME_longdata %>%
  group_by(group, timepoint) %>%
  summarize(
    lillie_stat = lillie.test(hadsA_scores)$statistic,
    lillie_p_value = lillie.test(hadsA_scores)$p.value,
    .groups = "drop"
  )
print(lillie_results)


# ------------------------------------------------------------------------------
# [OTHER METHODS] IQR METHOD (Boxplot Rule)
# ------------------------------------------------------------------------------
# Best for: Skewed data (clinical data) where Mean/SD might be biased.

iqr_check <- LME_longdata %>%
  group_by(timepoint) %>%
  mutate(
    Q1 = quantile(hadsA_scores, 0.25, na.rm = TRUE),
    Q3 = quantile(hadsA_scores, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_limit = Q1 - 1.5 * IQR,
    upper_limit = Q3 + 1.5 * IQR,
    is_outlier = hadsA_scores < lower_limit | hadsA_scores > upper_limit
  ) %>%
  dplyr::filter(is_outlier == TRUE) %>%
  select(id_demografik, timepoint, hadsA_scores, lower_limit, upper_limit)

print("--- OUTLIERS DETECTED VIA IQR METHOD ---")
print(iqr_check)
