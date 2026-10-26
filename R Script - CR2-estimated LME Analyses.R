### CRSE-CORRECTED LME ANALYSES (Robust Standard Errors) ###
# ---------------------------------------------------------
# If standard LME model residuals STILL show heteroscedasticity, you can keep the LME model 
# exactly as it is, but calculate "Robust Standard Errors". 
# This mathematically corrects your p-values so they are accurate even if homoscedasticity is violated. 

# References:
# Pustejovsky, J. E. (2026). clubSandwich: Cluster-Robust (Sandwich) Variance Estimators with Small-Sample Corrections. In R package version 0.6.2.9999. https://github.com/jepusto/clubsandwich
# Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. Journal of Business & Economic Statistics, 36(4), 672-683. 

# Install package if you haven't: install.packages("clubSandwich")
library(clubSandwich)

# We're still using HADS-A as example outcome score here. 
# 1. Create a "clean" dataset without NAs for this specific outcome
# This permanently prevents the clubSandwich row-mismatch errors
clean_LME_data <- subset(LME_longdata, !is.na(hadsA_scores))

# 1. Fit your standard model (even if assumptions are violated)
hadsA_lme4_cov_slope <- lmer(
  hadsA_scores ~ group*timepoint + study_site + gender + relationship + (1 + time_in_months | id_demografik),
  data = clean_LME_data, REML=TRUE,
  na.action = na.exclude
)
summary(hadsA_lme4_cov_slope)
print("--- Model B ('lme4') 95% Confidence Intervals ---")
print(confint(hadsA_lme4_cov_slope, method = "Wald")) 

# Get the "Type III Tests of Fixed Effects" table
# Note: We use 'anova()' (lowercase a) for lmerTest models
print("--- Model B ('lme4') Type III Tests of Fixed Effects ---")
print(anova(hadsA_lme4_cov_slope, type = 3, ddf = "Satterthwaite"))

# Step 6: Get EMMeans, SDs, SEs, and Run Post-Hoc Tests
# 1. Generate the emmeans object (without printing it immediately)
emm_lme <- emmeans(hadsA_lme4_cov_slope, ~ timepoint | group)

# 2. Convert the emmeans output into a workable data frame
emm_df <- as.data.frame(emm_lme)

# 3. Calculate the EXACT sample size (n) for each group at each timepoint
n_counts <- clean_LME_data %>%
  group_by(group, timepoint) %>%
  summarise(n = n(), .groups = "drop")

# 4. Merge the counts, calculate the SD, and format to exactly 2 decimal places
emm_formatted <- emm_df %>%
  left_join(n_counts, by = c("group", "timepoint")) %>%   # Join the n counts to the emmeans dataframe
  mutate(SD = SE * sqrt(n)) %>%   # Calculate Model-Implied SD (SE * sqrt(n))
  select(group, timepoint, n, emmean, SD, SE, df, lower.CL, upper.CL) %>%   # Reorganize columns so it's easy to read
  mutate(across(c(emmean, SD, SE, df, lower.CL, upper.CL), ~ sprintf("%.2f", .))) # Force all numeric columns to 2 d.p.

# 5. Print the beautifully formatted table!
print("--- Model-Implied Means and SDs (Formatted to 2 Decimal Places) ---")
print(emm_formatted)

# 6. Run your post-hoc comparisons as usual
print("--- Post-hoc Pairwise Comparisons ---")
emm_pairs_time <- emmeans(hadsA_lme4_cov_slope, pairwise ~ timepoint | group)
print(summary(emm_pairs_time, adjust = "bonferroni", infer = c(TRUE, TRUE))) 
# The 'infer = c(TRUE, TRUE)' argument explicitly requests both CIs and tests.

emm_pairs_group <- emmeans(hadsA_lme4_cov_slope, pairwise ~ group | timepoint)
print(summary(emm_pairs_group, adjust = "bonferroni", infer = c(TRUE, TRUE)))


#### CR2 CLUSTER-ROBUST ESTIMATION ####

# 2. Use the 'coef_test' function with "CR2" (Cluster-Robust) estimation.
# This clusters the variance around 'id_demografik' and corrects the p-values 
# for heteroscedasticity. Report THESE p-values instead of the standard summary() ones.

# By omitting the 'cluster' argument, clubSandwich automatically extracts 
# the 'id_demografik' clusters directly from the lmer model, ensuring a perfect match!
print("--- Robust Fixed Effects (Adjusted for Heteroscedasticity) ---")
robust_fe <- coef_test(hadsA_lme4_cov_slope, vcov = "CR2")
print(robust_fe)

# Generate Robust Type III ANOVA Table (Omnibus Tests)
print("--- Robust Type III ANOVA Table (Main Effects & Interactions) ---")
V_robust <- vcovCR(hadsA_lme4_cov_slope, type = "CR2") # Create the robust variance-covariance matrix
robust_anova <- Anova(hadsA_lme4_cov_slope, type = 3, test.statistic = "Chisq", vcov. = V_robust) # Feed the robust matrix into the car::Anova function
print(robust_anova)

# 5. Generate Post-Hoc Comparisons using the Robust standard errors
print("--- Robust Post-Hoc Comparisons (Web vs eBook over Time) ---")
emm_CR2_time <- emmeans(hadsA_lme4_cov_slope, pairwise ~ timepoint | group, vcov. = V_robust) # Feed robust matrix into emmeans
print(summary(emm_CR2_time, adjust = "bonferroni", infer = c(TRUE, TRUE))) # Print results with Bonferroni adjustment

emm_CR2_group <- emmeans(hadsA_lme4_cov_slope, pairwise ~ group | timepoint, vcov. = V_robust)
print(summary(emm_CR2_group, adjust = "bonferroni", infer = c(TRUE, TRUE)))

