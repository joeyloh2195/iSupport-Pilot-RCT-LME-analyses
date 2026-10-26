# R Script for standard LME Analysis (Fixed Effects, Random Effects, Plots)

#load packages to read .sav file (SPSS file)
library(tidyverse)
#a conflict command will come up. load the "devtools" package "install.packages("devtools")"
library(conflicted)
library(dplyr)
conflicts_prefer(dplyr::filter())
conflicts_prefer(stats::filter())
library(haven)
library(nlme)
library(emmeans)
library(lme4)
library(lmerTest) #choose this to show p-values 
library(car)
conflicts_prefer(lmerTest::lmer) #conflicts for lme test - select lmerTest as this can show p-values

#Load for plots and graphs
library(ggplot2)     
library(broom.mixed) #for diagnostic plots
library(lattice)  #for the random effects plot

# Step 2: Load Your SPSS Data (ensure the data is already in long data format)
LME_longdata <- read_sav(file.choose())
View(LME_longdata)

# Step 3: Pre-process the Data (CRITICAL)
# LME models in R work best when categorical variables are set as "factors"
# and time variables for correlation are numeric.

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

# Step 4: Run the LME Model using 'lme4'
# This model is:
# - Fixed Effects: group, timepoint, group*timepoint
# - Random Effects: Random Intercept (id_demografik) & Slope (time_in_months)
# - Covariates: study_site, gender, relationship
# - Residual Correlation: NONE (assumed to be independent)

print("--- Running Model (using 'lme4' package) ---")
# We wrap the lmer() call in 'as_lmerModLmerTest' to get p-values
#example variable using HADS-A
hadsA_lme4_cov_slope <- lmer(
  hadsA_scores ~ group*timepoint + study_site + gender + relationship + (1 + time_in_months | id_demografik),
  data = LME_longdata_outliers, REML=TRUE,
  na.action = na.exclude
)
summary(hadsA_lme4_cov_slope)

# Step 5: Get Results for 'lme4' model
print("--- Model B ('lme4') Summary (Coefficients, Random Effects) ---")
summary(hadsA_lme4_cov_slope) # This summary includes p-values

#Get 95% Confidence Intervals for 'lme4' model estimates
print("--- Model B ('lme4') 95% Confidence Intervals ---")
print(confint(hadsA_lme4_cov_slope, method = "Wald")) 

# Get the "Type III Tests of Fixed Effects" table
# Note: We use 'anova()' (lowercase a) for lmerTest models
print("--- Model B ('lme4') Type III Tests of Fixed Effects ---")
print(anova(hadsA_lme4_cov_slope, type = 3, ddf = "Satterthwaite"))

# Step 6: Get Post-Hoc Tests for Main Effects

# 6a. Time Effect (Simple Main Effects): Comparing t0 vs t1 vs t2 WITHIN each group
print("--- Post-hoc for Timepoint, Separated by Group (Simple Main Effects) ---")
em_group_time <- emmeans(hadsA_lme4_cov_slope, ~ timepoint | group)
# Run the pairwise comparisons
pairs_group_time <- pairs(em_group_time, adjust = "bonferroni")
print(summary(pairs_group_time, infer = c(TRUE, TRUE))) # The 'infer = c(TRUE, TRUE)' argument explicitly requests both CIs and tests.

# 6b. Group Effect: Comparing Web vs eBook WITHIN each timepoint
print("--- 6b: Post-Hoc Comparisons (Group within Timepoint) ---")
# We swap the formula to '~ group | timepoint'
em_time_group <- emmeans(hadsA_lme4_cov_slope, ~ group | timepoint)
pairs_time_group <- pairs(em_time_group, adjust = "bonferroni")
# Print with CIs and p-values
print(summary(pairs_time_group, infer = c(TRUE, TRUE)))

# 6c. Raw Summary Statistics (Mean and SD)
# Calculate the raw mean and SD for each group at each timepoint
print("--- 6c: Raw Summary Statistics (Mean & SD) ---")
summary_stats <- LME_longdata %>%
  group_by(group, timepoint) %>%
  summarise(
    # Use sprintf("%.2f", x) to format as a string with exactly 2 decimals
    # This forces R to show "18.50" instead of "18.5"
    Mean = sprintf("%.2f", mean(hadsA_scores, na.rm = TRUE)),
    SD = sprintf("%.2f", sd(hadsA_scores, na.rm = TRUE)),
    N = n(),
    .groups = 'drop' # Drop grouping after calculation
  )
print(as.data.frame(summary_stats))


##### RANDOM EFFECTS #####

## RANDOM INTERCEPT ##
print("--- Generating Plot 3: Random Intercepts Plot ---")

# 1. Extract the random effects
#    'condVar = TRUE' gets the standard errors for the error bars (optional but nice)
re <- ranef(hadsA_lme4_cov_slope, condVar = TRUE)
re_df <- as.data.frame(re$id_demografik)
re_df$id_demografik <- rownames(re_df)
colnames(re_df)[1] <- "intercept_deviation"

# 2. Add the 'Group' variable to the plotting data
#    We need to find out which group each participant belongs to from the original data.
group_lookup <- unique(LME_longdata[, c("id_demografik", "group")])
re_df <- merge(re_df, group_lookup, by = "id_demografik")

# 3. Create the plot with COLOR
plot_random_intercepts <- ggplot(re_df, aes(x = reorder(id_demografik, intercept_deviation), 
                                            y = intercept_deviation, 
                                            color = group)) + # This adds the color!
  
  # Add the points (random intercepts)
  geom_point(size = 2) +
  
  # Add the error bars (optional, shows uncertainty of the random effect)
  # This requires extracting the variances, which is complex. 
  # For a simple visual, points are sufficient.
  
  # Add the zero line (grand mean)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Add labels and formatting
  labs(
    title = "Random Intercepts Plot (HADS-A)",
    x = "Participant ID (Ordered by Intercept)",
    y = "Deviation from Mean (HADS-A Score)",
    color = "Intervention Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Hide the crowded ID labels
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

print(plot_random_intercepts)


## RANDOM SLOPE ##
print("--- Step 7: Generating Random Slope Plot... ---")

# 1. Bypass 'augment()' to avoid the NA row-mismatch error.
# Because you used na.action = na.exclude in your model, fitted() will correctly
# output 189 values (padding the missing dropout follow-ups with NAs).
# We can attach these directly to your original dataset!
LME_longdata$.fitted <- fitted(hadsA_lme4_cov_slope)

# 2. Filter out the NAs just for the plot so ggplot doesn't complain
#plot_data <- subset(LME_longdata, !is.na(.fitted))
plot_data <- subset(LME_longdata, !is.na(.fitted))

# 3. Create the spaghetti plot
plot_real_world_spaghetti <- ggplot(plot_data, aes(x = time_in_months, y = .fitted, 
                                                   group = id_demografik, color = group)) +
  
  # Plot all individual participant lines (the "spaghetti")
  geom_line(alpha = 0.4) + 
  # Add the "average" overall trends for Web and eBook groups
  geom_smooth(aes(group = group), method = "lm", se = FALSE, linewidth = 2, linetype = "solid") +
  
  # Add labels
  labs(
    title = "Participant Trajectories (Random Intercepts & Slopes)",
    subtitle = "Model-predicted HADS-A trajectories for each caregiver (N=63)",
    x = "Time in Months (0 = Baseline, 1 = 1-Month, 3 = 3-Month)",
    y = "Predicted HADS-A Scores",
    color = "Intervention Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# 4. Print the plot
print(plot_real_world_spaghetti)

#### EMMEANS INTERACTION PLOT (LINE GRAPH) ####
print("--- Generating Plot 1b: EMMeans Interaction Plot (Line Graph) ---")

# 1. Get the estimated marginal means from your model
# We are getting the means for the "group by timepoint" interaction
emmeans_obj <- emmeans(hadsA_lme4_cov_slope, ~ timepoint * group)
emmeans_data <- as.data.frame(emmeans_obj) #Convert the emmeans object to a data frame for plotting with ggplot

# 2. Create the plot using ggplot2
plot_emmeans_line <- ggplot(emmeans_data, aes(x = timepoint, y = emmean, color = group, group = group)) +
  # Add the lines
  geom_line(linewidth = 1) + 
  # Add the points (the means)
  geom_point(size = 3) +
  # Add the confidence intervals as error bars
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 0.8) +
  # Add labels
  labs(
    title = "Fixed Effects: Model-Predicted HADS-A Scores",
    subtitle = "ITT Sample (N=63)",
    x = "Timepoint",
    y = "Estimated Marginal Mean (hadsA_scores)",
    color = "Group"
  ) +
  # Use a clean theme
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# 3. Print the plot to the 'Plots' pane
print(plot_emmeans_line)
