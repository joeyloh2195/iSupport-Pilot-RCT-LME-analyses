# R Script for standard LME Analysis: Model Diagnostic Checks
# NOTE: Ensure that you have already run the LME model in the previous script - the current script is about  diagnostic checks on the LME model.

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

# influence.ME diagnostic tests for mixed effects model
# CRAN Package details: https://cran.r-project.org/web/packages/influence.ME/refman/influence.ME.html
install.packages("influence.ME")
library(influence.ME) # to check Cook's Distance, DFBETAS, significance change upon exclusion
library(HLMdiag) #to check leverage, Mahanobolis distance

#Load for Lilliefors Test (Kolmogorov-Smirnov Lilliefor Correction)
library(nortest)

# Handle the na.exclude behavior
# Because you used na.action = na.exclude in your lmer() model, resid() and fitted() 
# return 189 rows (padded with NAs for dropouts), matching your original data.
# So, we attach them to the original dataset first:
LME_longdata$.resid <- resid(hadsA_lme4_cov_slope)
LME_longdata$.fitted <- fitted(hadsA_lme4_cov_slope)

# Now we filter out the NAs to get our clean diagnostic dataset (168 rows)
model_diag <- subset(LME_longdata, !is.na(.resid))

# Plot 1: Q-Q Plot (for normality)
# Dots should fall on the line
qq_plot <- ggplot(model_diag, aes(sample = .resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals (HADS-A)", 
       x = "Theoretical Quantiles", y = "Residuals") +
  theme_minimal()
print(qq_plot)
shapiro.test(model_diag$.resid) # Shapiro-Wilk's test on residuals (n<50) - more sensitive in detecting deviations in tails
lillie.test(model_diag$.resid)
#ks.test(model_diag$.resid, "pnorm") # Kolmogrov-Smirnov's test on residuals (n>50) - more sensitive to deviations in the middle

# Plot 2: Histogram (for normality)
# Should look like a bell curve
hist_plot <- ggplot(model_diag, aes(x = .resid)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Residuals (HADS-A)", x = "Residuals", y = "Frequency") +
  theme_minimal()
print(hist_plot)

# Plot 3: Homoscedasticity (Residuals vs. Fitted)
# Should be a random "cloud" of points, NO funnel shape
cor(LME_longdata$.resid, LME_longdata$.fitted, method = "pearson", use = "complete.obs") #compute Pearson's R to check significance of regression line

rvf_plot <- ggplot(model_diag, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Residuals vs. Fitted Values (HADS-A)", 
       subtitle = "Pearson's r between .resid &.fitted = (insert r value)", 
       x = "Fitted Values", y = "Residuals") +
  theme_minimal()
print(rvf_plot)

# Plot 4: Linearity of residuals against outcomes
# should show a linear relationship
cor(LME_longdata$.resid, LME_longdata$hadsA_scores, method = "pearson", use = "complete.obs") #compute Pearson's R to check significance of regression line

rvc_plot <- ggplot(model_diag, aes(x = .resid, y = hadsA_scores)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Residuals vs. Outcome (HADS-A)", 
       subtitle = "Pearson's r between residuals & outcome = (insert r value)",
       x = "Residuals", y = "Outcome (hadsA_scores)") +
  theme_minimal()
print(rvc_plot)

#### DIAGNOSTIC CHECKS FOR INFLUENTIAL OBSERVATIONS #### 
## Detecting influential data using "influence" (iteratively deleting each outlier from the data)
## the estex object is essentially an object that will be used as an input to the subsequent functions
estex.hadsA_lme <- influence(hadsA_lme4_cov_slope, "id_demografik")

## 1. DFBETAS is a standardized measure that indicates level of influence observations have on single parameter estimates.
dfbetas(estex.hadsA_lme, parameters=c(2,3,4)) #parameters 2 = groupeBook; 3 = timepoint1; 4 = timepoint2

plot(estex.hadsA_lme,
     which="dfbetas",
     parameters=c(2,3,4),
     xlab="DFbetas",
     ylab="Participant ID")

## 2. Cook's Distance (more intuitive & useful, as it provides a summary measure for the influence a higher level unit exerts on all parameter estimates simultaneously, or a selection thereof)
# 'parameter=0 (default) = Cook's Distance calculated based on ALL parameters in the model
cooks.distance(estex.hadsA_lme,
               parameter=0, sort=TRUE)

plot(estex.hadsA_lme, which="cook",
     cutoff=.063, sort=TRUE, # cutoff calculation: 4/n (n = number of groups in the grouping factor under evaluation)
     xlab="Cook´s Distance",
     ylab="Participant ID")

## Test for changes in statistical significance (SigTest)
sigtest(estex.hadsA_lme, test = 1.96, parameters = 0, sort = FALSE, to.sort = NA)

## Use "exclude.influence" function to return an updated LME model from which the specified participant ID was deleted
## Recommend to delete influential cases one at the time and then re-evaluate the model each time
hadsA_exc <- exclude.influence(hadsA_lme4_cov_slope,"id_demografik", "SHC003") # example: SHC003 is the outlier
summary(hadsA_exc)


# ---------------------------------------------------------
# 3. LEVERAGE (Using HLMdiag)
# ---------------------------------------------------------
# Calculate leverage for the participant level (Level-2)
# This checks if any specific participant has extreme predictor values
lev_participant <- leverage(hadsA_lme4_cov_slope, level = "id_demografik")

# Convert to a data frame for easier viewing and plotting
# Extracting IDs directly from the model's random effects to ensure they match
lev_df <- data.frame(
  id_demografik = rownames(ranef(hadsA_lme4_cov_slope)$id_demografik),
  Leverage = as.numeric(lev_participant[,1])
)

# Sort to see who has the highest leverage
lev_df <- lev_df[order(-lev_df$Leverage), ]
print("--- Top 5 Participants with Highest Leverage ---")
head(lev_df, 5)

# Plot Leverage using base R (similar to your previous plots)
plot(lev_df$Leverage, 
     type = "h", 
     main = "Participant-Level Leverage",
     ylab = "Leverage Value", 
     xlab = "Index (Sorted Participants)",
     lwd = 2, col = "blue")
abline(h = 2 * mean(lev_df$Leverage), col = "red", lty = 2) # Common rule-of-thumb cutoff


# 4. MAHALANOBIS DISTANCE (Checking Random Effects Outliers)
# ---------------------------------------------------------
# Mahalanobis distance is best used here to check if a participant's combined Intercept and Time Slope are highly unusual.

# Step A: Extract the random effects (intercepts and slopes) for each participant
re_data <- ranef(hadsA_lme4_cov_slope)$id_demografik
colnames(re_data) <- c("Intercept", "Slope") # Rename for easier plotting

# Step B: Calculate the center (means) and covariance matrix of these random effects
re_means <- colMeans(re_data)
re_cov <- cov(re_data)

# Step C: Calculate Mahalanobis distance for each participant
mahal_dist <- mahalanobis(x = re_data, center = re_means, cov = re_cov)

# Step D: Put it in a dataframe
mahal_df <- data.frame(
  id_demografik = rownames(ranef(hadsA_lme4_cov_slope)$id_demografik),
  Mahalanobis = as.numeric(mahal_dist),
  Intercept = re_data$Intercept,
  Slope = re_data$Slope
)

# Step E: Determine the cutoff. 
# Mahalanobis distances roughly follow a Chi-Square distribution.
# Degrees of freedom (df) = number of random effects (2: intercept + slope)
cutoff_mahal <- qchisq(p = 0.95, df = 2) # 95% confidence cutoff
mahal_df$Outlier <- mahal_df$Mahalanobis > cutoff_mahal # Add a column to flag if they are an outlier based on the 95% cutoff

# Sort and view
mahal_df <- mahal_df[order(-mahal_df$Mahalanobis), ]
print("--- Top 5 Participants with Highest Mahalanobis Distance ---")
head(mahal_df, 5)

# Ellipse Plot of Mahalanobis Distance (using ggplot)
ggplot(mahal_df, aes(x = Intercept, y = Slope)) +
  # Draw the 95% confidence ellipse
  stat_ellipse(type = "norm", level = 0.95, fill = "orange", geom = "polygon", alpha = 0.3) +
  # Plot the individual participants
  geom_point(aes(color = Outlier), size = 3) +
  scale_color_manual(values = c("black", "red")) +
  # Label only the outliers with their participant ID
  geom_text(data = subset(mahal_df, Outlier == TRUE),
            aes(label = id_demografik), vjust = -1, color = "red", fontface = "bold") +
  labs(title = "Mahalanobis Distance: Random Effects (95% Confidence Ellipse)",
       subtitle = "Red points outside the ellipse represent participants with unusual trajectories",
       x = "Random Intercept (Individual Baseline Differences)",
       y = "Random Slope (Individual Rate of Change)") +
  theme_minimal() +
  theme(legend.position = "none")


# ==============================================================================
####  ---- Removal/ Replacement of outliers ----  ####
clean_LME_data_outl <- clean_LME_data #create a copy to be safe

# If your data is in long format (multiple rows per participant), you might need
# to specify both the Participant AND the specific Condition.

# Let's say we only want to remove it if Condition is "A"
# The condition inside the brackets uses & (AND) logic
clean_LME_data_outl[clean_LME_data_outl$id_demografik == "UM033" & 
                      clean_LME_data_outl$group == "Web", "hadsA_t1"] <- NA

### RERUN LME ANALYSES UPON REMOVING OUTLIERS ###
hadsA_lme4_cov_slope_outl <- lmer(
  hadsA_scores ~ group*timepoint + study_site + gender + relationship + (1 + time_in_months | id_demografik),
  data = clean_LME_data_outl, REML=TRUE,
  na.action = na.exclude
)
summary(hadsA_lme4_cov_slope_outl)

print("--- Model B ('lme4') 95% Confidence Intervals ---")
print(confint(hadsA_lme4_cov_slope_outl, method = "Wald")) 

# Get the "Type III Tests of Fixed Effects" table
# Note: We use 'anova()' (lowercase a) for lmerTest models
print("--- Model B ('lme4') Type III Tests of Fixed Effects ---")
print(anova(hadsA_lme4_cov_slope_outl, type = 3, ddf = "Satterthwaite"))
