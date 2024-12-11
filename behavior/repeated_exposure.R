setwd("~/Documents/behavior-water-experiments")

library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(readxl)
library(dplyr) 
library(survival)
library(survminer)
library(ggpubr)

# Read in the data
data <- read_excel("repeated_exposure.xlsx")
data <- as.data.frame(data)

# Check the structure of the dataset
str(data)

# Fit a linear mixed-effects model
# Fixed effects: Ecotype, Trial
# Random effects: Individual ID
model <- lmer(time_seconds ~ ecotype * trial + (1 | individual_id), data = data)

# RANDOM EFFECTS
# Variance: 315.7 -> There is some variation in baseline emergence time across individuals.

# FIXED EFFECTS
# SW beetles take 126.62 seconds longer to emerge than LW beetles (significant p = 0.0458)
# Each subsequent trial adds 16.08 seconds to emergence time on average (not significant p = 0.4106).
# The effect of trial does not differ significantly between ecotypes.

# The interaction term is not significant (p = 0.7776), this means that the change in emergence time across trials does not differ between SW and LW beetles.

# There is no significant effect of trial or interaction between ecotype and trial, suggesting no strong evidence that the beetles habituate or change their behavior over time, regardless of ecotype

# Summary of the model
summary(model)

# Perform ANOVA to test significance of fixed effects
anova(model)

# Post-hoc pairwise comparisons (if needed)
emmeans_results <- emmeans(model, pairwise ~ ecotype * trial)
emmeans_results

# Check residuals for model fit
plot(residuals(model))

# Try another model, this time considering both trial and ID as random effects
model_var <- lmer(time_seconds ~ trial + (1 + trial | individual_id), data = data)
summary(model_var)

#	Beetles exhibit variability in their responses, but on average, repeated exposure to flooding does not significantly 
# affect emerge times.

# On average, the emerge time increases by 9.26 seconds per trial, but this effect is not significant (p = 0.54668), meaning there is no strong evidence that trials systematically affect emerge time across all beetles.

# The residual variance is high, meaning that many factors influencing emerge time remain unexplained by this model (e.g., stress, physical condition, or other environmental variables).

# This model highlights significant individual differences in responses to trials, with the fixed effect of trial showing no significant population-wide trend. 
# However, the random effects indicate that trial responses are not uniform across individuals, suggesting that beetle-specific factors strongly influence emerge times.

# Extract random effects
ranef_effects <- ranef(model_var)

# View random slopes for trial
ranef_effects$individual_id

ggplot(filtered_data, aes(x = trial, y = time_seconds, group = individual_id, color = as.factor(individual_id))) +
  geom_point() +
  geom_line(aes(y = predict(model_var))) +
  labs(title = "Individual Responses to Repeated Exposure to Water",
       x = "Trial",
       y = "Emerge Time (seconds)",
       color = "Individual ID") +
  theme_minimal() + facet_wrap(~ ecotype)

AIC(model, model_var)  # model is better than model_var
BIC(model, model_var)  # Compare BIC for both models

# Boxplot of emerge time by ecotype
ggplot(data, aes(x = ecotype, y = time_seconds, fill = ecotype)) +
  geom_boxplot() +
  facet_wrap(~ trial) +  # Separate by trial
  labs(
    title = "Response to repeated exposure to flooding",
    x = "Ecotype",
    y = "Emerge Time (seconds)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("SW" = "#4E95D9", "LW" = "#E97132"))

# Line plot of emerge time over trials
ggplot(data, aes(x = trial, y = time_seconds, color = ecotype, group = individual_id)) +
  geom_line(alpha = 0.5) +  # Individual trends
  stat_summary(fun = mean, geom = "line", aes(group = ecotype), size = 1) +
  labs(
    title = "Emerge Time Trends Over Trials",
    x = "Trial",
    y = "Emerge Time (seconds)",
    color = "Ecotype"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("SW" = "#4E95D9", "LW" = "#E97132"))

# Violins + boxplots + p-values
ggplot(data, aes(x = ecotype, y = time_seconds, fill = ecotype)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
  stat_compare_means(
    aes(label = ..p.signif..), # Display significance levels (e.g., *, **, ***)
    method = "wilcox.test",        
    label.y = 320            
  ) +
  facet_wrap(~ trial) +  # Separate by trial
  labs(
    title = "Response to Repeated Exposure to Flooding",
    x = "Ecotype",
    y = "Emerge Time (seconds)"
  ) +
  theme_minimal() +
  # Custom fill colors for ecotypes
  scale_fill_manual(values = c("SW" = "#4E95D9", "LW" = "#E97132")) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),               
    axis.title = element_text(size = 12)              
  )

# Add model predictions

filtered_data <- na.omit(data)
filtered_data$pred_model <- predict(model)
filtered_data$pred_model_var <- predict(model_var)

# Visualization
ggplot(filtered_data, aes(x = trial, y = time_seconds, group = individual_id, color = ecotype)) +
  geom_point() +
  geom_line(aes(y = pred_model), linetype = "dashed") +         # Predictions from model
  geom_line(aes(y = pred_model_var), linetype = "solid") +      # Predictions from model_var
  facet_wrap(~ ecotype) +
  theme_minimal()

# Add a grouping variable for line type
ggplot(filtered_data, aes(x = trial, y = time_seconds, group = individual_id, color = ecotype)) +
  geom_point() +
  geom_line(aes(y = pred_model, linetype = "model")) +         # Dashed line for model
  geom_line(aes(y = pred_model_var, linetype = "model_var")) + # Solid line for model_var
  facet_wrap(~ ecotype) +
  scale_linetype_manual(
    name = "Model Type",                    # Legend title
    values = c("model" = "dashed", "model_var" = "solid")  # Map line types
  ) +
  theme_minimal()
