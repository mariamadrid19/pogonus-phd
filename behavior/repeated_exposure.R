setwd("~/Documents/behavior-water-experiments")

# Load necessary libraries
library(lme4)  # For linear mixed-effects models
library(lmerTest)  # For p-values in mixed models
library(ggplot2)  # For visualization
library(emmeans)  # For pairwise comparisons
library(readxl)

data <- read_excel("repeated_exposure.xlsx")
data <- as.data.frame(data)

str(data)

# Fit a linear mixed-effects model
# Fixed effects: Ecotype, Trial
# Random effects: Individual ID
model <- lmer(time_seconds ~ ecotype * trial + (1 | individual_id), data = data)

# Summary of the model
summary(model)

# Perform ANOVA to test significance of fixed effects
anova(model)

# Post-hoc pairwise comparisons (if needed)
emmeans_results <- emmeans(model, pairwise ~ ecotype * trial)
emmeans_results

# Check residuals for model fit
plot(residuals(model))

# --- PLOTS ---

# Boxplot of emerge time by ecotype
ggplot(data, aes(x = ecotype, y = time_seconds, fill = ecotype)) +
  geom_boxplot() +
  facet_wrap(~ trial) +  # Separate by trial
  labs(
    title = "Emerge Time by Ecotype and Trial",
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
