setwd("~/Documents/behavior-water-experiments")

library(dplyr)
library(ggplot2)
library(car)

# Load the data
df <- read.csv("large-scale.csv", header = TRUE)
df <- as.data.frame(df)
colnames(df)[2] <- "habitat"
df$time_h <- as.numeric(df$time_h)

# Simple model without interaction or random effects
model_simple <- lm(prop.dry ~ time_h + habitat, data = df)
summary(model_simple)
# This p-value when comparing habitat type (SW vs LW) is significant (p < 0.01), indicating that the difference between the 
# SW and LW habitats is highly significant (p = 0.002897).

# Model without random effects + interaction between time and habitat
model_interaction <- lm(prop.dry ~ time_h * habitat, data = df)
summary(model_interaction)
# the interaction between time and habitat is highly significant (p = 0.00576)
# This means that the interaction between time and habitat has a significant effect 
# on the proportion of individuals on the dry side. 

# Calculate AIC for each model
AIC(model_simple)     # AIC for the simpler model, 9.171824
AIC(model_interaction)  # AIC for the model with interaction, 1.351694 -> better model 

plot(residuals(model_interaction))
hist(residuals(model_interaction))
qqnorm(residuals(model_interaction))
qqline(residuals(model_interaction), col = "red")
# model seems to fit the data quite well

# Create predictions from the model
df$predicted <- predict(model_interaction)
# Plot observed vs predicted values
ggplot(df, aes(x = time_h, y = prop.dry, color = habitat)) +
  geom_point() +
  geom_line(aes(y = predicted), size = 1, linetype = "dashed") +
  labs(
    title = "Observed vs predicted proportions (lm(prop.dry ~ time_h * habitat)",
    x = "Time (hours)",
    y = "Proportion on Dry Side",
    color = "Habitat"
  ) +
  theme_minimal()

# Line plot
ggplot(df, aes(x = time_h, y = prop.dry, group = interaction(habitat, label), color = habitat)) +
  geom_vline(xintercept = c(0, 24, 48, 72, 96), color = "lightblue", size = 3, alpha = 0.3) +
  geom_vline(xintercept = c(2, 26, 50, 74, 98), color = "lightcoral", size = 1, alpha = 0.2) +
  # Main plot elements
  geom_line(size = 1) +
  geom_point(size = 2) +
  # Labels and theme
  labs(
    title = "Spatial sorting of ecotypes in response to flooding",
    x = "Time (hours)",
    y = "Proportion on Dry Side",
    color = "Habitat",
    linetype = "Label"
  ) +
  theme_minimal()


# I don't think the data is parametric, let's check!
ggplot(df, aes(x = prop.dry)) + 
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  theme_minimal()
qqnorm(df$prop.dry)
qqline(df$prop.dry, col = "red")
shapiro.test(df$prop.dry)
# THE DATA DOES NOT FOLLOW A NORMAL DISTRIBUTION, USE NON-PARAMETRIC TESTS

# Perform Wilcoxon rank-sum test (Mann-Whitney U test) since we just determined that the data is non-parametric
wilcox_test_result <- wilcox.test(prop.dry ~ habitat, data = df)
p_value <- wilcox_test_result$p.value

# Line plot with stats
ggplot(df, aes(x = time_h, y = prop.dry, group = interaction(habitat, label), color = habitat)) +
  geom_vline(xintercept = c(0, 24, 48, 72, 96), color = "lightblue", size = 3, alpha = 0.3) +
  geom_vline(xintercept = c(2, 26, 50, 74, 98), color = "lightcoral", size = 1, alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Spatial sorting of ecotypes in response to flooding",
    x = "Time (hours)",
    y = "Proportion of individuals on dry side",
    color = "Habitat",
    linetype = "Habitat"
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = 14,  # Adjust x-position as needed
    y = 1,   # Adjust y-position as needed
    label = paste("p-value =", signif(p_value, 3)),
    color = "black",
    size = 3
  )


# Create a table of counts for each habitat and each condition (dry vs. wet)
table_data <- table(df$habitat, df$prop.dry > 0.5)
# Proportions test
prop_test <- prop.test(table_data)
print(prop_test)
# statistically significant difference in the proportions between the two groups (p-value = 0.02248)
