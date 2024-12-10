# Load libraries
library(car) # For Type III Anova
library(dplyr) # For data manipulation
library(e1071) # For skewness
library(ggplot2) # For plotting
library(ggpubr) # For stats on the plots
library(survival) # For survival analysis
library(survminer) # For survival analysis

setwd("~/Documents/behavior-water-experiments")

data <- read.csv("inundation.csv")
data$escape <- ifelse(data$time_s == 1200, "no", "yes")
data$ecotype <- as.factor(data$ecotype)
data$region <- as.factor(data$label)
data$escape <- as.factor(data$escape) # Binomial response (0 = no, 1 = yes)

# Visualize the time data
hist(data$time_s, main = "Histogram of Time to Resurface", xlab = "Time (seconds)", breaks = 30)

# Remove the rows that cotain the NA
clean_data <- na.omit(data)

# Check for skewness
skewness(clean_data$time_s)
# Positive skew (> 0) means the data has a right tail, where the majority of the data points are concentrated on the left side
# Okay to proceed with the analysis using models like GLM, which can handle slightly skewed data

# Fit a Gamma GLM with a log link function
# I use a gamma GLM since I have continuous time data (binomial would make sense if i was looking at the binary escape or no escape)
# The log link function is used to ensure that the predicted values stay positive, which is important since the time cannot be negative
glm_model <- glm(time_s ~ ecotype, family = Gamma(link = "log"), data = data)

# Summary of the model
summary(glm_model) #best model, 2449

# Perform Type III Anova to test the significance of the factors
# LR, likelihood ratio
# looks at how well adding a specific factor improves the fit of the model
anova_glm <- Anova(glm_model, type = 3, test = "LR")
print(anova_glm)
#  ecotype significantly affects the time to resurface -> two ecotypes (SW vs LW) differ in how long it takes them to resurface after being submerged.
# SW beetles took significantly longer to emerge compared to LW beetles

## TRY OUT OTHER MODELS!
# Fit a model with both ecotype and region as predictors
glm_model_region <- glm(time_s ~ ecotype + region, family = Gamma(link = "log"), data = data)
AIC(glm_model_region)

# Fit a model with an interaction between ecotype and region
glm_model_interaction <- glm(time_s ~ ecotype * region, family = Gamma(link = "log"), data = data)
AIC(glm_model_interaction)

# Fit a Gaussian GLM
glm_model_gaussian <- glm(time_s ~ ecotype, family = gaussian(link = "identity"), data = data)
AIC(glm_model_gaussian)

# None of these models are better than the first one (according to their AIC)

# Get the coefficient for 'ecotype' from the GLM model
# A coefficient represents the strength and direction of the effect of that the ecotype has on the outcome (time)
coef_estimate <- coef(glm_model)["ecotypeSW"]

# Exp the coefficient to get the relative effect size (transforms it back to the original scale, since we logged it before)
# the effect size tells me how much the ecotype impacts the response variable (time, in this case)
effect_size <- exp(coef_estimate)
print(effect_size) #2.05, greater than 1 means SW resurfaces later

# Comparing mean values directly, test significance in a simpler way
# Data is not normal, use non-parametric test
wilcox_test_result <- wilcox.test(time_s ~ ecotype, data = data)
print(wilcox_test_result)

### PLOTS!!! make them pretty, add stats when necessary (in the case of boxplots a simple non-par t test is enough)
# boxplot with jitter points
ggplot(clean_data, aes(x = ecotype, y = time_s, fill = ecotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot without outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") +  # Add jittered points
  stat_compare_means(label.y = max(clean_data$time_s) + 50) +  # Add significance
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Ecotype") +
  labs(title = "Emergence Time by Ecotype") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("lightcoral", "lightblue"))

# violin plot with boxplot on top of it
ggplot(clean_data, aes(x = ecotype, y = time_s, fill = ecotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot for density
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Overlay boxplot
  stat_compare_means(label.y = max(clean_data$time_s) + 650) +
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Ecotype") +
  labs(title = "Emergence Time by Ecotype") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("lightcoral", "lightblue"))  # Custom colors

# density plot (time is a continuous variable, doesn't need to be shown as steps)
ggplot(clean_data, aes(x = time_s, fill = ecotype)) +
  geom_density(alpha = 0.6) +  # Density plot
  scale_x_continuous(name = "Emergence Time (seconds)") +
  scale_y_continuous(name = "Density") +
  labs(title = "Density of Emergence Time by Ecotype") +
  theme_minimal() +
  scale_fill_manual(values = c("lightcoral", "lightblue")) +
  theme(legend.title = element_blank())

# combined (box + violin + jitter)
ggplot(clean_data, aes(x = ecotype, y = time_s, fill = ecotype)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.3, color = "black") +
  stat_compare_means(label.y = max(clean_data$time_s) + 650) +
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Ecotype") +
  labs(title = "Emergence Time by Ecotype") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("lightcoral", "lightblue"))


## SURVIVAL ANALYSIS !!!

# Survival analysis is done when the time to event is important, in this case, if beetles have a value of 1200, this means they did not emerge
# If some beetles did not resurface within 20 minutes, their times are "censored"
clean_data <- clean_data %>%
  mutate(event = ifelse(time_s < 1200, 1, 0))  # 1 = Resurfaced, 0 = Censored

# Create a survival object (to use in the model)
surv_obj <- Surv(time = clean_data$time_s, event = clean_data$event)

# Kaplan-Meier model (NON parametric), estimates the survival probability over time for each ecotype
# It provides an unbiased estimate of the proportion of individuals that remain at risk of resurfacing at each time point
km_fit <- survfit(surv_obj ~ ecotype, data = clean_data)

# Plot the survival curve (km plot)
ggsurvplot(km_fit, data = clean_data, pval = TRUE,
           xlab = "Time (seconds)", ylab = "Probability of Staying Underwater",
           title = "Kaplan-Meier Survival Curve",
           legend.title = "Ecotype", palette = c("lightcoral", "lightblue"))

# Cox proportional hazards model
# tests the significance of ecotype on the hazard rate (rate of emergence over time)
cox_model <- coxph(surv_obj ~ ecotype, data = clean_data)
summary(cox_model)

# calculate Hazard Ratio (HR) and confidence intervals for easier interpretation 
# Hazard is the immediate "risk" of coming out of the water for a beetle at a given time
# A higher hazard means a beetle is more likely to resurface quickly
exp(coef(cox_model))  # HR
# ecotypeSW, 0.4182992

# The hazard ratio for ecotypeSW is 0.4183. This value represents the likelihood of resurfacing for SW beetles compared to 
# LW beetles. Since the HR is less than 1, it indicates that SW beetles are less likely to resurface
# The hazard of resurfacing for SW beetles is 58% (1 - 0.4183 x 100) lower than that of LW beetles
