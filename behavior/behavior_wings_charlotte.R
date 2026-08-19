setwd("C:/Users/lotje/Documents/OneDrive - KU Leuven/Charlotte/Articles postdoc/1. Pogonus physio")
data <- read.csv(file = "phenotypes_belgium_june_2026.csv", sep = ";")
str(data)
data$population <- factor(data$population, c("Nieuwpoort", "Dudzele"))
data$species <- factor(data$species, c("Pchalceus"))
data$ecotype <- factor(data$ecotype, c("SW", "LW"))
data$sex <- factor(data$sex, c("M", "F"))

# View(data)

# Load libraries #####
library(car) # For Type III Anova
library(e1071) # For skewness
library(ggplot2) # For plotting
library(ggpubr) # For stats on the plots
library(survival) # For survival analysis
library(survminer) # For survival analysis
library(rptR)
library(tidyverse)
library(dplyr)
library(coxme)
library(lme4)
library(lmerTest)
library(tidyr)
#####

# Calculate %MRWS #####

# calculate %MRWS (Maximum Realizable Wing Size) 
# group 2: + and - functional flight muscles
# A_females = 1.2146 
# A_males = 1.1924 
# B_females = 0.8497 
# B_males = 0.8516 
# group 5: standard group (+ functional flight muscles and regularly flying)
A_females = 1.4654 
A_males = 1.4879 
B_females = 0.8454 
B_males = 0.8385 
data$Esize <- data$EL * data$EW
data$Wsize <- data$WL * data$WW

data$maxWsize <- NA

  for (i in 1:nrow(data)) {
  ifelse(data$sex[i]=='F',data$maxWsize[i]<-exp(A_females) * data$Esize[i]^B_females,
         ifelse(data$sex[i]=='M',data$maxWsize[i]<-exp(A_males) * data$Esize[i]^B_males,NA))
}

data$relMRWS <- data$Wsize/data$maxWsize

# sex ratio
m <- count(subset(data, data$sex == "M" & relMRWS>0)) # males with %MRWS measurement
f <- count(subset(data, data$sex == "F" & relMRWS>0)) # females with %MRWS measurement
m   # 51
f   # 57
m/f # 0.8947368
#####

# Models ####
mod_ws <- lm(relMRWS ~ population * sex + mass, data = data)
qqnorm(resid(mod_ws)); qqline(resid(mod_ws))
plot(fitted(mod_ws), resid(mod_ws)); abline(h = 0, col = "red")
hist(data$relMRWS)
shapiro.test(resid(mod_ws)) # W > 0.9, P > 0.01
leveneTest(relMRWS ~ population * sex, data = data) # Sign --> white.adjust = TRUE
Anova(mod_ws, type = "III", white.adjust = TRUE) # white.adjust = TRUE --> no requirement of variances to be constant (when leveneTest not oke)
summary(mod_ws)
emmeans(mod_ws, pairwise ~ sex)

library(emmeans)
SEplot_ws <- summary(emmeans(mod_ws, ~ population * sex)); SEplot_ws
data$pop_sex <- factor(interaction(data$population, data$sex),
                       levels = c("Nieuwpoort.M", "Nieuwpoort.F",
                                  "Dudzele.M", "Dudzele.F"))
ggplot(data, aes(x = pop_sex, y = relMRWS, fill = population)) +
  geom_boxplot(alpha = 0.4, width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(colour = population), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Nieuwpoort" = "dodgerblue1", "Dudzele" = "deeppink1")) +
  scale_colour_manual(values = c("Nieuwpoort" = "dodgerblue4", "Dudzele" = "deeppink4")) +
  xlab("Population and sex") +
  ylab("Maximal realizable wing size (%MRWS)") +
  scale_x_discrete(labels = c("Nieuwpoort.M" = "Nieuwpoort Male", "Nieuwpoort.F" = "Nieuwpoort Female", "Dudzele.M" = "Dudzele Male", "Dudzele.F" = "Dudzele Female")) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

SEplot_ws <- summary(emmeans(mod_ws, ~ population)); SEplot_ws
ggplot(data, aes(x = population, y = relMRWS, fill = population)) +
  geom_boxplot(alpha = 0.4, width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(colour = population), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Nieuwpoort" = "dodgerblue1", "Dudzele" = "deeppink1")) +
  scale_colour_manual(values = c("Nieuwpoort" = "dodgerblue4", "Dudzele" = "deeppink4")) +
  xlab("Population") +
  ylab("Maximal realizable wing size (%MRWS)") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

ggplot(data, aes(x = relMRWS, fill = population, colour = population)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_fill_manual(values = c(
    "Nieuwpoort" = "dodgerblue1",
    "Dudzele" = "deeppink1"
  )) +
  scale_colour_manual(values = c(
    "Nieuwpoort" = "dodgerblue4",
    "Dudzele" = "deeppink4"
  )) +
  labs(
    x = "Maximal realizable wing size (%MRWS)",
    y = "Density"
  ) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

SEplot_ws <- summary(emmeans(mod_ws, ~ sex)); SEplot_ws
ggplot(data, aes(x = sex, y = relMRWS, fill = sex)) +
  geom_boxplot(alpha = 0.4, width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(colour = sex), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("M" = "dodgerblue1", "F" = "deeppink1")) +
  scale_colour_manual(values = c("M" = "dodgerblue4", "F" = "deeppink4")) +
  xlab("Sex") +
  ylab("Maximal realizable wing size (%MRWS)") +
  scale_x_discrete(labels = c("M" = "Male", "F" = "Female")) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 20),
    axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2, "pt"), colour = "black"),
    axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2, "pt"), colour = "black"),
    axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
    axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
    panel.grid.major = element_line(colour = "grey90", size = 1),
    panel.grid.minor = element_line(colour = "grey95", size = 0.25),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

library(ggeffects)
pred_mass <- ggpredict(mod_ws, terms = "mass")

ggplot(pred_mass, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey80", alpha = 0.4) +
  geom_line(size = 1.2, colour = "black") +
  xlab("Mass") +
  ylab("Maximal realizable wing size (%MRWS)") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

# long data format
data <- data %>%
  mutate(across(c(t1, t2, t3), ~replace(., . == -30, 0)))
   
# TIME
time_long <- data %>%
  pivot_longer(
    cols = c(t1, t2, t3),
    names_to = "trial",
    values_to = "time"
  )

# mod_time <- lmer(time ~ population * sex + mass + trial + (1|labelSeq), data = time_long)
# ranova(mod_time) 
# qqnorm(resid(mod_time)); qqline(resid(mod_time))
# plot(fitted(mod_time), resid(mod_time)); abline(h = 0, col = "red")
# hist(time_long$time)
# shapiro.test(resid(mod_time)) # W > 0.9, P > 0.01
# leveneTest(time ~ population * sex, data = time_long) # Sign --> white.adjust = TRUE
# Anova(mod_time, type = "III", white.adjust = TRUE) # white.adjust = TRUE --> no requirement of variances to be constant (when leveneTest not oke)
# summary(mod_time)

# cox survival model
library(survival)
time_long$status <- ifelse(time_long$time == 300, 0, 1)
cox_model <- coxme(Surv(time, status) ~ population * sex + mass + trial + (1|labelSeq), data = time_long)
# cox_model <- coxme(Surv(time, status) ~ population + trial + (1|labelSeq), data = time_long)
summary(cox_model)
emmeans(cox_model, pairwise ~ population * sex)
# proportional hazards assumption
cox_fixed <- coxph(Surv(time, status) ~ population * sex + mass + trial, data = time_long)
cox.zph(cox_fixed) # p > 0.05 → proportional hazards assumption OK; p < 0.05 → variable changes effect over time
plot(cox.zph(cox_fixed))
# influential observations
library(survminer)
cox_diag <- coxph(Surv(time,status) ~ population * sex + mass + trial, data = time_long)
ggcoxdiagnostics(cox_diag, type="dfbeta")
# plot (kaplan-meier)
km <- survfit(Surv(time, status) ~ population + sex, data = time_long)
ggsurvplot(
  km, data = time_long, risk.table = F, 
  conf.int = TRUE, conf.int.alpha = 0.1, 
  fun = "event", # switched y-axis
  ylab = "Probability to emerge (%)",
  palette = c("dodgerblue1", "dodgerblue4", "deeppink1", "deeppink4"),
  size = 1.3)
km <- survfit(Surv(time, status) ~ population, data = time_long)
ggsurvplot(
  km, data = time_long, risk.table = F, 
  conf.int = TRUE, conf.int.alpha = 0.1, 
  fun = "event", # switched y-axis
  ylab = "Probability to emerge (%)",
  palette = c("dodgerblue4", "deeppink4"),
  size = 1.3)
# plot (kaplan-meier) per trial
km_trial <- survfit(Surv(time, status) ~ population + trial, data = time_long)
ggsurvplot(
  km_trial, data = time_long, risk.table = F,
  conf.int = TRUE, conf.int.alpha = 0.1,
  fun = "event",
  ylab = "Probability to emerge (%)",
  palette = c("lightblue1", "dodgerblue1", "dodgerblue4", "hotpink1", "deeppink1","deeppink4"),
  size = 1.3)

# relMRWS --> divide in groups
time_long$wing_group <- cut(time_long$relMRWS, breaks = quantile(time_long$relMRWS, probs = c(0, 1/3, 2/3, 1)),
  include.lowest = TRUE, labels = c("Low", "Medium", "High"))
time_long$wing_group <- cut(time_long$relMRWS, breaks = quantile(time_long$relMRWS, probs = seq(0, 1, by = 0.25), na.rm = TRUE), 
                            include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
km <- survfit(Surv(time, status) ~ wing_group, data = time_long)
ggsurvplot(
  km, data = time_long, risk.table = F, 
  conf.int = TRUE, conf.int.alpha = 0.1, 
  fun = "event", # switched y-axis
  ylab = "Probability to emerge (%)",
  palette = c("dodgerblue1", "dodgerblue4", "deeppink1", "deeppink4"),
  size = 1.3)

# Time per trial
time_long <- data %>%
  pivot_longer(
    cols = c(t1, t2, t3),
    names_to = "trial",
    values_to = "time"
  ) %>%
  mutate(
    trial = factor(
      trial,
      levels = c("t1", "t2", "t3"),
      labels = c("Trial 1", "Trial 2", "Trial 3")
    )
  )
plot_time <- ggplot(time_long, aes(x = trial, y = time, fill = trial)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.15, color = "black", alpha = 0.7) +
  facet_wrap(~population) +
  scale_x_discrete(labels = c("t1" = "Trial 1", "t2" = "Trial 2", "t3" = "Trial 3")) +
  scale_fill_manual(values = c("Trial 1" = "forestgreen", "Trial 2" = "steelblue", "Trial 3" = "firebrick")) +
  labs(x = "Trial", y = "Time until emerging (s)") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))
plot_time

# POSITION
pos_long <- data %>%
  pivot_longer(
    cols = c(p1, p2, p3),
    names_to = "trial",
    values_to = "position"
  ) %>%
  mutate(
    emerged = ifelse(position == "E",1,0)
  )

mod_pos <- glmer(emerged ~ population * sex + mass + trial + (1|labelSeq), family = binomial, data = pos_long)
library(DHARMa)
res_pos <- simulateResiduals(mod_pos); plot(res_pos) # check overdispersion
testDispersion(res_pos) # check overdispersion
testResiduals(res_pos) # check residuals
VarCorr(mod_pos) # check random effect
Anova(mod_pos, type = "III", white.adjust = TRUE) # white.adjust = TRUE --> no requirement of variances to be constant (when leveneTest not oke)
summary(mod_pos)

library(ggeffects)
pred <- ggpredict(mod_pos, terms = c("population", "sex"))
plot(pred)

SEplot_pos <- summary(emmeans(mod_pos, ~ population * sex, type = "response")); SEplot_pos
SEplot_pos$pop_sex <- factor(interaction(SEplot_pos$population, SEplot_pos$sex),
                       levels = c("Nieuwpoort.M", "Nieuwpoort.F",
                                  "Dudzele.M", "Dudzele.F"))
ggplot(SEplot_pos, aes(x = pop_sex, y = prob*100, colour = population,)) +
  geom_errorbar(aes(ymin = asymp.LCL*100, ymax = asymp.UCL*100), width = 0.2, size = 1, position = position_dodge(0.4)) +
  geom_point(size = 5, position = position_dodge(0.4)) +
  scale_colour_manual(values = c("Nieuwpoort" = "dodgerblue4", "Dudzele" = "deeppink4")) +
  # scale_shape_manual(values = c("M" = 16, "F" = 17)) +
  xlab("Population and sex") +
  ylab("Probability to emerge (%)") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  scale_x_discrete(labels = c("Nieuwpoort.M" = "Nieuwpoort Male", "Nieuwpoort.F" = "Nieuwpoort Female", "Dudzele.M" = "Dudzele Male", "Dudzele.F" = "Dudzele Female")) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

SEplot_pos <- summary(emmeans(mod_pos, ~ population, type = "response")); SEplot_pos
ggplot(SEplot_pos, aes(x = population, y = prob*100, colour = population,)) +
  geom_errorbar(aes(ymin = asymp.LCL*100, ymax = asymp.UCL*100), width = 0.2, size = 1, position = position_dodge(0.4)) +
  geom_point(size = 5, position = position_dodge(0.4)) +
  scale_colour_manual(values = c("Nieuwpoort" = "dodgerblue4", "Dudzele" = "deeppink4")) +
  # scale_shape_manual(values = c("M" = 16, "F" = 17)) +
  xlab("Population") +
  ylab("Probability to emerge (%)") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))


library(ggeffects)
pred_mass <- ggpredict(mod_pos, terms = "mass")
plot(pred_mass)

ggplot(pred_mass, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey80", alpha = 0.4) +
  geom_line(size = 1.2, colour = "black") +
  xlab("Mass (mg)") +
  ylab("Probability to emerge (%)") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

# plot barplot per trial
pos_summary <- data %>%
  pivot_longer(cols = c(p1, p2, p3), names_to = "trial", values_to = "position") %>%
  group_by(population, trial, position) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(population, trial) %>%
  mutate(prop = n / sum(n))
ggplot(pos_summary, aes(x = trial, y = prop, fill = position)) +
  geom_col(color = "black") +
  facet_wrap(~population) +
  scale_fill_manual(values = c("E" = "forestgreen", "S" = "steelblue"), labels = c("E" = "Emerged", "S" = "Submerged")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Trial", y = "Proportion of beetles", fill = "") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

# immediate emergence
time_long$immediate <- ifelse(time_long$time == 0, 1, 0)
table(time_long$immediate)
mod_immediate <- glmer(immediate ~ population * sex + mass + trial + (1|labelSeq), family = binomial, data = time_long)
library(DHARMa)
res_immediate <- simulateResiduals(mod_immediate); plot(res_immediate)
testDispersion(res_immediate)
testResiduals(res_immediate)
VarCorr(mod_immediate)
summary(mod_immediate)
exp(fixef(mod_immediate))
SEplot_immediate <- summary(emmeans(mod_immediate, ~ population * sex, type = "response"))
SEplot_immediate$pop_sex <- factor(interaction(SEplot_immediate$population, SEplot_immediate$sex),
                             levels = c("Nieuwpoort.M", "Nieuwpoort.F",
                                        "Dudzele.M", "Dudzele.F"))
ggplot(SEplot_immediate, aes(x = pop_sex, y = prob*100, colour = population,)) +
  geom_errorbar(aes(ymin = asymp.LCL*100, ymax = asymp.UCL*100), width = 0.2, size = 1, position = position_dodge(0.4)) +
  geom_point(size = 5, position = position_dodge(0.4)) +
  scale_colour_manual(values = c("Nieuwpoort" = "dodgerblue4", "Dudzele" = "deeppink4")) +
  # scale_shape_manual(values = c("M" = 16, "F" = 17)) +
  xlab("Population and sex") +
  ylab("Probability of immediate emergence (%)") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  scale_x_discrete(labels = c("Nieuwpoort.M" = "Nieuwpoort Male", "Nieuwpoort.F" = "Nieuwpoort Female", "Dudzele.M" = "Dudzele Male", "Dudzele.F" = "Dudzele Female")) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

# plot per trial
SEplot_immediate_trial <- summary(emmeans(mod_immediate, ~ trial, type = "response")); SEplot_immediate_trial
ggplot(SEplot_immediate_trial, aes(x = trial, y = prob * 100)) +
  geom_errorbar(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100), width = 0.15, size = 1) +
  geom_point(size = 5) +
  xlab("Trial") +
  ylab("Probability of immediate emergence (%)") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

#####



# Physiology ####

table <- xtabs(~population + sex, data); table

# CEA

ggplot(Data, aes(x = population, y = CEA, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$CEA)
fit <- lm(CEA ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(CEA ~ population * sex, data = Data)

# Anova(fit, type = "III")
Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "CEA.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("CEA") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# Ea

ggplot(Data, aes(x = population, y = Ea, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$Ea)
fit <- lm(Ea ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(Ea ~ population * sex, data = Data)

Anova(fit, type = "III")
# Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "Ea.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("Ea") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# EcETS

ggplot(Data, aes(x = population, y = EcETS, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$EcETS)
fit <- lm(EcETS ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(EcETS ~ population * sex, data = Data)

Anova(fit, type = "III")
# Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "EcETS.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("EcETS") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# Protein

ggplot(Data, aes(x = population, y = Protein, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$Protein)
fit <- lm(Protein ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(Protein ~ population * sex, data = Data)

# Anova(fit, type = "III")
Anova(fit, type = "II")
summary(fit)

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "Protein.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("Protein") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# Fat

ggplot(Data, aes(x = population, y = Fat, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$Fat)
fit <- lm(Fat ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(Fat ~ population * sex, data = Data)

Anova(fit, type = "III")
# Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "Fat.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("Fat") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# Sugar

ggplot(Data, aes(x = population, y = Sugar, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$Sugar)
fit <- lm(Sugar ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(Sugar ~ population * sex, data = Data)

# Anova(fit, type = "III")
Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "Sugar.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("Sugar") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# ETS

ggplot(Data, aes(x = population, y = ETS, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$ETS)
fit <- lm(ETS ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(ETS ~ population * sex, data = Data)

Anova(fit, type = "III")
# Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "ETS.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("ETS") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



# PO

ggplot(Data, aes(x = population, y = PO, fill = sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), color = "black", size = 2) +
  theme_bw()
plot(Data$PO)
fit <- lm(PO ~ population * sex, data = Data)
qqnorm(resid(fit)); qqline(resid(fit))
shapiro.test(resid(fit))
leveneTest(PO ~ population * sex, data = Data)

Anova(fit, type = "III")
# Anova(fit, type = "II")
summary(fit)
emmeans(fit, pairwise ~ population * sex, adjust = "fdr")

SEfit <- summary(emmeans(fit, ~ population * sex, type = "response")); SEfit
# svg(filename = "PO.svg", width = 8, height = 6)
ggplot(SEfit, aes(x = population, y = emmean, group = sex)) +
  scale_shape_manual(values = c(21, 21), labels = c("F", "M")) +
  scale_fill_manual(values = c("white", "black"), labels = c("F", "M")) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.2, position = position_dodge(.5)) +
  geom_point(aes(shape = sex, fill = sex), size = 5, position = position_dodge(.5)) +
  xlab("Ecotype") +
  ylab("PO") +
  # scale_x_discrete(breaks = c("E", "P"), labels = c(expression(italic("I. elegans")), expression(italic("I. pumilio")))) +
  # scale_y_continuous(limits = (c(0.013, 0.017)), breaks = c(0.013, 0.014, 0.015, 0.016, 0.017)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")), text = element_text(size = 20), axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"), axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"), axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")), axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")), panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill = NA, colour = NA), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", size = 20))
# dev.off()



#####



# Repeatability #####

# BEHAVIOURAL REPEATABILITY - 3 TRIALS

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# HOW OFTEN EACH BEETLE EMERGED (0-3 TRIALS)

data$emergeCount <- rowSums(
  data[c("p1", "p2", "p3")] == "E",
  na.rm = TRUE
)

# Frequency distribution
emergeDistribution <- data %>%
  count(
    emergeCount = factor(
      emergeCount,
      levels = 0:3
    )
  ) %>%
  mutate(
    proportion = n / sum(n),
    se = sqrt((proportion * (1 - proportion)) / sum(n)),
    ymin = pmax(0, proportion - se),
    ymax = pmin(1, proportion + se),
    emergeLabel = paste0(emergeCount, "/3")
  )


# PLOT ACTUAL EMERGENCE FREQUENCY

plot_emerge_frequency <- ggplot(
  emergeDistribution,
  aes(x = emergeLabel, y = proportion)
) +
  geom_bar(
    stat = "identity",
    width = 0.6
  ) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Emergence frequency per individual",
    y = "Beetles (%)"
  ) +
  theme_minimal(base_size = 18)

plot_emerge_frequency


# RANDOMIZATION TEST
# Randomize E/S within each trial while preserving the
# observed number of E and S individuals in each trial.
# Performed separately for each population.

set.seed(123)

randomize_population <- function(dat, nperm = 1000) {
  
  # observed emergence frequency
  observed_count <- rowSums(
    dat[c("p1", "p2", "p3")] == "E",
    na.rm = TRUE
  )
  
  observed_variance <- var(observed_count)
  
  # store randomized emergence frequencies
  random_counts <- matrix(
    NA,
    nrow = nrow(dat),
    ncol = nperm
  )
  
  for (i in 1:nperm) {
    
    s1 <- dat$p1
    s2 <- dat$p2
    s3 <- dat$p3
    
    s1[!is.na(s1)] <- sample(
      na.omit(s1),
      replace = FALSE
    )
    
    s2[!is.na(s2)] <- sample(
      na.omit(s2),
      replace = FALSE
    )
    
    s3[!is.na(s3)] <- sample(
      na.omit(s3),
      replace = FALSE
    )
    
    randomized <- data.frame(
      s1 = s1,
      s2 = s2,
      s3 = s3
    )
    
    random_counts[, i] <- rowSums(
      randomized == "E",
      na.rm = TRUE
    )
  }
  
  # variance for each permutation
  null_variance <- apply(
    random_counts,
    2,
    var
  )
  
  # one-sided permutation P-value
  p_value <- (
    sum(null_variance >= observed_variance) + 1
  ) / (nperm + 1)
  
  # frequency distribution for randomized data
  random_distribution <- matrix(
    0,
    nrow = 4,
    ncol = nperm
  )
  
  rownames(random_distribution) <- 0:3
  
  for (i in 1:nperm) {
    
    random_distribution[, i] <-
      table(
        factor(
          random_counts[, i],
          levels = 0:3
        )
      )
  }
  
  random_distribution <- data.frame(
    emergeCount = 0:3,
    mean = rowMeans(random_distribution),
    se = apply(
      random_distribution,
      1,
      function(x) sd(x) / sqrt(length(x))
    )
  )
  
  random_distribution <- random_distribution %>%
    mutate(
      mean_prop = mean / nrow(dat),
      ymin = pmax(
        0,
        mean_prop - se / nrow(dat)
      ),
      ymax = pmin(
        1,
        mean_prop + se / nrow(dat)
      ),
      emergeLabel = paste0(emergeCount, "/3")
    )
  
  list(
    observed_variance = observed_variance,
    null_variance = null_variance,
    p_value = p_value,
    distribution = random_distribution
  )
}


# RUN RANDOMIZATION SEPARATELY FOR EACH POPULATION

Nieuwpoort <- data %>%
  filter(population == "Nieuwpoort")

Dudzele <- data %>%
  filter(population == "Dudzele")


random_Nieuwpoort <- randomize_population(
  Nieuwpoort,
  nperm = 1000
)

random_Dudzele <- randomize_population(
  Dudzele,
  nperm = 1000
)


# OBSERVED VARIANCE AND PERMUTATION P-VALUES

random_Nieuwpoort$observed_variance
random_Nieuwpoort$p_value

random_Dudzele$observed_variance
random_Dudzele$p_value


# ACTUAL EMERGENCE DISTRIBUTION SEPARATELY BY POPULATION

actual_Nieuwpoort <- Nieuwpoort %>%
  count(
    emergeCount = factor(
      rowSums(
        cbind(p1, p2, p3) == "E",
        na.rm = TRUE
      ),
      levels = 0:3
    )
  ) %>%
  mutate(
    proportion = n / sum(n),
    se = sqrt(
      (proportion * (1 - proportion)) / sum(n)
    ),
    ymin = pmax(0, proportion - se),
    ymax = pmin(1, proportion + se),
    label = paste0(emergeCount, "/3"),
    type = "Actual",
    population = "Nieuwpoort"
  )


actual_Dudzele <- Dudzele %>%
  count(
    emergeCount = factor(
      rowSums(
        cbind(p1, p2, p3) == "E",
        na.rm = TRUE
      ),
      levels = 0:3
    )
  ) %>%
  mutate(
    proportion = n / sum(n),
    se = sqrt(
      (proportion * (1 - proportion)) / sum(n)
    ),
    ymin = pmax(0, proportion - se),
    ymax = pmin(1, proportion + se),
    label = paste0(emergeCount, "/3"),
    type = "Actual",
    population = "Dudzele"
  )


# RANDOMIZED DISTRIBUTIONS

random_plot_Nieuwpoort <-
  random_Nieuwpoort$distribution %>%
  transmute(
    label = emergeLabel,
    percent = mean_prop * 100,
    ymin = ymin * 100,
    ymax = ymax * 100,
    type = "Randomized",
    population = "Nieuwpoort"
  )


random_plot_Dudzele <-
  random_Dudzele$distribution %>%
  transmute(
    label = emergeLabel,
    percent = mean_prop * 100,
    ymin = ymin * 100,
    ymax = ymax * 100,
    type = "Randomized",
    population = "Dudzele"
  )


# COMBINE ACTUAL + RANDOMIZED DATA

plot_data_Nieuwpoort <- bind_rows(
  
  actual_Nieuwpoort %>%
    transmute(
      label,
      percent = proportion * 100,
      ymin = ymin * 100,
      ymax = ymax * 100,
      type,
      population
    ),
  
  random_plot_Nieuwpoort
)


plot_data_Dudzele <- bind_rows(
  
  actual_Dudzele %>%
    transmute(
      label,
      percent = proportion * 100,
      ymin = ymin * 100,
      ymax = ymax * 100,
      type,
      population
    ),
  
  random_plot_Dudzele
)


# PLOT: NIEUWPOORT

plot_action_randomized_Nieuwpoort <-
  ggplot(
    plot_data_Nieuwpoort,
    aes(
      x = label,
      y = percent,
      fill = type
    )
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    x = "Emergence frequency per individual",
    y = "Proportion of beetles (%)",
    fill = ""
  ) +
  # scale_fill_manual(
  #   values = c(
  #     "Actual" = "darkseagreen4",
  #     "Randomized" = "darkseagreen2"
  #   )
  # ) +
  scale_fill_manual(
    values = c(
      "Actual" = "dodgerblue1",
      "Randomized" = "dodgerblue4"
    )
  ) +
  ggtitle("Nieuwpoort") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

plot_action_randomized_Nieuwpoort


# PLOT: DUDZELE

plot_action_randomized_Dudzele <-
  ggplot(
    plot_data_Dudzele,
    aes(
      x = label,
      y = percent,
      fill = type
    )
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    x = "Emergence frequency per individual",
    y = "Proportion of beetles (%)",
    fill = ""
  ) +
  scale_fill_manual(
    values = c(
      "Actual" = "deeppink1",
      "Randomized" = "deeppink4"
    )
  ) +
  ggtitle("Dudzele") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.text = element_text(size = 20, margin = margin(2, 2, 10, 2,"pt")),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
        axis.text.y = element_text(size = 20, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
        axis.title.x = element_text(margin = margin(25, 2, 2, 2, "pt")),
        axis.title.y = element_text(margin = margin(2, 30, 2, 2, "pt")),
        panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_line(colour = "grey95", size = 0.25),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black", size = 20))

plot_action_randomized_Dudzele


# NULL DISTRIBUTION OF VARIANCE

hist(
  random_Nieuwpoort$null_variance,
  breaks = 30,
  main = "Nieuwpoort",
  xlab = "Variance in emergence frequency"
)

abline(
  v = random_Nieuwpoort$observed_variance,
  lwd = 3
)


hist(
  random_Dudzele$null_variance,
  breaks = 30,
  main = "Dudzele",
  xlab = "Variance in emergence frequency"
)

abline(
  v = random_Dudzele$observed_variance,
  lwd = 3
)

#####

# COMPARE BEHAVIOURAL REPEATABILITY BETWEEN POPULATIONS #####
# Fleiss' kappa + permutation test

library(dplyr)

# Function to calculate Fleiss' kappa

fleiss_kappa <- function(x) {
  
  # remove incomplete individuals
  x <- x[complete.cases(x), , drop = FALSE]
  
  # number of individuals and trials
  N <- nrow(x)
  n <- ncol(x)
  
  # categories
  categories <- sort(unique(as.vector(as.matrix(x))))
  
  # number of individuals in each category per individual
  counts <- matrix(
    0,
    nrow = N,
    ncol = length(categories)
  )
  
  colnames(counts) <- categories
  
  for (i in 1:N) {
    counts[i, ] <- table(
      factor(x[i, ], levels = categories)
    )
  }
  
  # agreement for each individual
  P_i <- (
    rowSums(counts^2) - n
  ) / (n * (n - 1))
  
  # mean observed agreement
  P_bar <- mean(P_i)
  
  # overall proportion assigned to each category
  p_j <- colSums(counts) / (N * n)
  
  # expected agreement
  P_e <- sum(p_j^2)
  
  # Fleiss' kappa
  kappa <- (P_bar - P_e) / (1 - P_e)
  
  return(kappa)
}


# Prepare data

repeatability_data <- data %>%
  select(labelSeq, population, p1, p2, p3) %>%
  filter(
    !is.na(p1),
    !is.na(p2),
    !is.na(p3)
  )


# Calculate Fleiss' kappa for each population

kappa_Nieuwpoort <- fleiss_kappa(
  repeatability_data %>%
    filter(population == "Nieuwpoort") %>%
    select(p1, p2, p3)
)

kappa_Dudzele <- fleiss_kappa(
  repeatability_data %>%
    filter(population == "Dudzele") %>%
    select(p1, p2, p3)
)


# Print kappas
kappa_Nieuwpoort
kappa_Dudzele


# Observed difference in Fleiss' kappa

observed_difference <- kappa_Nieuwpoort - kappa_Dudzele

observed_difference


# PERMUTATION TEST

set.seed(123)

nperm <- 10000

# population labels
population_labels <- repeatability_data$population

# number of individuals in each population
n_Nieuwpoort <- sum(population_labels == "Nieuwpoort")
n_Dudzele <- sum(population_labels == "Dudzele")


# vector to store differences
kappa_difference_random <- numeric(nperm)


for (i in 1:nperm) {
  
  # randomly shuffle population labels
  shuffled_labels <- sample(population_labels)
  
  # randomly assigned Nieuwpoort individuals
  random_Nieuwpoort <- repeatability_data[
    shuffled_labels == "Nieuwpoort",
    c("p1", "p2", "p3")
  ]
  
  # randomly assigned Dudzele individuals
  random_Dudzele <- repeatability_data[
    shuffled_labels == "Dudzele",
    c("p1", "p2", "p3")
  ]
  
  # calculate kappa for each randomized group
  kappa_random_Nieuwpoort <- fleiss_kappa(
    random_Nieuwpoort
  )
  
  kappa_random_Dudzele <- fleiss_kappa(
    random_Dudzele
  )
  
  # difference between randomized kappas
  kappa_difference_random[i] <-
    kappa_random_Nieuwpoort -
    kappa_random_Dudzele
}


# Two-sided permutation P-value

p_value_kappa <- (
  sum(
    abs(kappa_difference_random) >=
      abs(observed_difference)
  ) + 1
) / (nperm + 1)

p_value_kappa



# 95% permutation confidence interval


CI_kappa <- quantile(
  kappa_difference_random,
  probs = c(0.025, 0.975)
)

CI_kappa



# Final results


cat(
  "\nFleiss' kappa:\n",
  "Nieuwpoort =", round(kappa_Nieuwpoort, 3), "\n",
  "Dudzele    =", round(kappa_Dudzele, 3), "\n\n",
  "Difference (Nieuwpoort - Dudzele) =",
  round(observed_difference, 3), "\n",
  "95% permutation CI =",
  round(CI_kappa[1], 3), "to",
  round(CI_kappa[2], 3), "\n",
  "Permutation P =",
  round(p_value_kappa, 4), "\n"
)


# VISUALIZE THE NULL DISTRIBUTION

hist(
  kappa_difference_random,
  breaks = 50,
  main = "Permutation distribution of difference in Fleiss' kappa",
  xlab = "Difference in Fleiss' kappa (Nieuwpoort - Dudzele)"
)

abline(
  v = observed_difference,
  lwd = 3
)

abline(
  v = -observed_difference,
  lwd = 3,
  lty = 2
)



#####


# Summary inundation behaviour ####

summary(data[c('EL','EW','WL','WW','relMRWS')])

y_min <- 0.1675 # Adjust y-axis = relMRWS min
y_max <- 0.9560 # Adjust y-axis = relMRWS max

y_min <- 0
y_max <- 1

plot_relMWRS <- ggplot(
  data,
  aes(x = population, y = relMRWS, fill = population)
) +
  geom_violin(color = "black", alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.12, color = "black", alpha = 0.8) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.8,
    label.y = y_max * 0.98
  ) +
  labs(x = "", y = "%MRWS") +
  coord_cartesian(ylim = c(y_min, y_max)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Nieuwpoort" = "lightblue",
    "Dudzele" = "pink"
  )) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none"
  )
plot_relMWRS


time_long <- data %>%
  pivot_longer(
    cols = c(t1, t2, t3),
    names_to = "trial",
    values_to = "time"
  )

ggplot(time_long,
       aes(population, time, fill = population)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  geom_boxplot(width = 0.15, alpha = 0.8) +
  facet_wrap(~trial) +
  scale_fill_manual(values = c(
    "Nieuwpoort" = "lightblue",
    "Dudzele" = "pink"
  )) +
  theme_classic()


mod_time <- lmer(
  time ~ population + trial + (1|labelSeq),
  data = time_long
)

summary(mod_time)
anova(mod_time)


pos_long <- data %>%
  pivot_longer(
    cols = c(p1,p2,p3),
    names_to = "trial",
    values_to = "position"
  ) %>%
  mutate(
    emerged = ifelse(position=="E",1,0)
  )

pos_summary <- pos_long %>%
  group_by(population, trial) %>%
  summarise(
    prop = mean(emerged),
    .groups = "drop"
  )

ggplot(pos_summary,
       aes(trial, prop,
           colour = population,
           group = population)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0,1)
  ) +
  scale_colour_manual(values = c(
    "Nieuwpoort" = "lightblue",
    "Dudzele" = "pink"
  )) +
  theme_classic()

mod_pos <- glmer(
  emerged ~ population + trial +
    (1|labelSeq),
  family = binomial,
  data = pos_long
)

summary(mod_pos)
anova(mod_pos)

data$mean_time <- rowMeans(
  data[,c("t1","t2","t3")]
)

ggplot(data,
       aes(population,
           mean_time,
           fill = population)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_fill_manual(values = c(
    "Nieuwpoort" = "lightblue",
    "Dudzele" = "pink"
  )) +
  theme_classic()


#### models
mod_time_sex <- lmer(
  time ~ population + sex + trial + (1|labelSeq),
  data = time_long
)
summary(mod_time_sex)
anova(mod_time_sex)

mod_time_interaction <- lmer(
  time ~ population * sex + trial + (1|labelSeq),
  data = time_long
)
summary(mod_time_interaction)
anova(mod_time_interaction)

mod_pos_sex <- glmer(
  emerged ~ population + sex + trial +
    (1|labelSeq),
  family = binomial,
  data = pos_long
)
summary(mod_pos_sex)
anova(mod_pos_sex)

mod_pos_interaction <- glmer(
  emerged ~ population * sex + trial +
    (1|labelSeq),
  family = binomial,
  data = pos_long
)
summary(mod_pos_interaction)
anova(mod_pos_interaction)

anova(mod_pos, mod_pos_sex, mod_pos_interaction)
anova(mod_time, mod_time_sex, mod_time_interaction)

summary(mod_pos_sex)
anova(mod_pos_sex)

summary(mod_time_sex)
anova(mod_time_sex)

# Nieuwpoort beetles took on average ~74 seconds longer to emerge than Dudzele beetles.
# Dudzele (LW) beetles emerged significantly faster than Nieuwpoort (SW) beetles.
# Nieuwpoort beetles had only ~29% of the odds of being emerged compared to Dudzele beetles.
# Dudzele beetles were about 3.4 times more likely to be emerged than Nieuwpoort beetles.

pos_summary <- data %>%
  pivot_longer(
    cols = c(p1, p2, p3),
    names_to = "trial",
    values_to = "position"
  ) %>%
  group_by(population, trial, position) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(population, trial) %>%
  mutate(prop = n / sum(n))

ggplot(pos_summary,
       aes(x = trial,
           y = prop,
           fill = position)) +
  geom_col(color = "black") +
  facet_wrap(~population) +
  scale_fill_manual(
    values = c(
      "E" = "forestgreen",
      "S" = "steelblue"
    ),
    labels = c(
      "E" = "Emerged",
      "S" = "Submerged"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Trial",
    y = "Proportion of beetles",
    fill = ""
  ) +
  theme_classic(base_size = 16)


pos_emerged <- data %>%
  pivot_longer(
    cols = c(p1,p2,p3),
    names_to = "trial",
    values_to = "position"
  ) %>%
  mutate(emerged = position == "E") %>%
  group_by(population, trial) %>%
  summarise(
    prop_emerged = mean(emerged),
    .groups = "drop"
  )

ggplot(pos_emerged,
       aes(trial,
           prop_emerged,
           color = population,
           group = population)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c(
    "Nieuwpoort" = "lightblue",
    "Dudzele" = "pink"
  )) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,1)
  ) +
  labs(
    x = "Trial",
    y = "Proportion emerged"
  ) +
  theme_classic(base_size = 16)



surv_data <- data %>%
  pivot_longer(
    cols = c(t1, t2, t3),
    names_to = "trial",
    values_to = "time"
  ) %>%
  mutate(
    time = ifelse(time < 0, 0, time),
    event = ifelse(time < 300, 1, 0)
  )


surv_obj <- Surv(
  time = surv_data$time,
  event = surv_data$event
)

km_fit <- survfit(
  surv_obj ~ population,
  data = surv_data
)

ggsurvplot(
  km_fit,
  data = surv_data,
  pval = TRUE,
  risk.table = FALSE,
  conf.int = TRUE,
  xlab = "Time to emergence (s)",
  ylab = "Probability of remaining submerged",
  legend.title = "",
  palette = c(
    "population=Nieuwpoort" = "lightblue",
    "population=Dudzele" = "pink"
  )
)


cox_model <- coxph(
  surv_obj ~ population,
  data = surv_data
)


summary(cox_model)


exp(coef(cox_model))
exp(confint(cox_model))


cox_model_repeated <- coxph(
  surv_obj ~ population +
    strata(trial) +
    cluster(labelSeq),
  data = surv_data
)
summary(cox_model_repeated)
# Dudzele beetles emerge about 2.1 times faster than Nieuwpoort beetles (estimated to emerge between 1.2 and 3.7 times faster than Nieuwpoort beetles).
# Long-winged Dudzele beetles respond to flooding more rapidly and are more likely to emerge than short-winged Nieuwpoort beetles.
# Nieuwpoort beetles showed a significantly lower emergence rate than Dudzele beetles (HR = 0.48, 95% CI = 0.27–0.83, p = 0.010).

#####


