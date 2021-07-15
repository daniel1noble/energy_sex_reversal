################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 
setwd("~/OneDrive - Australian National University/Respirometry")
# Packages
packages<-c("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp")
sapply(packages, require, character.only=T)

######### Bassiana  O2 data ######### 
# Load data and rename cols to make easy to follow
bassiana.data <- read.csv("./final.analysis.data/Bassiana.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min, 
         bmr_O2 = BMR_O2,
         sex = Geno.pheno,
         mass_g = mass) %>% 
  mutate(ztime = scale(time)) %>% 
  dplyr::select(-X, -Date.Hatched, -MR_O2_min)
# quick plot
# box/violin plot
	fig <- ggplot(bassiana.data, aes(x = sex, y = O2_min)) + 
	  geom_point() +
	  geom_violin() + 
	  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
	  geom_jitter(fill = "grey", alpha = 0.3) + 
	  labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "Sex")+ 
	  theme_bw()
	fig
# Model 1
# comparing o2 measurements across sex for bassiana
# accounting for random factor of lizard, sex day and time (maker_sample)
Bas_m1_lmer <- lmer(log(O2_min) ~ sex + scale(time) + (1 + ztime | id) + (1 | day), data = bassiana.data)
summary(Bas_m1_lmer)
anova(Bas_m1_lmer )
AIC(Bas_m1_lmer )
# check final model
anova(Bas_m1_lmer)
hist(residuals(Bas_m1_lmer))
plot(residuals(Bas_m1_lmer))

######### full data
Bas_all_model1 <- lmer(log(O2_min) ~ sex*log(mass_g) + ztime + log(mass_g) + (1 + ztime | id) + (1 | day), data = bassiana.data)
summary(Bas_all_model1)
lsmeans(Bas_all_model1, pairwise ~ sex)
# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(bassiana.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log Metabolic Rate") +
  ggtitle("Bassiana MR") 

# BMR
Bas_bmr_1 <- lmer(log(bmr_O2) ~ sex + log(mass_g) + (1 | id) + (1 | day), data = bassiana.data)
summary(Bas_bmr_1)
confint(Bas_bmr_1)
# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(bassiana.data, aes(log(mass_g), log(bmr_O2), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log BMR") +
  ggtitle("Bassiana BMR") 

# For repeatability
rep_bas <- lmer(log(O2_min) ~ 1 + (1 | id) + (1 | day), data = bassiana.data)
rep_bas


######### Pogona  O2 data  ######### 
# Load data and rename cols to make easy to follow
pogona.data <- read.csv("./final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min, 
         bmr_O2 = BMR_O2,
         sex = Geno.pheno,
         mass_g = mass) %>% 
  mutate(ztime = scale(time)) %>% 
  dplyr::select(-X, -Date.Hatched, -MR_O2_min)
# quick plot
# box/violin plot
fig <- ggplot(pogona.data, aes(x = sex, y = O2_min)) + 
  geom_point() +
  geom_violin() + 
  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
  geom_jitter(fill = "grey", alpha = 0.3) + 
  labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "Sex")+ 
  theme_bw()
fig

# Model 1
# comparing o2 measurements across sex for Pogona
# accounting for random factor of lizard, sex day and time (maker_sample)
Pog_m1_lmer <- lmer(log(O2_min) ~ sex + scale(time) + (1 + ztime | id) + (1 | day), data = pogona.data)
summary(Pog_m1_lmer)
anova(Pog_m1_lmer )
AIC(Pog_m1_lmer )
# check final model
anova(Pog_m1_lmer)
hist(residuals(Pog_m1_lmer))
plot(residuals(Pog_m1_lmer))

######### full data
Pog_all_model1 <- lmer(log(O2_min) ~ sex*log(mass_g) + ztime + log(mass_g) + (1 + ztime | id) + (1 | day), data = pogona.data)
summary(Pog_all_model1)
lsmeans(Pog_all_model1, pairwise ~ sex)
# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(pogona.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log Metabolic Rate") +
  ggtitle("Pogona") 

# BMR
Pog_bmr_1 <- lmer(log(bmr_O2) ~ sex + log(mass_g) + (1 | id) + (1 | day), data = pogona.data)
summary(Pog_bmr_1)
confint(Pog_bmr_1)
# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(pogona.data, aes(log(mass_g), log(bmr_O2), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log BMR") +
  ggtitle("Pogona BMR") 

# For repeatability
rep_pogona <- lmer(log(O2_min) ~ 1 + (1 | id) + (1 | day), data = pogona.data)
rep_pogona



