################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 

# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa")

#####################################
######### Bassiana  O2 data #########
#####################################
#########
#Load data
#########
  bassiana.data <- read.csv("./final.analysis.data/Bassiana.finalO2.sexreversal.analysis.data.clean.csv") %>% 
    rename(day =date.dd.mm.yy.,
           time = marker_sample,
           id = bd_liz_id, 
           O2_min = Final_MR_O2_min, 
           bmr_O2 = BMR_O2,
           sex = Geno.pheno,
           mass_g = mass) %>% 
    mutate(ztime = scale(time),
           zlogMass = scale(log(mass_g), scale = FALSE)) %>% 
    dplyr::select(-X, -Date.Hatched, -MR_O2_min)
  # quick plot
  # box/violin plot
  	fig <- ggplot(bassiana.data, aes(x = sex, y = log(O2_min))) + 
  	  geom_point() +
  	  geom_violin() + 
  	  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
  	  geom_jitter(fill = "grey", alpha = 0.3) + 
  	  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "Sex")+ 
  	  theme_bw()
  	fig

#####################################
# Bassiana Models - Full Data
#####################################
########
## Model 1
########
suppressMessages(
    Bas_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                   family = "gaussian", data = bassiana.data, iter= 2000, warmup = 1000, thin = 1, control = list(adapt_delta = 0.95)))
  	Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo", "waic"))
saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
summary(Bas_m1_brms)

# extract posteriors
post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")
b_XX_male_slope <- post_Bas_m1[,"b_zlogMass"] + post_Bas_m1[,"b_sexXX_Male:zlogMass"]
b_XY_male_slope <- post_Bas_m1[,"b_zlogMass"] + post_Bas_m1[,"b_sexXY_Male:zlogMass"]

# contarst slopes
contrastSlope <- as.mcmc(post_meanXX_male_slope - post_meanXY_male_slope)
mean(contrastSlope)
HPDinterval(contrastSlope)

#R2 of full model
bayes_R2(post_Bas_m1)

# Conditional effects
conditional_effects(post_Bas_m1)

# Model checks
plot(post_Bas_m1)

##############
# Model 2
##############
mod_bas <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
          sigma ~ zlogMass + ztime)
Bas_m2_brms <- brm(mod_bas, family = gaussian(), data = bassiana.data, iter= 2000, warmup = 1000, thin = 1)
Bas_m2_brms <- add_criterion(Bas_m2_brms, c("loo", "waic"))
saveRDS(Bas_m2_brms, "./models/Bas_m2_brms")
bayes_R2(Bas_m2_brms)

# Model checks
plot(Bas_m2_brms)
summary(Bas_m2_brms)


####################
# Model comparison
####################
loo_compare(Bas_m1_brms, Bas_m2_brms)


####################
# Plots
####################
# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(bassiana.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log Metabolic Rate") +
  ggtitle("Bassiana MR") 

####################################
######### Pogona  O2 data  ######### 
####################################

# Load data and rename cols to make easy to follow
pogona.data <- read.csv("./final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min, 
         bmr_O2 = BMR_O2,
         sex = Geno.pheno,
         mass_g = mass) %>% 
  mutate(ztime = scale(time),
         zlogMass = scale(log(mass_g), scale = FALSE)) %>% 
  dplyr::select(-X, -Date.Hatched, -MR_O2_min)

# quick plot
# box/violin plot
fig <- ggplot(pogona.data, aes(x = sex, y = log(O2_min))) + 
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
Pog_m1_brms <- brm(log(O2_min) ~ sex + ztime + (1 + ztime | id) + (1 | day),  
                   family = gaussian(), data = pogona.data, iter= 2000, warmup = 1000, thin = 1)
Pog_m1_brms <- add_criterion(Pog_m1_brms, c("loo", "waic"))
saveRDS(Pog_m1_brms, "./models/Pog_m1_brms")
summary(Pog_m1_brms)

# Model 2
Pog_m2_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                   family = gaussian(), data = pogona.data, iter= 2000, warmup = 1000, thin = 1)
Pog_m2_brms <- add_criterion(Pog_m2_brms, c("loo", "waic"))
saveRDS(Pog_m2_brms, "./models/Pog_m2_brms")
summary(Pog_m2_brms)
bayes_R2(Pog_m2_brms)

# Model 3
        mod <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                        sigma ~ zlogMass + ztime)
Pog_m3_brms <- brm(mod, family = gaussian(), data = pogona.data, iter= 2000, warmup = 1000, thin = 1)
Pog_m3_brms <- add_criterion(Pog_m3_brms, c("loo", "waic"))
saveRDS(Pog_m3_brms, "./models/Pog_m3_brms")
bayes_R2(Pog_m3_brms)

loo_compare(Pog_m2_brms,Pog_m3_brms)
## Posteriors
post_Pog_m2_brms <- posterior_samples(Pog_m2_brms, "^b")

# Slope contrast
b_sexZZf <- post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZf:zlogMass"]
b_sexZZm <- post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZm:zlogMass"]

contrast_ZZ <- as.mcmc(b_sexZZf - b_sexZZm)
mean(contrast_ZZ)
HPDinterval(contrast_ZZ)


# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(pogona.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log Metabolic Rate") +
  ggtitle("Pogona") 

