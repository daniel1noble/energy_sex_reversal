##########################################
######### SMR lower 10%  O2 data #########
##########################################
# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "patchwork", "ggExtra", "gridExtra", "cowplot", "reshape2", "bayestestR")

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
         sex = Geno.pheno,
         mass_g = mass) %>% 
  mutate(ztime = scale(time),
         logMass = scale(log(mass_g), scale = FALSE),
         zstartmass = scale(log(start_mass_g), scale = FALSE),
         zendmass = scale(log(end_mass_g), scale = FALSE))%>% 
  dplyr::select(-X, -Date.Hatched, -MR_O2_min)

################################################
# Calculating SMR - keeping bottom 10% of data #
################################################
bassiana.data <- bassiana.data %>% 
  group_by(match) %>% # match accounts for individual and date
  arrange(O2_min) %>% 
  slice(1:floor(n()*.1))


#####################################
# Bassiana Models - Full Data
#####################################
#####
## residual function - used to check residuals for all models
#####
residuals_brms <- function(model, data){
  fitted <- predict(model)[1:dim(data)[1], 1] # Predcit mean for each data
  e <- with(data, log(O2_min)) - fitted
  return(e)
}
##############
## Model 1 -  name:Bas_m1_brms_smr
##############
Bas_m1_brms_smr <- brm(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4)
Bas_m1_brms_smr <- add_criterion(Bas_m1_brms_smr, c("loo"))
saveRDS(Bas_m1_brms_smr, "./models/SMR_mods/Bas_m1_brms_smr")
Bas_m1_brms_smr <- readRDS("./models/SMR_mods/Bas_m1_brms_smr")
####################
# Model 1 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
draws <- as.array(Bas_m1_brms_smr)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_logMass"), lags =10)
# check residuals
e <- residuals_brms(Bas_m1_brms_smr, bassiana.data)
hist(e)
# plots
plot(Bas_m1_brms_smr)
#R2 and summary of full model
bayes_R2(Bas_m1_brms_smr)
summary(Bas_m1_brms_smr)

##############
# Model 2 - Bas_m2_brms_smr (heteroskedastic model)
##############
mod_bas <- bf(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                sigma ~ logMass + ztime)
Bas_m2_brms_smr <- brm(mod_bas, family = "gaussian", 
                     data = bassiana.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
Bas_m2_brms_smr <- add_criterion(Bas_m2_brms_smr, c("loo"))
saveRDS(Bas_m2_brms_smr, "./models/SMR_mods/Bas_m2_brms_smr")
Bas_m2_brms_smr <- readRDS(file = "models/SMR_mods/Bas_m2_brms_smr")

####################
# Model 2 checks: lags, residuals, r2, summary
####################
draws <- as.array(Bas_m2_brms_smr)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_logMass"), lags =10)
# Check residuals
e <- residuals_brms(Bas_m2_brms_smr, bassiana.data)
hist(e)
# plots
plot(Bas_m2_brms_smr)
#R2 and summary
bayes_R2(Bas_m2_brms_smr)
summary(Bas_m2_brms_smr)

####################
# Model 1 & Model 2 comparison
####################
# SE is large and difference is not sufficient to warrant a heteroskedastic model 
loo_compare(Bas_m1_brms_smr, Bas_m2_brms_smr)

####################  
# extract posteriors for plotting and hypothesis testing for interaction (sex*mass) 
####################
summary(Bas_m1_brms_smr)
post_Bas_m1 <- posterior_samples(Bas_m1_brms_smr, pars = "^b")
variable.names(post_Bas_m1)

## extracting posteriors for interaction of sex and mass
XXf.mass.posterior <- as.array(post_Bas_m1[,"b_logMass"])
XXm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXXm:logMass"])
XYm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXYm:logMass"])

## pMCMC Function
# Calculates the, p-value or pMCMC value for a posterior distribution
# x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
pmcmc <- function(x){
  2*(1 - max(table(x<0) / nrow(x)))
}

# H1: Like Phenotype Hypothesis - pMCMC
Bass.RslopeDiff.Pheno.mass <- XYm.mass.posterior - XXm.mass.posterior
Bass.pMCMC_phenotype_metabolism <- pmcmc(Bass.RslopeDiff.Pheno.mass)
Bass.pMCMC_phenotype_metabolism

# H2: Like Genotype Hypothesis - pMCMC
Bass.RslopeDiff.Geno.mass <- XXf.mass.posterior - XYm.mass.posterior
Bass.pMCMC_genotype_metabolism <- pmcmc(Bass.RslopeDiff.Geno.mass)
Bass.pMCMC_genotype_metabolism

## plotting posteriors accounting for sex*mass interaction
Bass.mass.o2 <- cbind(XXf.mass.posterior, XXm.mass.posterior, XYm.mass.posterior)
mcmc_areas(Bass.mass.o2, 
           pars = c("XXf.mass.posterior", "XXm.mass.posterior", "XYm.mass.posterior"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope Differences") 



####################################
######### Pogona  O2 data  ######### 
####################################
#########
#Load data
########
pogona.data <- read.csv("./final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min,
         sex = Geno.pheno,
         mass_g = mass) %>% 
  group_by(sex) %>% 
  mutate(ztime = scale(time),
         logmass = log(mass_g), 
         logMass = scale(log(mass_g), scale = FALSE),
         zstartmass = scale(log(start_mass_g), scale = FALSE),
         zendmass = scale(log(end_mass_g), scale = FALSE)) %>% 
  dplyr::select(-X, -Date.Hatched, -MR_O2_min)

################################################
# Calculating SMR - keeping bottom 10% of data #
################################################
pogona.data <- pogona.data %>% 
  group_by(match) %>% # match accounts for individual and date
  arrange(O2_min) %>% 
  slice(1:floor(n()*.1))

##############
## Model 1 - name:Pog_m1_brms
##############
Pog_m1_brms_smr <- brm(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4)
Pog_m1_brms_smr <- add_criterion(Pog_m1_brms_smr, c("loo"))
saveRDS(Pog_m1_brms_smr, "./models/SMR_mods/Pog_m1_brms_smr")
Pog_m1_brms_smr <- readRDS("./models/SMR_mods/Pog_m1_brms_smr")
####################
# Model 1 checks: lags, residuals, r2, summary
####################
draws <- as.array(Pog_m1_brms_smr)
mcmc_acf(draws,  pars = c("b_Intercept", "b_sexZZf", "b_sexZZm"), lags =10)
# check residuals
e <- residuals_brms(Pog_m1_brms_smr, pogona.data)
hist(e)
# plots
plot(Pog_m1_brms_smr)
#R2 of and summary
bayes_R2(Pog_m1_brms_smr)
summary(Pog_m1_brms_smr)

##############
## Model 2 - Pog_m2_brms_smr (heteroskedastic model)
##############
mod <- bf(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
            sigma ~ logMass + ztime)
Pog_m2_brms_smr <- brm(mod, family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
Pog_m2_brms_smr <- add_criterion(Pog_m2_brms_smr, c("loo"))
saveRDS(Pog_m2_brms_smr, "./models/SMR_mods/Pog_m2_brms_smr")
Pog_m2_brms_smr <- readRDS(file = "./models/SMR_mods/Pog_m2_brms_smr")
####################
# Model 2 - checks: lags, residuals, r2, summary
####################
draws <- as.array(Pog_m2_brms_smr)
mcmc_acf(draws,  pars = c("b_Intercept", "b_sexZZf", "b_sexZZm"), lags =10)
# check residuals
e <- residuals_brms(Pog_m2_brms_smr, pogona.data)
hist(e)
# plots
plot(Pog_m2_brms_smr)
#R2 of and summary
bayes_R2(Pog_m2_brms_smr)
summary(Pog_m2_brms_smr)

####################  
# extract posteriors for plotting and hypothesis testing for interaction (sex*mass) 
####################
post_pog_m2 <- posterior_samples(Pog_m2_brms_smr, pars = "^b")
Pog_m2_brms_smr
dimnames(post_pog_m2)

## extracting posteriors for interaction of sex and mass
ZWf.mass.posterior <- as.array(post_pog_m2[,"b_logMass"])
ZZf.mass.posterior <- as.array(post_pog_m2[,"b_logMass"] + 
                                 post_pog_m2[,"b_sexZZf:logMass"])
ZZm.mass.posterior <- as.array(post_pog_m2[,"b_logMass"] + 
                                 post_pog_m2[,"b_sexZZm:logMass"])

## H1: Like Phenotype Hypothesis - pMCMC
Pog.RslopeDiff.Pheno.mass <- ZWf.mass.posterior - ZZf.mass.posterior
Pog.pMCMC_phenotype_metabolism <- pmcmc(Pog.RslopeDiff.Pheno.mass)
Pog.pMCMC_phenotype_metabolism

## H2: Like Genotype Hypothesis - pMCMC
Pog.RslopeDiff.Geno.mass <- ZZf.mass.posterior - ZZm.mass.posterior
Pog.pMCMC_genotype_metabolism <- pmcmc(Pog.RslopeDiff.Geno.mass)
Pog.pMCMC_genotype_metabolism

## plotting posteriors accounting for sex*mass interaction
pog.dat <- cbind(ZWf.mass.posterior, ZZf.mass.posterior, ZZm.mass.posterior)
mcmc_areas(pog.dat, 
           pars = c("ZWf.mass.posterior", "ZZf.mass.posterior", "ZZm.mass.posterior"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=10)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 

