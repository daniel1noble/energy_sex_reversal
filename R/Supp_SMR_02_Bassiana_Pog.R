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
Bass.RslopeDiff.Geno.mass <- XXf.mass.posterior - XXm.mass.posterior
Bass.pMCMC_genotype_metabolism <- pmcmc(Bass.RslopeDiff.Geno.mass)
Bass.pMCMC_genotype_metabolism

## plotting posteriors accounting for sex*mass interaction
Bass.mass.o2 <- cbind(XXf.mass.posterior, XXm.mass.posterior, XYm.mass.posterior)
color_scheme_set("darkgray")
mcmc_areas(Bass.mass.o2, 
           pars = c("XXf.mass.posterior", "XXm.mass.posterior", "XYm.mass.posterior"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope Differences") 

####################  
# manual predict values for regression lines for Figure 2A
####################
# XX females
a_female <- post_Bas_m1$b_Intercept
b_mass_female <- post_Bas_m1$b_logMass
logMass = seq(-0.3856391, 0.3710889, 
              length.out = 100)
y_female <- a_female + b_mass_female %*% t(logMass)
colnames(y_female) <- logMass
Estimate <- colMeans(y_female)
Q2.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_female, 2, function(x) sd(x))
XXf.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5) %>% 
  mutate(sex = "XXf")
# XX males (sex-reversed)
a_XX_male_xx <- a_female + post_Bas_m1$b_sexXXm 
b_XX_male_xx <- b_mass_female + post_Bas_m1$`b_sexXXm:logMass`
y_XX_male <- a_XX_male_xx + b_XX_male_xx %*% t(logMass)
colnames(y_XX_male) <- logMass
Estimate <- colMeans(y_XX_male)
Q2.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XX_male, 2, function(x) sd(x))
XXm.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XXm")
# XY males
a_XY_male <- a_female + post_Bas_m1$b_sexXYm 
b_XY_male <- b_mass_female + post_Bas_m1$`b_sexXYm:logMass`
y_XY_male <- a_XY_male + b_XY_male %*% t(logMass)
colnames(y_XY_male) <- logMass
Estimate <- colMeans(y_XY_male)
Q2.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XY_male, 2, function(x) sd(x))
XYm.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XYm")
# combine predictions
bass.mannual.pred <- rbind(XXf.man.pred, XXm.man.pred, XYm.man.pred)


############
# XX females posterior predictions
############
head(post_Bas_m1)
a_female_SD <- post_Bas_m1$b_Intercept
b_mass_female_SD <- post_Bas_m1$b_logMass
logMass = c(-0.28, 0, 0.28)
y_female_SD <- a_female_SD + b_mass_female_SD %*% t(logMass)
XXf.SD <- data.frame(y_female_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 ) %>% 
  mutate(sex = "XXf")
############
# XX posterior predictions
############
a_XX_male_xx_SD <- a_female_SD + post_Bas_m1$b_sexXXm 
b_XX_male_xx_SD <- b_mass_female_SD + post_Bas_m1$`b_sexXXm:logMass`
y_XX_male_SD <- a_XX_male_xx_SD + b_XX_male_xx_SD %*% t(logMass)
XXm.SD <- data.frame(y_XX_male_SD)%>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "XXm")
############
# XY male SD 
############
a_XY_male_SD <- a_female_SD + post_Bas_m1$b_sexXYm 
b_XY_male_SD <- b_mass_female_SD + post_Bas_m1$`b_sexXYm:logMass`
y_XY_male_SD <- a_XY_male_SD + b_XY_male_SD %*% t(logMass)
XYm.SD <- data.frame(y_XY_male_SD)%>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 ) %>% 
  mutate(sex = "XYm")
# combine predictions and save for analysis & Figures
bass.SD <- rbind(XXf.SD, XXm.SD, XYm.SD) 
# rename df so that tests are not numbers
bass.SD.Pred <- bass.SD %>% 
  dplyr::rename("lower" = 1,
                "mean" = 2 ,
                "upper" = 3)


############## 
# df for summarizing raw data points for Figure 2A 
############## 
# add in the fitted value from the model.
bassiana.data2<- data.frame(cbind(bassiana.data, predict(Bas_m1_brms_smr)))
#Summarise raw data for figure 2 
bass.raw.summary <- bassiana.data2 %>% 
  group_by(day, id, sex) %>% 
  summarise(MR = mean(Estimate),
            MR.se = std.error(Estimate), 
            logMass = mean(logMass),
            zstartmass = mean(zstartmass), 
            zendmass = mean(zendmass)) %>% 
  mutate(a = ((zstartmass - logMass)^2),
         b = ((zstartmass - logMass)^2), 
         c = (a+b)/2,
         d = c/2,
         sd = sqrt(d),
         mass.se = sd/(sqrt(2))) %>% 
  dplyr::select(-a,-b,-c,-d)

####################
# ALL Bassiana Plots
#######################
# Regression plot with predicted line and body mass
mycolors <- c("#333333", "#990000", "#3399FF")
bass.reg<-  ggplot(data =bassiana.data2 , aes(logMass, Estimate, group = sex, color= sex )) +
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = bass.mannual.pred, aes(x=logMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = bass.mannual.pred, aes(x=logMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) + 
  geom_smooth(data = bass.mannual.pred, aes(x=logMass, y=Estimate))+
  geom_vline(xintercept = 0.28, colour="black", linetype = "dotdash") + # Adding SD
  geom_text(aes(x=.28, label="+1.5 SD\n", y=-3.5), colour="black", angle=90, text=element_text(size=5))+
  geom_vline(xintercept = -0.28, colour="black", linetype = "dotdash") + # adding SD
  geom_text(aes(x=-.28, label="-1.5 SD\n", y=-3.5), colour="black", angle=90, text=element_text(size=5))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=10),axis.title = element_text(face="bold"))+
  scale_y_continuous(breaks = seq (-6.0, -3.0, by = .50), limits=c(-6.0, -3.0))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
bass.reg.plot <- ggMarginal(bass.reg, margins = "x", groupColour = TRUE, groupFill = TRUE)

############
# density plot of regression +- SD of mean for bassiana
############
# SD Plot
# stacking SD predict value category for plot
# now 
Bass.Sd.plot.data <-bass.SD %>% # renaming for figure 
  rename("-1.5 SD" = 1,
         "Mean" = 2, 
         "+1.5 SD" = 3 ) %>% 
  melt(id.vars="sex") %>% # this stacks data by predictions from each test by sex 
  rename(test = variable,
         Estimate = value)
Bass.Sd.plot.data$test <- factor(Bass.Sd.plot.data$test, c("+1.5 SD", "Mean", "-1.5 SD"))
mycolors <- c("#333333", "#990000", "#3399FF")
bass.sd.plot <- ggplot(Bass.Sd.plot.data, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., switch="y")+
  scale_y_continuous(position = "right")+
  xlab("Predicted Mean Metabolic Rate")+
  scale_x_continuous(name="Predicted Metabolic Rate", breaks = seq (-6.0, -3.6, by=0.2), limits=c(-6.0, -3.6))+
  theme_bw()+
  theme(axis.text = element_text(size=10), axis.title = element_text(face="bold"))
# combined plots
bassiana.final.fig <- cowplot::plot_grid(bass.reg.plot, bass.sd.plot, labels = c('A', 'B'), ncol=2)



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
# Model 1 & Model 2 comparison
####################
# SE is large and difference is not sufficient to warrant a heteroskedastic model 
loo_compare(Pog_m1_brms_smr, Pog_m2_brms_smr)

####################  
# extract posteriors for plotting and hypothesis testing for interaction (sex*mass) 
####################
post_pog_m2 <- posterior_samples(Pog_m2_brms_smr, pars = "^b")
Pog_m2_brms_smr
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
color_scheme_set("darkgray")
mcmc_areas(pog.dat, 
           pars = c("ZWf.mass.posterior", "ZZf.mass.posterior", "ZZm.mass.posterior"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=10)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 

####################  
# manual predict values for regression lines - Figure 2C
####################
# ZW females
a_female <- post_pog_m2$b_Intercept
b_mass_female <- post_pog_m2$b_logMass
logMass = seq(-0.5015452, 1.262014, 
              length.out = 100)
y_female <- a_female + b_mass_female %*% t(logMass)
colnames(y_female) <- logMass
Estimate <- colMeans(y_female)
Q2.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_female, 2, function(x) sd(x))
ZWf.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5) %>% 
  mutate(sex = "ZWf")
# ZZf
a_ZZ_female <- a_female + post_pog_m2$b_sexZZf 
b_ZZ_female <- b_mass_female + post_pog_m2$`b_sexZZf:logMass`
y_ZZ_female <- a_ZZ_female + b_ZZ_female %*% t(logMass)
colnames(y_ZZ_female) <- logMass
Estimate <- colMeans(y_ZZ_female)
Q2.5 <- apply(y_ZZ_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_ZZ_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_ZZ_female, 2, function(x) sd(x))
ZZf.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "ZZf")
# ZZm
a_ZZ_male <- a_female + post_pog_m2$b_sexZZm 
b_ZZ_male <- b_mass_female + post_pog_m2$`b_sexZZm:logMass`
y_ZZ_male <- a_ZZ_male + b_ZZ_male %*% t(logMass)
colnames(y_ZZ_male) <- logMass
Estimate <- colMeans(y_ZZ_male)
Q2.5 <- apply(y_ZZ_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_ZZ_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_ZZ_male, 2, function(x) sd(x))
ZZm.man.pred <- data.frame(logMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "ZZm")
# combined predictions
pog.mannual.pred <- rbind(ZZm.man.pred, ZZf.man.pred, ZWf.man.pred)

############## 
# df for summarizing raw datapoints for Figure 2C
############## 
pog.pre<- data.frame(predict(Pog_m2_brms_smr))
pogona.data2<- cbind(pogona.data, pog.pre)
#Summarise 
pogona.raw.summary <- pogona.data %>% 
  group_by(day, id, sex) %>% 
  summarise(MR = mean(log(O2_min)),
            MR.se = std.error(log(O2_min)), 
            logMass = mean(logMass),
            zstartmass = mean(zstartmass), 
            zendmass = mean(zendmass)) %>% 
  mutate(a = ((zstartmass - logMass)^2),
         b = ((zstartmass - logMass)^2), 
         c = (a+b),
         d = c/2,
         sd = sqrt(d),
         mass.se = sd/(sqrt(2)))%>% 
  dplyr::select(-a,-b,-c,-d)

############
# ZW female SD 
############
a_ZWfemale_SD <- post_pog_m2$b_Intercept
b_ZW_mass_female_SD <- post_pog_m2$b_logMass
logMass = c(-0.45, 0, 0.45)
y_ZWfemale_SD <- a_ZWfemale_SD + b_ZW_mass_female_SD %*% t(logMass)
ZWf.SD <- data.frame(y_ZWfemale_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "ZW_Female") 
############
# ZZ female SD 
############
a_ZZ_female_SD <- a_ZWfemale_SD + post_pog_m2$b_sexZZf 
b_ZZ_female_SD <- b_ZW_mass_female_SD + post_pog_m2$`b_sexZZf:logMass`
y_ZZ_female_SD <- a_ZZ_female_SD + b_ZZ_female_SD %*% t(logMass)
ZZf.SD <- data.frame(y_ZZ_female_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "ZZ_Female")
############
# ZZ male SD 
############
a_ZZ_male_SD <- a_ZWfemale_SD + post_pog_m2$b_sexZZm 
b_ZZ_male_SD <- b_ZW_mass_female_SD + post_pog_m2$`b_sexZZm:logMass`
y_ZZ_male_SD <- a_ZZ_male_SD + b_ZZ_male_SD %*% t(logMass)
ZZm.SD <- data.frame(y_ZZ_male_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "ZZ_Male")

# combine predictions and save for analysis
Pog.SD <- rbind(ZWf.SD, ZZf.SD, ZZm.SD)
Pog.SD.Pred <- Pog.SD %>% 
  dplyr::rename("-1.5_SD_pedict" = 1,
                "Mean_predict" = 2 ,
                "+1.5_SD_predict" = 3)

####################
# All Pogona Plots
####################
# color
mycolors <- c("#333333", "#990000", "#3399FF")
pog.reg <-  ggplot(data =pogona.data2 , aes(logMass, Estimate, group = sex, color= sex ))+
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = pog.mannual.pred, aes(x=logMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = pog.mannual.pred, aes(x=logMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = pog.mannual.pred, aes(x=logMass, y=Estimate))+
  geom_vline(xintercept = 0.45, colour="black", linetype = "dotdash") + # Adding SD line
  geom_text(aes(x=0.45, label="+1.5 SD\n", y=-0.5), colour="black", angle=90, text=element_text(size=10))+
  geom_vline(xintercept = -0.45, colour="black", linetype = "dotdash") + # adding SD line
  geom_text(aes(x= -0.45, label="-1.5 SD\n", y=-0.5), colour="black", angle=90, text=element_text(size=10))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  scale_y_continuous(breaks=c(0,-0.5, -1, -1.5, -2, -2.5, -3, -3.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(face="bold")) +
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
pog.reg.plot <- ggMarginal(pog.reg, margins = "x", groupColour = TRUE, groupFill = TRUE)
############
# density plots 
############
# setting up data for plots
Pog.Sd.plot.data <-Pog.SD %>% # renaming for figure 
  rename("-1.5 SD" = 1,
         "Mean" = 2, 
         "+1.5 SD" = 3 ) %>% 
  melt(id.vars="sex") %>% # this stacks data by predictions from each test by sex 
  rename(test = variable,
         Estimate = value)
# reordering data for plots
Pog.Sd.plot.data$test <- factor(Pog.Sd.plot.data$test, levels = c("+1.5 SD", "Mean", "-1.5 SD"))
# plot
mycolors <- c("#333333", "#990000", "#3399FF")
pog.sd.plot <- ggplot(Pog.Sd.plot.data, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., switch="y", scales="free_y")+
  scale_y_continuous(position = "right")+
  scale_x_continuous(name="Predicted Metabolic Rate", breaks = seq (-3.5,-1.0, by=0.5), limits = c(-3.5,-1.0))+
  theme_bw() +
  theme(axis.text = element_text(size=10),axis.title = element_text(face="bold"))
# combined figure
pogona.final.fig <- cowplot::plot_grid(pog.reg.plot,pog.sd.plot, labels = c('C', 'D'), ncol=2)

#combined pogona and bassiana figure:
grid.arrange(bassiana.final.fig, pogona.final.fig, nrow = 2)


