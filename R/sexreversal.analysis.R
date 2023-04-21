################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 
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
# quick plot
# box/violin plot
fig <- ggplot(bassiana.data, aes(x = sex, y = (O2_min))) + 
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
#####
## residual function - used to check residuals for all models
#####
residuals_brms <- function(model, data){
  fitted <- predict(model)[1:dim(data)[1], 1] # Predcit mean for each data
  e <- with(data, log(O2_min)) - fitted
  return(e)
}
##############
## Model 1 -  name:Bas_m1_brms
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Bas_m1_brms <- brm(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4))
  Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo"))
  saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
} else {Bas_m1_brms <- readRDS("./models/Bas_m1_brms")}
####################
# Model 1 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
draws <- as.array(Bas_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_logMass"), lags =10)
# check residuals
e <- residuals_brms(Bas_m1_brms, bassiana.data)
hist(e)
# plots
plot(Bas_m1_brms)
#R2 and summary of full model
bayes_R2(Bas_m1_brms)
summary(Bas_m1_brms)


##############
# Model 2 - Bas_m2_brms (heteroskedastic model)
##############
rerun2=FALSE
if(rerun2){
  mod_bas <- bf(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                sigma ~ logMass + ztime)
  Bas_m2_brms <- brm(mod_bas, family = "gaussian", 
                     data = bassiana.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
  Bas_m2_brms <- add_criterion(Bas_m2_brms, c("loo"))
  saveRDS(Bas_m2_brms, "./models/Bas_m2_brms")
} else {
  # read file in
  Bas_m2_brms <- readRDS(file = "models/Bas_m2_brms")
}
####################
# Model 2 checks: lags, residuals, r2, summary
####################
draws <- as.array(Bas_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_logMass"), lags =10)
# Check residuals
e <- residuals_brms(Bas_m2_brms, bassiana.data)
hist(e)
# plots
plot(Bas_m2_brms)
#R2 and summary
bayes_R2(Bas_m2_brms)
summary(Bas_m2_brms)


####################
# Model 1 & Model 2 comparison
####################
# SE is large and difference is not sufficient to warrant a heteroskedastic model 
loo_compare(Bas_m1_brms, Bas_m2_brms)


####################  
# Model 2 - pior table
####################  
prior_summary(Bas_m1_brms)

####################  
# extract posteriors for plotting and hypothesis testing for interaction (sex*mass) 
####################
summary(Bas_m1_brms)
post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")
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
Bass.RslopeDiff.Pheno.mass

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


#############
# MR posterior predictions at 1SD +/1 mean for density plots for Figure 2B
#############
SD_values <- bassiana.data %>% 
  summarise(mean = mean(logMass),
            sd =  sd(logMass),
            SD1.5_above = mean + (sd*1.5),
            SD1.5_below = mean - (sd*1.5)) %>% 
  mutate(across(1:4, round, 2))

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
saveRDS(bass.SD.Pred, "final.analysis.data/Bassiana.SD.mod.dat.RDS")


# LIKE Phenotype SD DIFFERENCES
# +1.5 SD contrasts
contrats_1.5SD = XXm.SD[,3] - XYm.SD[,3]
contrats_1.5SD_male <- as.array(contrats_1.5SD)
bass_mean_upper <- mean(contrats_1.5SD)
bass_error_upper <- round(quantile(contrats_1.5SD, c(0.025, 0.975)), digits = 2)
bass_LikeP_pmcmc_upper <- Bass_likeP_upper<-pmcmc(contrats_1.5SD_male)
bass_likeP_upper_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_likeP_upper_tbl[1,1] <- "Like Phenotype +1.5"
bass_likeP_upper_tbl[1,2] <- round(bass_mean_upper, digits = 2)
bass_likeP_upper_tbl[1,3] <- "(-0.19 - 0.23)"
bass_likeP_upper_tbl[1,4] <- round(bass_LikeP_pmcmc_upper, digits = 2)
# mean contrasts
bass_contrats_mean = XXm.SD[,2] - XYm.SD[,2]
contrats_mean_male <- as.array(bass_contrats_mean)
bass_error_mean <- mean(bass_contrats_mean)
bass_error_upper <- round(quantile(bass_contrats_mean, c(0.025, 0.975)), digits = 2)
bass_LikeP_pmcmc_mean <- round(pmcmc(contrats_mean_male), digits = 2)
bass_likeP_mean_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_likeP_mean_tbl[1,1] <- "Like Phenotype mean"
bass_likeP_mean_tbl[1,2] <- round(bass_error_mean, digits = 2)
bass_likeP_mean_tbl[1,3] <- "(-0.21 - 0.15)"
bass_likeP_mean_tbl[1,4] <- round(bass_LikeP_pmcmc_mean, digits = 2)
# -1.5 SD contrasts
contrats_lower = XXm.SD[,1] - XYm.SD[,1]
contrats_lower_male <- as.array(contrats_lower)
bass_error_lower_mean <- mean(contrats_lower)
bass_error_lower <- round(quantile(contrats_lower, c(0.025, 0.975)), digits = 2)
bass_error_lower_pmcmc <-  pmcmc(contrats_lower_male)
bass_likeP_lower_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_likeP_lower_tbl[1,1] <- "Like Phenotype -1.5"
bass_likeP_lower_tbl[1,2] <- round(bass_error_lower_mean, digits = 2)
bass_likeP_lower_tbl[1,3] <- "(-0.29 - 0.13)"
bass_likeP_lower_tbl[1,4] <- round(bass_error_lower_pmcmc, digits = 2)
# rbind bass like Phenotype
bass_contrast_errors_LikeP <- rbind(bass_likeP_lower_tbl, 
                              bass_likeP_mean_tbl, 
                              bass_likeP_upper_tbl) %>% 
  dplyr::rename(Hypothesis = X1,
                Estimate = X2,
                "Estimate error" = X3,
                "Estimate pMCMC" = X4)

# LIKE Genotype SD DIFFERENCES
# +1.5 SD contrasts
bass_contrats_upper = XXm.SD[,3] - XXf.SD[,3]
bass_contrats_upper_genotype <- as.array(bass_contrats_upper)
bass_LikeG_upper_mean <- mean(bass_contrats_upper)
bass_LikeG_upper_mean_error <- round(quantile(bass_contrats_upper, c(0.025, 0.975)), digits = 2)
bass_likeG_upper_pmcmc<-pmcmc(bass_contrats_upper_genotype)
bass_LikeG_upper_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_LikeG_upper_tbl[1,1] <- "Like Genotype +1.5"
bass_LikeG_upper_tbl[1,2] <- round(bass_LikeG_upper_mean, digits = 2)
bass_LikeG_upper_tbl[1,3] <- "(-0.50 - -0.11)"
bass_LikeG_upper_tbl[1,4] <- round(bass_likeG_upper_pmcmc, digits = 2)
# mean contrasts
bass_contrats_mean = XXm.SD[,2] - XXf.SD[,2]
bass_contrats_mean_genotype <- as.array(bass_contrats_mean)
bass_LikeG_mean<- mean(bass_contrats_mean)
bass_LikeG_mean_error <- round(quantile(bass_contrats_mean, c(0.025, 0.975)), digits = 2)
bass_likeG_mean_pmcmc<-pmcmc(bass_contrats_mean_genotype)
bass_LikeG_mean_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_LikeG_mean_tbl[1,1] <- "Like Genotype mean"
bass_LikeG_mean_tbl[1,2] <- round(bass_LikeG_mean, digits = 2)
bass_LikeG_mean_tbl[1,3] <- "(-0.32 - 0.02)"
bass_LikeG_mean_tbl[1,4] <- round(bass_likeG_mean_pmcmc, digits = 2)
# -1.5 SD contrasts
bass_contrats_lower = XXm.SD[,1] - XXf.SD[,1]
bass_contrats_lower_genotype <- as.array(bass_contrats_lower)
bass_LikeG_lower <- mean(bass_contrats_lower)
bass_LikeG_lower_error <- round(quantile(bass_contrats_lower, c(0.025, 0.975)), digits = 2)
bass_likeG_lower_pmcmc<-pmcmc(bass_contrats_lower_genotype)
bass_LikeG_lower_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
bass_LikeG_lower_tbl[1,1] <- "Like Genotype -1.5"
bass_LikeG_lower_tbl[1,2] <- round(bass_LikeG_lower, digits = 2)
bass_LikeG_lower_tbl[1,3] <- "(-0.19  - 0.21)"
bass_LikeG_lower_tbl[1,4] <- round(bass_likeG_lower_pmcmc, digits = 2)
# rbind
# rbind bass like Phenotype
bass_contrast_errors_LikeG <- rbind(bass_LikeG_lower_tbl, 
                                    bass_LikeG_mean_tbl, 
                                    bass_LikeG_upper_tbl) %>% 
  dplyr::rename(Hypothesis = X1,
                Estimate = X2,
                "Estimate error" = X3,
                "Estimate pMCMC" = X4)

# final contrast df for likeG or LikeP 
bass_error_contrast_fig2 <- rbind(bass_contrast_errors_LikeP,
                                  bass_contrast_errors_LikeG) %>% 
  mutate(Species = "Bassiana duperreyi")


############## 
# df for summarizing raw data points for Figure 2A 
############## 
# add in the fitted value from the model.
bassiana.data2<- data.frame(cbind(bassiana.data, predict(Bas_m1_brms)))
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

#############
# body mass differences test
#############
bass.bodymass <- brm(logMass ~ sex , data = bass.raw.summary)
saveRDS(bass.bodymass, "./models/Bas_bodymass")
summary(bass.bodymass)


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
  scale_x_continuous(name="Predicted Metabolic Rate", breaks = seq (-5.8, -3.2, by=0.2), limits=c(-5.8, -3.2))+
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
# quick plot
# box/violin plot
fig <- ggplot(pogona.data.smr, aes(x = sex, y = log(O2_min))) + 
  geom_point() +
  geom_violin() + 
  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
  geom_jitter(fill = "grey", alpha = 0.3) + 
  labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "Sex")+ 
  theme_bw()
fig

##############
## Model 1 - name:Pog_m1_brms
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Pog_m1_brms <- brm(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4))
  Pog_m1_brms <- add_criterion(Pog_m1_brms, c("loo"))
  saveRDS(Pog_m1_brms, "./models/Pog_m1_brms")
} else {Pog_m1_brms <- readRDS("./models/Pog_m1_brms")}
####################
# Model 1 checks: lags, residuals, r2, summary
####################
draws <- as.array(Pog_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept", "b_sexZZf", "b_sexZZm"), lags =10)
# check residuals
e <- residuals_brms(Pog_m1_brms, pogona.data)
hist(e)
# plots
plot(Pog_m1_brms)
#R2 of and summary
bayes_R2(Pog_m1_brms)
summary(Pog_m1_brms)


##############
## Model 2 - Pog_m2_brms (heteroskedastic model)
##############
rerun2=FALSE
if(rerun2){
  mod <- bf(log(O2_min) ~ sex*logMass + ztime + (1 + ztime | id) + (1 | day),  
            sigma ~ logMass + ztime)
  Pog_m2_brms <- brm(mod, family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
  Pog_m2_brms <- add_criterion(Pog_m2_brms, c("loo"))
  saveRDS(Pog_m2_brms, "./models/Pog_m2_brms")
} else {
  # read file in
  Pog_m2_brms <- readRDS(file = "./models/Pog_m2_brms")
}
####################
# Model 2 - checks: lags, residuals, r2, summary
####################
draws <- as.array(Pog_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept", "b_sexZZf", "b_sexZZm"), lags =10)
# check residuals
e <- residuals_brms(Pog_m2_brms, pogona.data)
hist(e)
# plots
plot(Pog_m2_brms)
#R2 of and summary
bayes_R2(Pog_m2_brms)
summary(Pog_m2_brms)


####################
# Model 1 & Model 2 comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)

####################  
# Model 2 - pior table
####################  
prior_summary(Pog_m2_brms)

####################  
# extract posteriors for plotting and hypothesis testing for interaction (sex*mass) 
####################
post_pog_m2 <- posterior_samples(Pog_m2_brms, pars = "^b")
Pog_m2_brms
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

# combined predictions
pog.mannual.pred <- rbind(ZZm.man.pred, ZZf.man.pred, ZWf.man.pred)


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

############## 
# df for summarizing raw datapoints for Figure 2C
############## 
pog.pre<- data.frame(predict(Pog_m2_brms))
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


#############
# test for differences in body mass
#############
pog.bodymass <- brm(logMass ~ sex , data = pogona.raw.summary)
saveRDS(pog.bodymass, "./models/Pog_bodymass")

#############
# Calculating 1SD +/1 mean figure 2D
#############
Pog_SD_values <- pogona.data2 %>% 
  ungroup() %>% 
  summarise(mean = mean(logMass),
            sd =  sd(logMass),
            SD1.5_above = mean + (sd*1.5),
            SD1.5_below = mean - (sd*1.5)) %>% 
  mutate(across(1:4, round, 2))
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
saveRDS(Pog.SD.Pred, "final.analysis.data/PogonaSD.mod.dat.RDS")

# LIKE Phenotype SD DIFFERENCES
# +1.5 SD contrasts
pog_contrats_upper = ZZf.SD[,3] - ZWf.SD[,3]
pog_contrats_upper_genotype <- as.array(pog_contrats_upper)
pog_LikeP_upper_mean <- mean(pog_contrats_upper)
pog_LikeP_upper_mean_error <- round(quantile(pog_contrats_upper, c(0.025, 0.975)), digits = 2)
pog_LikeP_upper_pmcmc<-pmcmc(pog_contrats_upper_genotype)
pog_LikeP_upper_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeP_upper_tbl[1,1] <- "Like Phenotype +1.5"
pog_LikeP_upper_tbl[1,2] <- round(pog_LikeP_upper_mean, digits = 2)
pog_LikeP_upper_tbl[1,3] <- "(-0.37 - -0.03)"
pog_LikeP_upper_tbl[1,4] <- round(pog_LikeP_upper_pmcmc, digits = 2)
# mean contrasts
pog_contrats_mean = ZZf.SD[,2] - ZWf.SD[,2]
pog_contrats_mean_genotype <- as.array(pog_contrats_mean)
pog_LikeP_mean<- mean(pog_contrats_mean)
pog_LikeP_mean_error <- round(quantile(pog_contrats_mean, c(0.025, 0.975)), digits = 2)
pog_LikeP_mean_pmcmc<-pmcmc(pog_contrats_mean_genotype)
pog_LikeP_mean_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeP_mean_tbl[1,1] <- "Like Phenotype mean"
pog_LikeP_mean_tbl[1,2] <- round(pog_LikeP_mean, digits = 2)
pog_LikeP_mean_tbl[1,3] <- "(-0.28 - 0.03)"
pog_LikeP_mean_tbl[1,4] <- round(pog_LikeP_mean_pmcmc, digits = 2)
# -1.5 SD contrasts
pog_contrats_lower = ZZf.SD[,1] - ZWf.SD[,1]
pog_contrats_lower_genotype <- as.array(pog_contrats_lower)
pog_LikeP_lower <- mean(pog_contrats_lower)
pog_LikeP_lower_error <- round(quantile(pog_contrats_lower, c(0.025, 0.975)), digits = 2)
pog_LikeP_lower_pmcmc<-pmcmc(pog_contrats_lower_genotype)
pog_LikeP_lower_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeP_lower_tbl[1,1] <- "Like Phenotype -1.5"
pog_LikeP_lower_tbl[1,2] <- round(pog_LikeP_lower, digits = 2)
pog_LikeP_lower_tbl[1,3] <- "(-0.23  - 0.11)"
pog_LikeP_lower_tbl[1,4] <- round(pog_LikeP_lower_pmcmc, digits = 2)
# rbind bass like Phenotype
pog_contrast_errors_LikeP <- rbind(pog_LikeP_lower_tbl, 
                                    pog_LikeP_mean_tbl, 
                                    pog_LikeP_upper_tbl) %>% 
  dplyr::rename(Hypothesis = X1,
                Estimate = X2,
                "Estimate error" = X3,
                "Estimate pMCMC" = X4)

# LIKE Genotype SD DIFFERENCES
# +1.5 SD contrasts
pog_contrats_upper = ZZf.SD[,3] - ZZm.SD[,3]
pog_contrats_upper_genotype <- as.array(pog_contrats_upper)
pog_LikeG_upper_mean <- mean(pog_contrats_upper)
pog_LikeG_upper_mean_error <- round(quantile(pog_contrats_upper, c(0.025, 0.975)), digits = 2)
pog_LikeG_upper_pmcmc<-pmcmc(pog_contrats_upper_genotype)
pog_LikeG_upper_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeG_upper_tbl[1,1] <- "Like Genotype +1.5"
pog_LikeG_upper_tbl[1,2] <- round(pog_LikeG_upper_mean, digits = 2)
pog_LikeG_upper_tbl[1,3] <- "(-0.12 - 0.19)"
pog_LikeG_upper_tbl[1,4] <- round(pog_LikeG_upper_pmcmc, digits = 2)
# mean contrasts
pog_contrats_mean = ZZf.SD[,2] - ZZm.SD[,2]
pog_contrats_mean_genotype <- as.array(pog_contrats_mean)
pog_LikeG_mean<- mean(pog_contrats_mean)
pog_LikeG_mean_error <- round(quantile(pog_contrats_mean, c(0.025, 0.975)), digits = 2)
pog_LikeG_mean_pmcmc<-pmcmc(pog_contrats_mean_genotype)
pog_LikeG_mean_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeG_mean_tbl[1,1] <- "Like Genotype mean"
pog_LikeG_mean_tbl[1,2] <- round(pog_LikeG_mean, digits = 2)
pog_LikeG_mean_tbl[1,3] <- "(-0.20 - 0.09)"
pog_LikeG_mean_tbl[1,4] <- round(pog_LikeG_mean_pmcmc, digits = 2)
# -1.5 SD contrasts
pog_contrats_lower = ZZf.SD[,1] - ZZm.SD[,1]
pog_contrats_lower_genotype <- as.array(pog_contrats_lower)
pog_LikeG_lower <- mean(pog_contrats_lower)
pog_LikeG_lower_error <- round(quantile(pog_contrats_lower, c(0.025, 0.975)), digits = 2)
pog_LikeG_lower_pmcmc<-pmcmc(pog_contrats_lower_genotype)
pog_LikeG_lower_tbl <- data.frame(matrix(ncol = 4, nrow = 1))
pog_LikeG_lower_tbl[1,1] <- "Like Genotype -1.5"
pog_LikeG_lower_tbl[1,2] <- round(pog_LikeG_lower, digits = 2)
pog_LikeG_lower_tbl[1,3] <- "(-0.31  - 0.00)"
pog_LikeG_lower_tbl[1,4] <- round(pog_LikeG_lower_pmcmc, digits = 2)
# rbind Pog like Genotype
pog_contrast_errors_LikeG <- rbind(pog_LikeG_lower_tbl, 
                                   pog_LikeG_mean_tbl, 
                                   pog_LikeG_upper_tbl) %>% 
  dplyr::rename(Hypothesis = X1,
                Estimate = X2,
                "Estimate error" = X3,
                "Estimate pMCMC" = X4)

# final contrast df for likeG or LikeP 
pog_error_contrast_fig2 <- rbind(pog_contrast_errors_LikeP,
                                  pog_contrast_errors_LikeG) %>% 
  mutate(Species = "Pogona vitticeps")

# rbind error for results section
figure_2B_2D_contrast_tbl <- rbind(bass_error_contrast_fig2,
                                   pog_error_contrast_fig2)

saveRDS(figure_2B_2D_contrast_tbl, file = "final.analysis.data/figure_2B_2D_contrast_tbl")


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
  scale_x_continuous(name="Predicted Metabolic Rate", breaks = seq (-3,-0.5, by=0.5), limits = c(-3,-0.5))+
  theme_bw() +
  theme(axis.text = element_text(size=10),axis.title = element_text(face="bold"))
# combined figure
pogona.final.fig <- cowplot::plot_grid(pog.reg.plot,pog.sd.plot, labels = c('C', 'D'), ncol=2)

#combined pogona and bassiana figure:
grid.arrange(bassiana.final.fig, pogona.final.fig, nrow = 2)


####################
# growth  comparisons when accounting 
####################
##########
# MR and growth rate
#########
# filter data with growth and MR and summarize MR
# O2 - Pogona
pogona.data <- read.csv("~/Dropbox/energy_sex_reversal/final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min,
         sex = Geno.pheno,
         mass_g = mass, 
         ) %>% 
  group_by(id, sex) %>% 
  summarise(mean_02 = mean(O2_min))
# O2 - Bassiana
bassiana.data <- read.csv("~/Dropbox/energy_sex_reversal/final.analysis.data/Bassiana.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min,
         sex = Geno.pheno,
         mass_g = mass) %>% 
  group_by(id, sex) %>% 
  summarise(mean_02 = mean(O2_min))
# growth & survival data both spp
growth <- read.csv(file = "~/Dropbox/energy_sex_reversal/final.analysis.data/growth.bassiana.pogona.csv") %>% 
  dplyr::rename(sex = Geno.pheno, 
                id = ID)
growth$status <- ifelse(growth$Dead.Y.N. == "Alive", 1, 0)
pogona_growth <- growth %>% 
  filter(Species == "Pogona")
bassiana_growth <- growth %>% 
  filter(Species == "Bassiana")
# merge data
pogona_final <- merge(x = pogona.data, y = pogona_growth, by = "id", all = TRUE) %>% 
  rename(sex = sex.x)
bassiana_final <- merge(x = bassiana.data, y = bassiana_growth, by = "id", all = TRUE) %>% 
  rename(sex = sex.x) %>% 
  mutate(Growth.rate.mass.cg = Growth.rate.mass.*100)


######
# pogona MR and growth analysis mass 
######
pog_o2_growth_mass <- brm(log(mean_02)~Growth.rate.mass.*sex , family = "gaussian", 
                     data = pogona_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 5, cores = 8) # warning for NA's are animals that died
# assumption plot checks
plot(pog_o2_growth_mass)
# summary & r2
summary(pog_o2_growth_mass)
bayes_R2(pog_o2_growth_mass)
# save model
saveRDS(pog_o2_growth_mass, "./models/Pog_growth_mass_metabolism")

# pogona MR and growth analysis SVL 
pog_o2_growth_svl <- brm(log(mean_02)~Growth.rate.SVL.*sex , family = "gaussian", 
                          data = pogona_final, 
                          iter= 5000, warmup = 1000, 
                          thin = 5, cores = 8) # warning for NA's are animals that died
# assumption plot checks
plot(pog_o2_growth_svl)
# summary & r2
summary(pog_o2_growth_svl)
bayes_R2(pog_o2_growth_svl)
# save model
saveRDS(pog_o2_growth_svl, "./models/Pog_growth_svl_metabolism")


######
# Bassiana
######
# mass growth
bass_o2_growth_mass <- brm(log(mean_02)~ Growth.rate.mass.cg  * sex,
                      family = "gaussian", 
                      data = bassiana_final, 
                      iter= 5000, warmup = 1000, 
                      thin = 5, cores = 8) # warning for NA's are animals that died

# assumption plot checks
plot(bass_o2_growth_mass)
# summary & r2
summary(bass_o2_growth_mass)
bayes_R2(bass_o2_growth_mass)
# save model
saveRDS(bass_o2_growth_mass, "./models/Bass_growth_mass_metabolism")
# SVL
bass_o2_growth_svl <- brm(log(mean_02)~Growth.rate.SVL.  * sex,
                           family = "gaussian", 
                           data = bassiana_final, 
                           iter= 5000, warmup = 1000, 
                           thin = 5, cores = 8) # warning for NA's are animals that died

# assumption plot checks
plot(bass_o2_growth_svl)
# summary & r2
summary(bass_o2_growth_svl)
bayes_R2(bass_o2_growth_svl)

# save model
saveRDS(bass_o2_growth_svl, "./models/Bass_growth_svl_metabolism")

####################
# Chisqure for mortality comparisons: summary,  analysis and figures for each species
####################
# mortality frequency summary Bassiana
Bassiana_chi_sq_final <- bassiana.growth %>% 
  group_by(Dead.Y.N., sex) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))
# mortality chi_square Pogona
Pogona_chi_sq_final <- pogona.growth %>% 
  group_by(Dead.Y.N., sex) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))

# analysis - Bassiana
bassiana.table<- table(bassiana.growth$Dead.Y.N., bassiana.growth$sex)
fisher.test(bassiana.table)
# analysis - Pogona
pogona.table <- table(pogona.growth$Dead.Y.N., pogona.growth$sex)
fisher.test(pogona.table)

## Survival Figures
# proportions figure - Bassiana
bass.mortality <- ggplot(Bassiana_chi_sq_final, aes(fill=sex, y=n, x=Dead.Y.N.)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("darkgrey", "red"))+
  scale_fill_manual(values=c("darkgrey", "red", "#56B4E9"))+
  ylab("% of survival | mortality") +
  xlab("Dead or Alive")+
  theme_bw(base_size = 18, base_family = "Times New Roman") 
# proportions figure - Pogona
pog.mortality <- ggplot(Pogona_chi_sq_final, aes(fill=sex, y=n, x=Dead.Y.N.)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("darkgrey", "red"))+
  scale_fill_manual(values=c("darkgrey", "red", "#56B4E9"))+
  ylab("% of survival | mortality") +
  xlab("Dead or Alive")+
  theme_bw(base_size = 18, base_family = "Times New Roman") 
# both figures
grid.arrange(bass.mortality, pog.mortality, nrow = 2)

















######## like phenotype hypothesis extract - OLD
# extract 
pog_post <- posterior_samples(pog_o2_growth)
pog_surv_brms_m1 <- posterior_samples(pog_o2_growth, pars = "^b")
dimnames(pog_surv_brms_m1)

## extracting posteriors for interaction of sex and mass
ZWf.mass.posterior <- as.array(pog_surv_brms_m1[,"b_Intercept"])
ZZf.mass.posterior <- as.array(ZWf.mass.posterior + pog_surv_brms_m1[,"b_sexZZf"])
ZZm.mass.posterior <- as.array(ZWf.mass.posterior + pog_surv_brms_m1[,"b_sexZZm"])

## H1: Like Phenotype Hypothesis - pMCMC
Pog.RslopeDiff.Pheno.mass <- ZWf.mass.posterior - ZZf.mass.posterior
Pog.pMCMC_phenotype_metabolism <- pmcmc(Pog.RslopeDiff.Pheno.mass)
Pog.pMCMC_phenotype_metabolism

## H2: Like Genotype Hypothesis - pMCMC
Pog.RslopeDiff.Geno.mass <- ZZf.mass.posterior - ZZm.mass.posterior
Pog.pMCMC_genotype_metabolism <- pmcmc(Pog.RslopeDiff.Geno.mass)
Pog.pMCMC_genotype_metabolism

# overall O2 on survival - high o2 high survival; dead = 0 alive = 1
plot(conditional_effects(pog_o2_growth, "mean_02"), ask = FALSE) 



