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
         zlogMass = scale(log(mass_g), scale = FALSE),
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
## residual function - used for all models
#####
residuals_brms <- function(model, data){
  fitted <- predict(model)[1:dim(data)[1], 1] # Predcit mean for each data
  e <- with(data, log(O2_min)) - fitted
  return(e)
}
##############
## Model 1 -  Bas_m1_brms
# @ dan updated the model so that ESS values were near 3000 so increased thinning to 5 and upped the iterations
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Bas_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4))
  Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo"))
  saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
} else {Bas_m1_brms <- readRDS("./models/Bas_m1_brms")}
####################
# Model 1 -  Bas_m1_brms checks
####################
# checking lags for this model 
# @dan looks  better after increasing thinning and iterations  
draws <- as.array(Bas_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_zlogMass"), lags =10)
# check residuals
e <- residuals_brms(Bas_m1_brms, bassiana.data)
hist(e)
# plots
plot(Bas_m1_brms)
#R2 and summary of full model
bayes_R2(Bas_m1_brms)
summary(Bas_m1_brms)


##############
# Model 2 - Bas_m2_brms
# @ dan updated the model so that ESS values were near 3000 so increased thinning to 5 and upped the iterations
##############
rerun2=FALSE
if(rerun2){
  mod_bas <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                sigma ~ zlogMass + ztime)
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
# Model 2 - Bas_m2_brms checks
####################
# checking lags for this model
# @ dan looks much better after increasing thinning and iterations  
draws <- as.array(Bas_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXXm", "b_sexXYm", "b_zlogMass"), lags =10)
# Check residuals
e <- residuals_brms(Bas_m2_brms, bassiana.data)
hist(e)
# plots
plot(Bas_m2_brms)
#R2 and summary
bayes_R2(Bas_m2_brms)
summary(Bas_m2_brms)

####################
# Model comparison
####################
# SE is large and difference is probably not sufficient to warrant a het model
loo_compare(Bas_m1_brms, Bas_m2_brms)

####################  
# extract posteriors for plotting and hypothesis testing
####################
summary(Bas_m1_brms)
post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")
variable.names(post_Bas_m1)

## Hypothesis testing - contrast OVERALL sex differences
## extracting posteriors for just sex
XXf.posterior <- as.array(post_Bas_m1[,"b_Intercept"])
XXm.posterior <- as.array(post_Bas_m1[,"b_sexXXm"] + XXf.posterior ) 
XYm.posterior <- as.array(post_Bas_m1[,"b_sexXYm"] + XXf.posterior)

# H1: Like Phenotype Hypothesis
RslopeDiff.Pheno <- XYm.posterior - XXm.posterior
rope(RslopeDiff.Pheno, ci = 0.95)
plot(p_direction(RslopeDiff.Pheno))
pd.Pheno <- p_direction(RslopeDiff.Pheno)
pd_to_p(pd.Pheno, direction = "two-sided")
# H2: Like Genotype Hypothesis
RslopeDiff.Geno <- XXf.posterior - XXm.posterior
rope(RslopeDiff.Geno, ci = 0.95)
plot(p_direction(RslopeDiff.Geno))
pd.geno <- p_direction(RslopeDiff.Geno)
pd_to_p(pd.geno, direction = "two-sided")


# Hypothesis testing - accounting for INTERACTION of sex and mass
## extracting posteriors for interaction of sex and mass
XXf.mass.posterior <- as.array(post_Bas_m1[,"b_Intercept"] + post_Bas_m1[,"b_zlogMass"])
XXm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXXm:zlogMass"])
XYm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXYm:zlogMass"])

# H1: Like Phenotype Hypothesis
RslopeDiff.Pheno.mass <- XYm.mass.posterior - XXm.mass.posterior
rope(RslopeDiff.Pheno.mass, ci = 0.95)
plot(p_direction(RslopeDiff.Pheno.mass))
pd.Pheno.mass <- p_direction(RslopeDiff.Pheno.mass)
pd_to_p(pd.Pheno.mass, direction = "two-sided")

# H2: Like Genotype Hypothesis
RslopeDiff.Geno.mass <- XXf.mass.posterior - XYm.mass.posterior
rope(RslopeDiff.Geno.mass, ci = 0.95)
plot(p_direction(RslopeDiff.Geno.mass))
pd.geno.mass <- p_direction(RslopeDiff.Geno.mass)
pd_to_p(pd.geno.mass, direction = "two-sided")

## plotting posteriors
bass.dat <- cbind(XXf.mass.posterior, XXm.mass.posterior, XYm.mass.posterior)
mcmc_areas(bass.dat, 
           pars = c("XXf.mass.posterior", "XXm.mass.posterior", "XYm.mass.posterior"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope Differences") 

####################  
# manual predict values for regression lines
####################
# XX females
a_female <- post_Bas_m1$b_Intercept
b_mass_female <- post_Bas_m1$b_zlogMass
zlogMass = seq(-0.3856391, 0.3710889, 
               length.out = 100)
y_female <- a_female + b_mass_female %*% t(zlogMass)
colnames(y_female) <- zlogMass
Estimate <- colMeans(y_female)
Q2.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_female, 2, function(x) sd(x))
XXf.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5) %>% 
  mutate(sex = "XXf")

# XXm
a_XX_male_xx <- a_female + post_Bas_m1$b_sexXXm 
b_XX_male_xx <- b_mass_female + post_Bas_m1$`b_sexXXm:zlogMass`
y_XX_male <- a_XX_male_xx + b_XX_male_xx %*% t(zlogMass)
colnames(y_XX_male) <- zlogMass
Estimate <- colMeans(y_XX_male)
Q2.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XX_male, 2, function(x) sd(x))
XXm.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XXm")

# XYm
a_XY_male <- a_female + post_Bas_m1$b_sexXYm 
b_XY_male <- b_mass_female + post_Bas_m1$`b_sexXYm:zlogMass`
y_XY_male <- a_XY_male + b_XY_male %*% t(zlogMass)
colnames(y_XY_male) <- zlogMass
Estimate <- colMeans(y_XY_male)
Q2.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XY_male, 2, function(x) sd(x))
XYm.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XYm")
# combine predictions
bass.mannual.pred <- rbind(XXf.man.pred, XXm.man.pred, XYm.man.pred)


#############
# MR posterior predictions at 1SD +/1 mean
#############
SD_values <- bassiana.data %>% 
  summarise(mean = mean(zlogMass),
            sd =  sd(zlogMass),
            SD1.5_above = mean + (sd*1.5),
            SD1.5_below = mean - (sd*1.5)) %>% 
  mutate(across(1:4, round, 2))
############
# XX females posterior predictions
############
head(post_Bas_m1)
a_female_SD <- post_Bas_m1$b_Intercept
b_mass_female_SD <- post_Bas_m1$b_zlogMass
zlogMass = c(-0.28, 0, 0.28)
y_female_SD <- a_female_SD + b_mass_female_SD %*% t(zlogMass)
XXf.SD <- data.frame(y_female_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 ) %>% 
  mutate(sex = "XXf")

############
# XX posterior predictions
############
a_XX_male_xx_SD <- a_female_SD + post_Bas_m1$b_sexXXm 
b_XX_male_xx_SD <- b_mass_female_SD + post_Bas_m1$`b_sexXXm:zlogMass`
y_XX_male_SD <- a_XX_male_xx_SD + b_XX_male_xx_SD %*% t(zlogMass)
XXm.SD <- data.frame(y_XX_male_SD)%>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "XXm")

############
# XY male SD 
############
a_XY_male_SD <- a_female_SD + post_Bas_m1$b_sexXYm 
b_XY_male_SD <- b_mass_female_SD + post_Bas_m1$`b_sexXYm:zlogMass`
y_XY_male_SD <- a_XY_male_SD + b_XY_male_SD %*% t(zlogMass)
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


############## 
# df for summarizing raw datapoints for ggplot 
############## 
# add in the fitted value from the model.
bassiana.data2<- data.frame(cbind(bassiana.data, predict(Bas_m1_brms)))
#Summarise 
bass.raw.summary <- bassiana.data2 %>% 
  group_by(day, id, sex) %>% 
  summarise(MR = mean(Estimate),
            MR.se = std.error(Estimate), 
            zlogMass = mean(zlogMass),
            zstartmass = mean(zstartmass), 
            zendmass = mean(zendmass)) %>% 
  mutate(a = ((zstartmass - zlogMass)^2),
         b = ((zstartmass - zlogMass)^2), 
         c = (a+b)/2,
         d = c/2,
         sd = sqrt(d),
         mass.se = sd/(sqrt(2))) %>% 
  dplyr::select(-a,-b,-c,-d)


#############
# test for differences in body mass from raw data file
#############
bass.bodymass <- brm(zlogMass ~ sex , data = bass.raw.summary)
saveRDS(bass.bodymass, "./models/Bas_bodymass")
summary(bass.bodymass)
bm_diff <- hypothesis(bass.bodymass, 'sexXXm + sexXYm = 0')
mean(bm_diff$samples$H1)
HPDinterval(mcmc(bm_diff$samples$H1))
# p-value equivalent for the full hypothesis test. 
pMCMC <- 1 - (table(bm_diff$samples$H1 > 0)[1] / (length(bm_diff$samples$H1) - 1))     
pMCMC

####################
# ALL Bassiana Plots
#######################
# Regression plot with predicted line and body mass
mycolors <- c("#333333", "#990000", "#3399FF")
bass.reg<-  ggplot(data =bassiana.data2 , aes(zlogMass, Estimate, group = sex, color= sex )) +
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate+Est.Error, colour = sex)) +
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate-Est.Error, colour = sex)) + # @Kris, to smooth just add in two more smoothed lines ontop
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate))+
  geom_vline(xintercept = 0.28, colour="black", linetype = "dotdash") + # Adding SD
  geom_text(aes(x=.28, label="+1.5 SD\n", y=-3.5), colour="black", angle=90, text=element_text(size=5))+
  geom_vline(xintercept = -0.28, colour="black", linetype = "dotdash") + # adding SD
  geom_text(aes(x=-.28, label="-1.5 SD\n", y=-3.5), colour="black", angle=90, text=element_text(size=5))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=10))+
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
  theme(axis.text = element_text(size=10))
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
         zlogMass = scale(log(mass_g), scale = FALSE),
         zstartmass = scale(log(start_mass_g), scale = FALSE),
         zendmass = scale(log(end_mass_g), scale = FALSE)) %>% 
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
##############
## Model 1 - Pog_m1_brms
# @dan control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 5000
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Pog_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4))
    # @dan got a warning about our waic values and it suggested to add moment_match = TRUE 
  Pog_m1_brms <- add_criterion(Pog_m1_brms, c("loo"))
  saveRDS(Pog_m1_brms, "./models/Pog_m1_brms")
} else {Pog_m1_brms <- readRDS("./models/Pog_m1_brms")}
####################
# Model 1 - Pog_m1_brms checks
####################
# checking lags for this model 
#@dan looks much better after increasing thinning and iterations  
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
## Model 2 - Pog_m2_brms
# @ dan same as above, control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 5000
##############
rerun2=FALSE
if(rerun2){
  mod <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
            sigma ~ zlogMass + ztime)
  Pog_m2_brms <- brm(mod, family = "gaussian", data = pogona.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
  Pog_m2_brms <- add_criterion(Pog_m2_brms, c("loo"))
  saveRDS(Pog_m2_brms, "./models/Pog_m2_brms")
} else {
  # read file in
  Pog_m2_brms <- readRDS(file = "./models/Pog_m2_brms")
}
####################
# Model 2 - Pog_m2_brms checks
####################
# checking lags for this model 
#@dan looks much better after increasing thinning and iterations  
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
# Model comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)
tab_model(Pog_m2_brms)

 ####################  
# extract posteriors for Pog_m2  + Plotting 
####################
post_pog_m2 <- posterior_samples(Pog_m2_brms, pars = "^b")
Pog_m2_brms
dimnames(post_pog_m2)
# extracting posteriors 
ZWf <- as.array(post_pog_m2[,"b_zlogMass"])
ZZf <- as.array(post_pog_m2[,"b_zlogMass"] + 
                  post_pog_m2[,"b_sexZZf:zlogMass"])
ZZm <- as.array(post_pog_m2[,"b_zlogMass"] + 
                  post_pog_m2[,"b_sexZZm:zlogMass"])
# combining to one df
pog.dat <- cbind(ZWf, ZZf, ZZm)
# plotting posteriors lot
mcmc_areas(pog.dat, 
           pars = c("ZWf", "ZZf", "ZZm"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=10)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 


####################  
# manual predict values for regression lines
####################
# ZW females
a_female <- post_pog_m2$b_Intercept
b_mass_female <- post_pog_m2$b_zlogMass
zlogMass = seq(-0.5015452, 1.262014, 
               length.out = 100)
y_female <- a_female + b_mass_female %*% t(zlogMass)
colnames(y_female) <- zlogMass
Estimate <- colMeans(y_female)
Q2.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_female, 2, function(x) sd(x))
ZWf.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5) %>% 
  mutate(sex = "ZWf")

# ZZf
a_ZZ_female <- a_female + post_pog_m2$b_sexZZf 
b_ZZ_female <- b_mass_female + post_pog_m2$`b_sexZZf:zlogMass`
y_ZZ_female <- a_ZZ_female + b_ZZ_female %*% t(zlogMass)
colnames(y_ZZ_female) <- zlogMass
Estimate <- colMeans(y_ZZ_female)
Q2.5 <- apply(y_ZZ_female, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_ZZ_female, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_ZZ_female, 2, function(x) sd(x))
ZZf.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "ZZf")

# ZZm
a_ZZ_male <- a_female + post_pog_m2$b_sexZZm 
b_ZZ_male <- b_mass_female + post_pog_m2$`b_sexZZm:zlogMass`
y_ZZ_male <- a_ZZ_male + b_ZZ_male %*% t(zlogMass)
colnames(y_ZZ_male) <- zlogMass
Estimate <- colMeans(y_ZZ_male)
Q2.5 <- apply(y_ZZ_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_ZZ_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_ZZ_male, 2, function(x) sd(x))
ZZm.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "ZZm")
# combined predictions
pog.mannual.pred <- rbind(ZZm.man.pred, ZZf.man.pred, ZWf.man.pred)


############## 
# df for summarizing raw datapoints for ggplot 
############## 
# add in the fitted value from the model.
pog.pre<- data.frame(predict(Pog_m2_brms))
pogona.data2<- cbind(pogona.data, pog.pre)
#Summarise 
pogona.raw.summary <- pogona.data %>% 
  group_by(day, id, sex) %>% 
  summarise(MR = mean(log(O2_min)),
            MR.se = std.error(log(O2_min)), 
            zlogMass = mean(zlogMass),
            zstartmass = mean(zstartmass), 
            zendmass = mean(zendmass)) %>% 
  mutate(a = ((zstartmass - zlogMass)^2),
         b = ((zstartmass - zlogMass)^2), 
         c = (a+b),
         d = c/2,
         sd = sqrt(d),
         mass.se = sd/(sqrt(2)))%>% 
  dplyr::select(-a,-b,-c,-d)


#############
# test for differences in body mass
#############
pog.bodymass <- brm(zlogMass ~ sex , data = pogona.raw.summary)
saveRDS(pog.bodymass, "./models/Pog_bodymass")
summary(pog.bodymass)
bm_diff <- hypothesis(pog.bodymass, 'sexZZf + sexZZm = 0')
mean(bm_diff$samples$H1)
HPDinterval(mcmc(bm_diff$samples$H1))
# p-value equivalent for the full hypothesis test. 
pMCMC <- 1 - (table(bm_diff$samples$H1 > 0)[1] / (length(bm_diff$samples$H1) - 1))   
pMCMC

#############
# Calculating 1SD +/1 mean
#############
Pog_SD_values <- pogona.data2 %>% 
  ungroup() %>% 
  summarise(mean = mean(zlogMass),
            sd =  sd(zlogMass),
            SD1.5_above = mean + (sd*1.5),
            SD1.5_below = mean - (sd*1.5)) %>% 
  mutate(across(1:4, round, 2))
############
# ZW female SD 
############
a_ZWfemale_SD <- post_pog_m2$b_Intercept
b_ZW_mass_female_SD <- post_pog_m2$b_zlogMass
zlogMass = c(-0.45, 0, 0.45)
y_ZWfemale_SD <- a_ZWfemale_SD + b_ZW_mass_female_SD %*% t(zlogMass)
ZWf.SD <- data.frame(y_ZWfemale_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "ZW_Female") 
############
# ZZ female SD 
############
a_ZZ_female_SD <- a_ZWfemale_SD + post_pog_m2$b_sexZZf 
b_ZZ_female_SD <- b_ZW_mass_female_SD + post_pog_m2$`b_sexZZf:zlogMass`
y_ZZ_female_SD <- a_ZZ_female_SD + b_ZZ_female_SD %*% t(zlogMass)
ZZf.SD <- data.frame(y_ZZ_female_SD) %>% 
  dplyr::rename("-1.5 SD" = 1,
                "Mean" = 2, 
                "+1.5 SD" = 3 )%>% 
  mutate(sex = "ZZ_Female")
############
# ZZ male SD 
############
a_ZZ_male_SD <- a_ZWfemale_SD + post_pog_m2$b_sexZZm 
b_ZZ_male_SD <- b_ZW_mass_female_SD + post_pog_m2$`b_sexZZm:zlogMass`
y_ZZ_male_SD <- a_ZZ_male_SD + b_ZZ_male_SD %*% t(zlogMass)
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


####################
# All Pogona Plots
####################
# color
mycolors <- c("#333333", "#990000", "#3399FF")
pog.reg <-  ggplot(data =pogona.data2 , aes(zlogMass, Estimate, group = sex, color= sex ))+
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate+Est.Error, colour = sex)) +
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate-Est.Error, colour = sex)) +
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate))+
  geom_vline(xintercept = 0.45, colour="black", linetype = "dotdash") + # Adding SD line
  geom_text(aes(x=0.45, label="+1.5 SD\n", y=-0.5), colour="black", angle=90, text=element_text(size=10))+
  geom_vline(xintercept = -0.45, colour="black", linetype = "dotdash") + # adding SD line
  geom_text(aes(x= -0.45, label="-1.5 SD\n", y=-0.5), colour="black", angle=90, text=element_text(size=10))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  scale_y_continuous(breaks=c(0,-0.5, -1, -1.5, -2, -2.5, -3, -3.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=10)) +
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
  theme(axis.text = element_text(size=10))
# combined figure
pogona.final.fig <- cowplot::plot_grid(pog.reg.plot,pog.sd.plot, labels = c('C', 'D'), ncol=2)

#combined pogona and bassiana figure:
grid.arrange(bassiana.final.fig, pogona.final.fig, nrow = 2)


####################
# growth ANCOVA comparison and posthoc test with emmeans
####################
growth <- read.csv(file = "final.analysis.data/growth.bassiana.pogona.csv") %>% 
  rename(sex = Geno.pheno)
# bassiana growth
bassiana.growth <- growth %>% 
  filter(Species == "Bassiana")
# 1) SVL
Bass.svl.mod <- lm(Growth.rate.SVL. ~sex +  SVL.1.mm, data = bassiana.growth)
Bass.svl.pc <- emmeans(Bass.svl.mod, pairwise ~ sex ) 
#2) MASS
Bass.mass.mod <- lm(Growth.rate.mass. ~sex + Mass.1.g, data = bassiana.growth)
Bass.mass.pc <- emmeans(Bass.mass.mod, pairwise ~ sex )

# pogona grwoth
##### 
pogona.growth <- growth %>% 
  filter(Species == "Pogona")
# 1) SVL
Pog.svl.mod <- lm(Growth.rate.SVL. ~ sex +  SVL.1.mm, data = pogona.growth)
Pog.svl.pairwise.comparision <- emmeans(Pog.svl.mod, pairwise ~ sex )
# 2) Mass
Pog.mass.mod <- lm(Growth.rate.mass. ~sex + Mass.1.g, data = pogona.growth)
Pog.mass.pairwise.comparision <- emmeans(Pog.mass.mod, pairwise ~ sex )

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
fisher.test(bassiana.test)
# analysis - Pogona
pogona.table <- table(pogona.growth$Dead.Y.N., pogona.growth$sex)
fisher.test(pogona.table)

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




############ 
# OLD DATA CHUNCKS - Not used
############ 

#############
##### SD prediction data and Analysis for SD plots
#############
# 1) Mean predictions
Bass_mean_pred <- Bass_SD_dat %>% 
  dplyr::select(c("mean", "sex"))
# brms test for predicted mean values
Bass_mean_Mod <- brm(mean ~ sex , data = Bass_mean_pred)
## extracting posteriors for just mean mass by sex
post.bass.mean.mass <-  posterior_samples(Bass_mean_Mod, pars = "^b")
XXf.posterior.mean.mass <- as.array(post.bass.mean.mass[,"b_Intercept"])
XXm.posterior.mean.mass <- as.array(post.bass.mean.mass[,"b_sexXXm"] + XXf.posterior.mean.mass ) 
XYm.posterior.mean.mass <- as.array(post.bass.mean.mass[,"b_sexXYm"] + XXf.posterior.mean.mass)

# H1: Like Phenotype Hypothesis
Diff.Pheno.mean.mass <- XYm.posterior.mean.mass - XXm.posterior.mean.mass
rope(Diff.Pheno.mean.mass, ci = 0.95)
plot(p_direction(Diff.Pheno.mean.mass))
pd.Pheno.mean.mass <- p_direction(Diff.Pheno.mean.mass)
pd_to_p(pd.Pheno.mean.mass, direction = "two-sided")
# H2: Like Genotype Hypothesis
Diff.Geno.mean.mass <- XXf.posterior.mean.mass - XXm.posterior.mean.mass
rope(Diff.Geno.mean.mass, ci = 0.95)
plot(p_direction(Diff.Geno.mean.mass))
pd.geno.mean.mass <- p_direction(Diff.Geno.mean.mass)
pd_to_p(pd.geno.mean.mass, direction = "two-sided")

#############
# 2) + 1.5 SD predictions
Bass_upper_pred <- Bass_SD_dat %>% 
  dplyr::select(c("upper", "sex"))
# brms test 
Bass_upper_Mod <- brm(upper ~ sex , data = Bass_upper_pred)
### Hypothesis testing 
## extracting posteriors for just mean mass by sex
post.bass.upper.mass <-  posterior_samples(Bass_upper_Mod, pars = "^b")
XXf.posterior.upper.mass <- as.array(post.bass.upper.mass[,"b_Intercept"])
XXm.posterior.upper.mass <- as.array(post.bass.upper.mass[,"b_sexXXm"] + XXf.posterior.upper.mass ) 
XYm.posterior.upper.mass <- as.array(post.bass.upper.mass[,"b_sexXYm"] + XXf.posterior.upper.mass)
# H1: Like Phenotype Hypothesis
Diff.Pheno.upper.mass <- XYm.posterior.upper.mass - XXm.posterior.upper.mass
rope(Diff.Pheno.upper.mass, ci = 0.95)
plot(p_direction(Diff.Pheno.upper.mass))
pd.Pheno.upper.mass <- p_direction(Diff.Pheno.upper.mass)
pd_to_p(pd.Pheno.upper.mass, direction = "two-sided")