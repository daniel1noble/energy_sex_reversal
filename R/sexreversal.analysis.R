################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 

# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "patchwork", "ggExtra", "gridExtra", "cowplot")

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
## Model 1
# @ dan updated the model so that ESS values were near 3000 so increased thinning to 5 and upped the iterations
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Bas_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 5000, warmup = 1000, 
                       thin = 5, cores = 4))
  Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo", "waic"))
  saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
} else {Bas_m1_brms <- readRDS("./models/Bas_m1_brms")}
summary(Bas_m1_brms)

# checking lags for this model 
# @dan looks  better after increasing thinning and iterations  
draws <- as.array(Bas_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXX_Male", "b_sexXY_Male", "b_zlogMass"), lags =10)

# check residuals
e <- residuals_brms(Bas_m1_brms, bassiana.data)
hist(e)
#R2 of full model
bayes_R2(Bas_m1_brms)


##############
## Model 2
# @ dan updated the model so that ESS values were near 3000 so increased thinning to 5 and upped the iterations
##############
rerun2=FALSE
if(rerun2){
  mod_bas <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                sigma ~ zlogMass + ztime)
  Bas_m2_brms <- brm(mod_bas, family = "gaussian", 
                     data = bassiana.data, iter= 5000, warmup = 1000, 
                     thin = 5, cores = 4)
  Bas_m2_brms <- add_criterion(Bas_m2_brms, c("loo", "waic"))
  saveRDS(Bas_m2_brms, "./models/Bas_m2_brms")
} else {
  # read file in
  Bas_m2_brms <- readRDS(file = "models/Bas_m2_brms")
}
plot(Bas_m2_brms)
summary(Bas_m2_brms)

# checking lags for this model
# @ dan looks much better after increasing thinning and iterations  
draws <- as.array(Bas_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexXX_Male", "b_sexXY_Male", "b_zlogMass"), lags =10)

# Check residuals
e <- residuals_brms(Bas_m2_brms, bassiana.data)
hist(e)
#R2 of full model
bayes_R2(Bas_m2_brms)


####################
# Model comparison
####################
# SE is large and difference is probably not sufficient to warrant a het model
# going with Bas_m1_brms
loo_compare(Bas_m1_brms, Bas_m2_brms)
tab_model(Bas_m1_brms, file = "R/Figs/Bas_m1_brms_table.html")

####################  
# extract posteriors for Bas_m1_brms model + Plotting 
####################
post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")
dimnames(Bas_m1_brms)

# extracting posteriors 
XXf <- as.array(post_Bas_m1[,"b_zlogMass"])
XXm <- as.array(post_Bas_m1[,"b_zlogMass"] + 
                post_Bas_m1[,"b_sexXX_Male:zlogMass"])
XYm <- as.array(post_Bas_m1[,"b_zlogMass"] + 
                post_Bas_m1[,"b_sexXY_Male:zlogMass"])
bass.dat <- cbind(XXf, XXm, XYm)

# plotting posteriors lot
mcmc_areas(bass.dat, 
           pars = c("XXf", "XXm", "XYm"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 
ggsave(filename ="figures/bassiana.mod2.posterior.pdf",  height = 5, width = 7)


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
  mutate(sex = "XX_Female")

# XXm
a_XX_male_xx <- a_female + post_Bas_m1$b_sexXX_Male 
b_XX_male_xx <- b_mass_female + post_Bas_m1$`b_sexXX_Male:zlogMass`
y_XX_male <- a_XX_male_xx + b_XX_male_xx %*% t(zlogMass)
colnames(y_XX_male) <- zlogMass
Estimate <- colMeans(y_XX_male)
Q2.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XX_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XX_male, 2, function(x) sd(x))
XXm.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XX_Male")

# XYm
a_XY_male <- a_female + post_Bas_m1$b_sexXY_Male 
b_XY_male <- b_mass_female + post_Bas_m1$`b_sexXY_Male:zlogMass`
y_XY_male <- a_XY_male + b_XY_male %*% t(zlogMass)
colnames(y_XY_male) <- zlogMass
Estimate <- colMeans(y_XY_male)
Q2.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.025))
Q97.5 <- apply(y_XY_male, 2, function(x) quantile(x, probs = 0.975))
Est.Error <- apply(y_XY_male, 2, function(x) sd(x))
XYm.man.pred <- data.frame(zlogMass, Estimate, Est.Error, Q2.5, Q97.5)%>% 
  mutate(sex = "XY_Male")
# combine predictions
bass.mannual.pred <- rbind(XXf.man.pred, XXm.man.pred, XYm.man.pred)


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
bm_diff <- hypothesis(bass.bodymass, 'sexXX_Male + sexXY_Male = 0')
mean(bm_diff$samples$H1)
HPDinterval(mcmc(bm_diff$samples$H1))
# p-value equivalent for the full hypothesis test. 
pMCMC <- 1 - (table(bm_diff$samples$H1 > 0)[1] / (length(bm_diff$samples$H1) - 1))     


#############
# # Calculating 1SD above & below predicted values for density plots
#############
# 1sd and mean zlogmass values for seq
# group by sex if this is what we want SD for each sex
SD_2_values <- bassiana.data %>% 
  mutate(max_mass = max(zlogMass),
         min_mass =  min(zlogMass)) %>% 
  summarise(mean = mean(zlogMass),
            sd =  sd(zlogMass),
            SD1_above = mean + sd,
            SD1_below = mean - sd,
            max_mass = unique(max_mass),
            min_mass = unique(min_mass))

############
#1SD above the mean
############
# XX females
newdata_XX_Female <- data.frame(sex = "XX_Female",
                      zlogMass = seq(0.1855218, 0.3710889, length.out = 100),
                      ztime = 0, day = NA, id = NA)

prXX_female <- data.frame(cbind(predict(Bas_m1_brms, 
                                        newdata = newdata_XX_Female, 
                                        summary = TRUE), 
                                zlogMass=newdata_XX_Female$zlogMass))

XX_female_above <- prXX_female %>% 
                    mutate(sex = "XX_Female")

## XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(0.1855218, 0.3710889, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))

XY_male_above <- prXY_male %>% 
  mutate(sex = "XY_Male")

# XX male
newdata <- data.frame(sex = "XX_Male", 
                      zlogMass = seq(0.1855218, 0.3710889, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass = newdata$zlogMass))
XX_male_above <- prXX_male %>% 
  mutate(sex = "XX_Male")

SD.bass.above <- rbind(XX_male_above, XY_male_above, XX_female_above) %>% 
  mutate(test = "+1SD")

############
#1SD below the mean
############
# XX female
newdata <- data.frame(sex = "XX_Female",
                      zlogMass = seq(-0.3856391, -0.1855218, length.out = 100),
                      ztime = 0)

prXX_female <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata,
                                        re_formula = NA, summary = TRUE), 
                                zlogMass=newdata$zlogMass))

XX_female_below <- prXX_female %>% 
  mutate(sex = "XX_Female")
## XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(-0.3856391, -0.1855218, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))
XY_male_below <- prXY_male %>% 
  mutate(sex = "XY_Male")
# XX male
newdata <- data.frame(sex = "XX_Male", 
                      zlogMass = seq(-0.3856391, -0.1855218, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))

XX_male_below <- prXX_male %>% 
  mutate(sex = "XX_Male")

SD.bass.below <- rbind(XX_male_below, XY_male_below, XX_female_below) %>% 
  mutate(test = "-1SD")

############
# MEAN 
############
# XX female
newdata <- data.frame(sex = "XX_Female",
                      zlogMass = seq(-0.1855218, 0.1855218, length.out = 100),
                      ztime = 0)

prXX_female <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                        re_formula = NA, summary = TRUE), 
                                zlogMass=newdata$zlogMass))

XX_female_within <- prXX_female %>% 
  mutate(sex = "XX_Female")
# XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(-0.1855218, 0.1855218, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))
XY_male_within <- prXY_male %>% 
  mutate(sex = "XY_Male")
# XX male
newdata <- data.frame(sex = "XX_Male",
                      zlogMass = seq(-0.1855218, 0.1855218, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m1_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE), 
                              zlogMass=newdata$zlogMass))
XX_male_within <- prXX_male %>% 
  mutate(sex = "XX_Male")

SD.bass.within <- rbind(XX_male_within, XY_male_within, XX_female_within) %>% 
  mutate(test = "Mean")

####################
# ALL Bassiana Plots
#######################
# Regression plot with predicted line and body mass
mycolors <- c("#333333", "#990000", "#3399FF")

p<-  ggplot(data =bassiana.data2 , aes(zlogMass, Estimate, group = sex, color= sex )) +
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate+Est.Error, colour = sex)) +
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate-Est.Error, colour = sex)) + # @Kris, to smooth just add in two more smoothed lines ontop
  geom_smooth(data = bass.mannual.pred, aes(x=zlogMass, y=Estimate))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  labs(tag = "A")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
reg.plot <- ggMarginal(p, margins = "x", groupColour = TRUE, groupFill = TRUE)

############
# density plot of regression +- SD of mean for bassiana
############
# combining data for plots
SD.bass.mod.dat <- rbind(SD.bass.above, SD.bass.below, SD.bass.within) %>% 
  group_by(test, sex)
saveRDS(SD.bass.mod.dat, "final.analysis.data/Bassiana.SD.mod.dat.RDS")
# SD Plot
SD.bass.mod.dat$test <- factor(SD.bass.mod.dat$test, levels = c("+1SD", "Mean", "-1SD"))
mycolors <- c("#333333", "#990000", "#3399FF")
sd.plot <- ggplot(SD.bass.mod.dat, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., scales = "free", switch="y")+
  scale_y_continuous(position = "right")+
  xlab("Predicted Mean Metabolic Rate")+
  scale_x_continuous(name="Predicted Mean Metabolic Rate", breaks = seq (-5.3, -3.9, by=0.3), limits=c(-5.3, -3.9))+
  theme_bw()

# combind plots
bassiana.final.fig <- plot_grid(reg.plot,sd.plot, ncol = 2)

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
## Model 1 @dan control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 5000
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


# model checks
# checking lags for this model, 
#@dan looks much better after increasing thinning and iterations  
summary(Pog_m1_brms)
draws <- as.array(Pog_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept", "b_sexZZf", "b_sexZZm"), lags =10)
#R2 of full model
bayes_R2(Pog_m1_brms)

####################
# Checking residuals
####################
e <- residuals_brms(Pog_m1_brms, pogona.data)
hist(e)


##############
## Model 2 @ dan same as above, control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 5000
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

# Model checks
# checking lags for this model, looks much better after increasing thinning and iterations  
draws <- as.array(Pog_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_sexZZf", "b_sexZZm"), lags =10)
# other checks
plot(Pog_m2_brms)
summary(Pog_m2_brms)
#R2 of full model
bayes_R2(Pog_m2_brms)

####################
# Checking residuals
####################
e <- residuals_brms(Pog_m2_brms, pogona.data)
hist(e)

####################
# Model comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)
tab_model(Pog_m2_brms, file = "R/Figs/Pog_m2_brms_table.html")

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
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 


####################  
# manual predict values for lines
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
bm_diff <- hypothesis(bodymass, 'sexZZf + sexZZm = 0')
mean(bm_diff$samples$H1)
HPDinterval(mcmc(bm_diff$samples$H1))
# p-value equivalent for the full hypothesis test. 
pMCMC <- 1 - (table(bm_diff$samples$H1 > 0)[1] / (length(bm_diff$samples$H1) - 1))   


#############
## Calculating 1SD above & below predicted values for density plots
# @dan sd for zlogmass seems really weird, the errors are very small. Any ideas? 
#############
# 1sd and mean zlogmass values for seq
Pog_SD_2_values <- pogona.data %>% 
  ungroup() %>% 
  mutate(max_mass = max(zlogMass),
         min_mass =  min(zlogMass)) %>% 
  summarise(mean = mean(zlogMass),
            sd =  sd(zlogMass),
            SD1_above = mean + sd,
            SD1_below = mean - sd,
            max_mass = unique(max_mass),
            min_mass = unique(min_mass))

############
#1SD above the mean
############
# ZWf
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(0.331, 1.26, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZWf_above <- prZWf %>% 
  mutate(sex = "ZWf")
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(0.331, 1.26, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_above <- prZZm %>% 
  mutate(sex = "ZZm")
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(0.331, 1.26, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_above <- prZZf %>% 
  mutate(sex = "ZZf")
SD.pog.above <- rbind(ZWf_above, ZZm_above, ZZf_above) %>% 
  mutate(test = "+1SD")

############
#1SD below the mean
############
# ZWf
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq( -0.502, -0.301, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZWf_below <- prZWf %>% 
  mutate(sex = "ZWf")
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.502, -0.301, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_below <- prZZm %>% 
  mutate(sex = "ZZm")
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.502, -0.301, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_below <- prZZf %>% 
  mutate(sex = "ZZf")
SD.pog.below <- rbind(ZWf_below, ZZm_below, ZZf_below) %>% 
  mutate(test = "-1SD")
############
#within 1sd the mean
############
# ZWf
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(-0.331, 0.331, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZWf_within <- prZWf %>% 
  mutate(sex = "ZWf")
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.331, 0.331, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_within <- prZZm %>% 
  mutate(sex = "ZZm")
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.331, 0.331, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_within <- prZZf %>% 
  mutate(sex = "ZZf")
SD.pog.within <- rbind(ZWf_within, ZZf_within, ZZm_within) %>% 
  mutate(test = "Mean")

####################
# All Pogona Plots
####################
# color
mycolors <- c("#333333", "#990000", "#3399FF")
p<-  ggplot(data =pogona.data2 , aes(zlogMass, Estimate, group = sex, color= sex ))+
  # Add in the predicted data given each rows data. 
  geom_point(alpha =.3)+
  # Now add in the model predictions
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate+Est.Error, colour = sex)) +
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate-Est.Error, colour = sex)) + # @Kris, to smooth just add in two more smoothed lines ontop
  geom_smooth(data = pog.mannual.pred, aes(x=zlogMass, y=Estimate))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  labs(tag = "B")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
reg.plot<-ggMarginal(p, margins = "x", groupColour = TRUE, groupFill = TRUE)

############
# density plots 
# @dan yeah the SD plots here look pretty crap at +1SD
############
# combinding data for plots
SD.pog.mod.dat <- rbind(SD.pog.above, SD.pog.below, SD.pog.within) %>% 
  group_by(test, sex)
saveRDS(SD.pog.mod.dat, "final.analysis.data/PogonaSD.mod.dat.RDS")
SD.pog.mod.dat$test <- factor(SD.pog.mod.dat$test, levels = c("+1SD", "Mean", "-1SD"))
mycolors <- c("#333333", "#990000", "#3399FF")
sd.plot <- ggplot(SD.pog.mod.dat, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., switch="y", scales="free_y")+
  scale_y_continuous(position = "right")+
  scale_x_continuous(name="Predicted Mean Metabolic Rate", breaks = seq (-3,0, by=0.5), limits = c(-3,0))+
  theme_bw()

# combined figure
pogona.final.fig <- plot_grid(reg.plot,sd.plot, ncol = 2)

#combined pogona and bassiana figure:
grid.arrange(bassiana.final.fig, pogona.final.fig, nrow = 2)
