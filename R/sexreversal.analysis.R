################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 

# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "patchwork", "ggExtra")

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
##############
## Model 1
##############
rerun1=FALSE
if(rerun1){
  
  suppressMessages(
    Bas_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 2000, warmup = 1000, 
                       thin = 1, control = list(adapt_delta=0.95), cores = 4))
  Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo", "waic"))
  saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
} else {Bas_m1_brms <- readRDS("./models/Bas_m1_brms")}
summary(Bas_m1_brms)

#####
## residual function
#####
residuals_brms <- function(model, data){
  fitted <- predict(model)[1:dim(data)[1], 1] # Predcit mean for each data
  e <- with(data, log(O2_min)) - fitted
  return(e)
}

# check residuals
e <- residuals_brms(Bas_m1_brms, bassiana.data)
hist(e)

# extract posteriors
post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")

post_meanXX_male_slope <- post_Bas_m1[,"b_zlogMass"] + 
                          post_Bas_m1[,"b_sexXX_Male:zlogMass"]
post_meanXY_male_slope <- post_Bas_m1[,"b_zlogMass"] + 
                          post_Bas_m1[,"b_sexXY_Male:zlogMass"]

# contarst slopes
contrastSlope <- as.mcmc(post_meanXX_male_slope - post_meanXY_male_slope)
mean(contrastSlope)
HPDinterval(contrastSlope)
#R2 of full model
bayes_R2(Bas_m1_brms)
# Model checks
plot(post_Bas_m1)

##############
## Model 2
##############
rerun2=FALSE
if(rerun2){
  mod_bas <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                sigma ~ zlogMass + ztime)
  Bas_m2_brms <- brm(mod_bas, family = gaussian(), 
                     data = bassiana.data, iter= 2000, warmup = 1000, 
                     thin = 1, control = list(adapt_delta=0.95), cores = 4)
  Bas_m2_brms <- add_criterion(Bas_m2_brms, c("loo", "waic"))
  saveRDS(Bas_m2_brms, "./models/Bas_m2_brms")
} else {
  # read file in
  Bas_m2_brms <- readRDS(file = "models/Bas_m2_brms")
}
plot(Bas_m2_brms)
summary(Bas_m2_brms)

# Check residuals. 
plot(Bas_m2_brms)
e <- residuals_brms(Bas_m2_brms, bassiana.data)
hist(e)

#R2 of full model
bayes_R2(Bas_m2_brms)


####################
# Model comparison
####################
# SE is large and difference is probably not sufficient to warrant a het model
loo_compare(Bas_m1_brms, Bas_m2_brms)


####################  
# extract posteriors for best model + Plotting 
####################
post_Bas_m2 <- posterior_samples(Bas_m2_brms, pars = "^b")
dimnames(post_Bas_m2)

# extracting posteriors 
XXf <- as.array(post_Bas_m2[,"b_zlogMass"])
XXm <- as.array(post_Bas_m2[,"b_zlogMass"] + 
                post_Bas_m2[,"b_sexXX_Male:zlogMass"])
XYm <- as.array(post_Bas_m2[,"b_zlogMass"] + 
                post_Bas_m2[,"b_sexXY_Male:zlogMass"])

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

# contrast phenotype slopes
bass.phenotype <- as.mcmc(XXm - XYm)
mean(bass.phenotype)
HPDinterval(bass.phenotype)

# contrast genotype slopes 
bass.genotype <- as.mcmc(XXf - XXm)
mean(bass.genotype)
HPDinterval(bass.genotype)


####################
# Sex predictions from Bas_m2_brms for regression ggplot
####################
#XY female 
newdata_XY_female <- data.frame( sex = "XX_Female",
                                 zlogMass = seq(-0.3856391, 0.3710889, 
                                                length.out = 100),
                                 ztime = 0, day = NA, id = NA) # If I understand this correctly, NA for day and id would choose an average deviation, so 0. This should allow you to use the re_formula in model predictions

prXX_female <- data.frame(cbind(predict(Bas_m2_brms, 
                                        newdata = newdata_XY_female, 
                                        summary = TRUE),
                                zlogMass = newdata_XY_female$zlogMass))
prXX_female <- prXX_female %>% 
  mutate(sex = "XX_Female") 
#XY male 
newdata_XY_male <- data.frame(sex = "XY_Male",
                              zlogMass = seq(-0.3856391, 0.3710889, 
                                             length.out = 100),
                              ztime = 0, day = NA, id = NA) 
prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata_XY_male,
                                      summary = TRUE), 
                              zlogMass = newdata_XY_male$zlogMass))
prXY_male <- prXY_male %>% 
              mutate(sex = "XY_Male")
#XX male
newdata_XX_male <- data.frame(sex = "XX_Male",
                              zlogMass = seq(-0.3856391, 0.3710889, 
                                             length.out = 100),
                              ztime = 0, day = NA, id = NA) 
prXX_male <- data.frame(cbind(predict(Bas_m2_brms, 
                                      newdata = newdata_XX_male), 
                                      zlogMass=newdata_XX_male$zlogMass))
prXX_male <- prXX_male %>% 
              mutate(sex = "XX_Male")

# setting up prediction data into one df 
bass.mod.dat <- rbind(prXX_female, prXX_male, prXY_male) %>% 
  group_by(sex)    

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
# test for differences in body mass
#############
bodymass <- brm(zlogMass ~ sex , data = bass.raw.summary)
summary(bodymass)
bm_diff <- hypothesis(bodymass, 'sexXX_Male + sexXY_Male = 0')
mean(bm_diff$samples$H1)
HPDinterval(mcmc(bm_diff$samples$H1))

# p-value equivalent for the full hypothesis test. 
pMCMC <- 1 - (table(bm_diff$samples$H1 > 0)[1] / (length(bm_diff$samples$H1) - 1))     


#############
# # Calculating 1SD above & below predicted values for density plots
#############
# 1sd and mean zlogmass values for seq
SD_2_values <- bassiana.data %>% 
  group_by(sex) %>% 
  summarise(mean = mean(zlogMass),
            sd =  sd(zlogMass),
            above = mean + sd,
            below = mean - sd)
############
#1SD above the mean
############
# XX females
newdata_XX_Female <- data.frame(sex = "XX_Female",
                      zlogMass = seq(0.212, 0.3710889, length.out = 100),
                      ztime = 0, day = NA, id = NA)

prXX_female <- data.frame(cbind(predict(Bas_m2_brms, 
                                        newdata = newdata_XX_Female, 
                                        summary = TRUE), 
                                zlogMass=newdata_XX_Female$zlogMass))

XX_female_above <- prXX_female %>% 
                    mutate(sex = "XX_Female",
                           se = std.error(Estimate))

## XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(0.165, 0.3710889, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))

XY_male_above <- prXY_male %>% 
  mutate(sex = "XY_Male",
         se = std.error(Estimate))
# XX male
newdata <- data.frame(sex = "XX_Male", 
                      zlogMass = seq(0.166, 0.3710889, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass = newdata$zlogMass))
XX_male_above <- prXX_male %>% 
  mutate(sex = "XX_Male",
         se = std.error(Estimate))

SD.bass.above <- rbind(XX_male_above, XY_male_above, XX_female_above) %>% 
  mutate(test = "+1SD")

############
#1SD below the mean
############
# XX female
newdata <- data.frame(sex = "XX_Female",
                      zlogMass = seq(-0.3856391, 0.171, length.out = 100),
                      ztime = 0)

prXX_female <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata,
                                        re_formula = NA, summary = TRUE), 
                                zlogMass=newdata$zlogMass))

XX_female_below <- prXX_female %>% 
  mutate(sex = "XX_Female", 
         se = std.error(Estimate))
## XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(-0.3856391, 0.156, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))
XY_male_below <- prXY_male %>% 
  mutate(sex = "XY_Male",
         se = std.error(Estimate))
# XX male
newdata <- data.frame(sex = "XX_Male", 
                      zlogMass = seq(-0.3856391, 0.228, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))

XX_male_below <- prXX_male %>% 
  mutate(sex = "XX_Male",
         se = std.error(Estimate))

SD.bass.below <- rbind(XX_male_below, XY_male_below, XX_female_below) %>% 
  mutate(test = "-1SD")

############
# MEAN 
# @dan found a mistake here and changed the seq values to reflect the right ranges. This error explains the wonky density plots
############
# XX female
newdata <- data.frame(sex = "XX_Female",
                      zlogMass = seq(-0.171, 0.212, length.out = 100),
                      ztime = 0)

prXX_female <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                        re_formula = NA, summary = TRUE), 
                                zlogMass=newdata$zlogMass))

XX_female_within <- prXX_female %>% 
  mutate(sex = "XX_Female", 
         se = std.error(Estimate))
# XY male 
newdata <- data.frame(sex = "XY_Male",
                      zlogMass = seq(-0.156, 0.165, length.out = 100),
                      ztime = 0)

prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE),
                              zlogMass=newdata$zlogMass))
XY_male_within <- prXY_male %>% 
  mutate(sex = "XY_Male",
         se = std.error(Estimate))
# XX male
newdata <- data.frame(sex = "XX_Male",
                      zlogMass = seq(-0.228, 0.166, length.out = 100),
                      ztime = 0)

prXX_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, 
                                      re_formula = NA, summary = TRUE), 
                              zlogMass=newdata$zlogMass))
XX_male_within <- prXX_male %>% 
  mutate(sex = "XX_Male",
         se = std.error(Estimate))

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
  geom_smooth(data = bass.mod.dat, aes(x=zlogMass, y=Estimate, colour = sex)) + 
  geom_ribbon(data = bass.mod.dat, aes(x=zlogMass, y=Estimate, ymin = Estimate-Est.Error, ymax = Estimate+Est.Error, fill = sex, colour = sex), alpha = 0.2) +
  geom_smooth(data = bass.mod.dat, aes(x=zlogMass, y=Estimate+Est.Error, colour = sex)) +
  geom_smooth(data = bass.mod.dat, aes(x=zlogMass, y=Estimate-Est.Error, colour = sex)) + # @Kris, to smooth just add in two more smoothed lines ontop
  geom_smooth(data = bass.mod.dat, aes(x=zlogMass, y=Estimate))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
ggExtra::ggMarginal(p, margins = "x", groupColour = TRUE, groupFill = TRUE)
### save plot ##
ggsave(filename ="figures/bassiana.mod2.regression.pdf",   height = 10, width = 16)

############
# density plot of regression +- SD of mean for bassiana
############
# combining data for plots
SD.bass.mod.dat <- rbind(SD.bass.above, SD.bass.below, SD.bass.within) %>% 
  group_by(test, sex)
# SD Plot
SD.bass.mod.dat$test <- factor(SD.bass.mod.dat$test, levels = c("+1SD", "Mean", "-1SD"))
mycolors <- c("#333333", "#990000", "#3399FF")
ggplot(SD.bass.mod.dat, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., scales = "free", switch="y")+
  scale_y_continuous(position = "right")+
  xlab("Predicted Mean Metabolic Rate")+
  scale_x_continuous(name="Predicted Mean Metabolic Rate", breaks = seq (-5.3, -3.9, by=0.3), limits=c(-5.3, -3.9))+
  theme_bw()

# save plot
ggsave(filename ="figures/bassiana.mod2.density.plot.pdf",   height = 10, width = 16)



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
## Model 1 @dan control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 4000
##############
Pog_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                   family = "gaussian", data = pogona.data, iter= 4000, warmup = 1000, 
                   thin = 5, cores = 4)
# @dan got a warning about our waic values and it suggested to add moment_match = TRUE 
Pog_m1_brms <- add_criterion(Pog_m1_brms, c("loo", "waic"), moment_match = TRUE)
plot(Pog_m1_brms)
summary(Pog_m1_brms)
bayes_R2(Pog_m1_brms)
saveRDS(Pog_m1_brms, "./models/Pog_m1_brms")

# check residuals
e <- residuals_brms(Pog_m1_brms, pogona.data)
hist(e)


##############
## Model 2 @ dan same as above, control = list(adapt_delta=0.95) in our model is slowing up my computer so I removed it. Also based off ESS values and ACF plots I changed the thin to 5 and the iterations to 4000
##############
mod <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
          sigma ~ zlogMass + ztime)
Pog_m2_brms <- brm(mod, family = "gaussian", data = pogona.data, iter= 4000, warmup = 1000, 
                   thin = 5, cores = 4)
# @dan same as above, I got a warning about our waic values and it suggested to add moment_match = TRUE 
Pog_m2_brms <- add_criterion(Pog_m2_brms, c("loo", "waic"), moment_match = TRUE)
saveRDS(Pog_m2_brms, "./models/Pog_m2_brms")

# checking lags for this model, looks much better after increasing thinning and iterations  
draws <- as.array(Pog_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept"), lags =10)

# Model checks
plot(Pog_m2_brms)
summary(Pog_m2_brms)

#R2 of full model
bayes_R2(Pog_m2_brms)

####################
# Model comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)

####################
# Checking residuals
####################
e <- residuals_brms(Pog_m2_brms, pogona.data)
hist(e)


####################  
# extract posteriors for best model + Plotting 
####################
post_pog_m2 <- posterior_samples(Pog_m2_brms, pars = "^b")
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
ggsave(filename ="figures/Pogona.mod2.posterior.pdf",  height = 5, width = 7)

# contrast phenotype slopes
pog.phenotype <- as.mcmc(ZWf - ZZf)
mean(pog.phenotype)
HPDinterval(pog.phenotype)

# contrast genotype slopes 
pog.genotype <- as.mcmc(ZZf - ZZm)
mean(pog.genotype)
HPDinterval(pog.genotype)

############################### FINISHED HERE KW 15 Sep 2021 Need to continue working 
####################
# Predictions from best model for plots
####################
#ZW female
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(-0.5029249, 1.23776, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZWf <- prZWf %>% 
  mutate(sex = "ZWf",
         se = std.error(Estimate))
#ZZ male 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.5029249, 1.23776, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZm <- prZZm %>% 
  mutate(sex = "ZZm", 
         se = std.error(Estimate))
#ZZ females
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.5029249, 1.23776, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZf <- prZZf %>% 
  mutate(sex = "ZZf", 
         se = std.error(Estimate))
# setting up prediction data into one df
mod.pog.dat <- rbind(prZWf, prZZm, prZZf) %>% 
  group_by(sex)

############## ############### ############## ##### 
#### df for summarizing raw datapoints ggplot #####
############## ############### ############## ####
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
         mass.se = sd/(sqrt(2)))

###########
# quick check for differences in body mass
###########
bodymass <- lm(zlogMass ~sex  , data = pogona.raw.summary)
summary(bodymass)
anova(bodymass)

######## ######## ######## ########  ########
######          Pogona    ########  ########           
######### Calculating 1SD above ########  #####
#########  & below predicted values ########  
######## ######## ######## ######## ########
# 1sd from mean plots
Pog_SD_2_values <- pogona.data  %>% 
  summarise(mean = mean(zlogMass),
            sd = sd(zlogMass),
            above = mean +sd,
            below = mean -sd) %>% 
  as.data.frame()
############
#1SD above the mean
############
# ZWf
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(0.331, 1.23776, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZWf_above <- prZWf %>% 
  mutate(sex = "ZWf", 
         se = std.error(Estimate))
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(0.275, 1.23776, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_above <- prZZm %>% 
  mutate(sex = "ZZm",
         se = std.error(Estimate))
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(0.305, 1.23776, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_above <- prZZf %>% 
  mutate(sex = "ZZf",
         se = std.error(Estimate))
SD.pog.above <- rbind(ZWf_above, ZZm_above, ZZf_above) %>% 
  mutate(test = "+1SD")
############
#1SD below the mean
############
# ZWf
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(-0.5029249, -0.331, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZWf_below <- prZWf %>% 
  mutate(sex = "ZWf", 
         se = std.error(Estimate))
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.5029249, -0.275, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_below <- prZZm %>% 
  mutate(sex = "ZZm",
         se = std.error(Estimate))
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.5029249, -0.305, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_below <- prZZf %>% 
  mutate(sex = "ZZf",
         se = std.error(Estimate))
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
  mutate(sex = "ZWf", 
         se = std.error(Estimate))
## ZZm 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.275, 0.275, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZm_within <- prZZm %>% 
  mutate(sex = "ZZm",
         se = std.error(Estimate))
# ZZf
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.305, 0.305, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
ZZf_within <- prZZf %>% 
  mutate(sex = "ZZf",
         se = std.error(Estimate))
SD.pog.within <- rbind(ZWf_within, ZZf_within, ZZm_within) %>% 
  mutate(test = "Mean")

####################
# All Pogona Plots
####################
mycolors <- c("#333333", "#990000", "#3399FF")
p <- ggplot(data = pogona.raw.summary, aes(zlogMass, MR, group = sex, color= sex )) +
  geom_point(alpha =.6)+
  geom_errorbar(aes(ymin = MR-MR.se, ymax = MR+MR.se)) + 
  geom_errorbarh(aes(xmin = zlogMass-mass.se, xmax = zlogMass+mass.se))+
  geom_ribbon(data = mod.pog.dat, aes(x=zlogMass, y=Estimate, ymin = (Estimate-se), ymax = (Estimate+se), fill = sex, colour = sex), stat = "identity", alpha = .5)+
  geom_smooth(data = mod.pog.dat, aes(x=zlogMass, y=Estimate))+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
reg <- ggExtra::ggMarginal(p, margins = "x", groupColour = TRUE, groupFill = TRUE)
reg

############
# dinsity plot pogona
############
# combinding data for plots
SD.pog.mod.dat <- rbind(SD.pog.above, SD.pog.below, SD.pog.within) %>% 
  group_by(test, sex)
SD.pog.mod.dat$test <- factor(SD.pog.mod.dat$test, levels = c("+1SD", "Mean", "-1SD"))
mycolors <- c("#333333", "#990000", "#3399FF")
ggplot(SD.pog.mod.dat, aes(x=Estimate, group = sex, fill = sex)) +
  geom_density(alpha = .4) +
  scale_fill_manual(values = mycolors, guide = FALSE)+
  facet_grid(test~., switch="y")+
  scale_y_continuous(position = "right")+
  scale_x_continuous(name="Predicted Mean Metabolic Rate", breaks = seq (-3,0, by=0.5), limits = c(-3,0))+
  theme_bw()

# save plot
ggsave(filename ="figures/pog.mod2.density.plot.pdf",   height = 10, width = 16)
