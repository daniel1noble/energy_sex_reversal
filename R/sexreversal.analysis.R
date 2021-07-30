################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 

# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa", "tidybayes", "bayesplot", "rstanarm", "plotrix")

#####################################
######### Bassiana  O2 data #########
#####################################
#########
#Load dataa
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
    suppressMessages(
        Bas_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                       family = "gaussian", data = bassiana.data, iter= 2000, warmup = 1000, thin = 1, control = list(adapt_delta=0.95), cores = 4))
      	Bas_m1_brms <- add_criterion(Bas_m1_brms, c("loo", "waic"))
    saveRDS(Bas_m1_brms, "./models/Bas_m1_brms")
    summary(Bas_m1_brms)

    # extract posteriors
    post_Bas_m1 <- posterior_samples(Bas_m1_brms, pars = "^b")
    post_meanXX_male_slope <- post_Bas_m1[,"b_zlogMass"] + post_Bas_m1[,"b_sexXX_Male:zlogMass"]
    post_meanXY_male_slope <- post_Bas_m1[,"b_zlogMass"] + post_Bas_m1[,"b_sexXY_Male:zlogMass"]

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
    mod_bas <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
              sigma ~ zlogMass + ztime)
    Bas_m2_brms <- brm(mod_bas, family = gaussian(), data = bassiana.data, iter= 2000, warmup = 1000, thin = 1, control = list(adapt_delta=0.95), cores = 4)
    Bas_m2_brms <- add_criterion(Bas_m2_brms, c("loo", "waic"))
    saveRDS(Bas_m2_brms, "./models/Bas_m2_brms")
    
    # read file in
    Bas_m2_brms <- readRDS(file = "models/Bas_m2_brms")

    # Model checks
    plot(Bas_m2_brms)
    summary(Bas_m2_brms)
    #R2 of full model
    bayes_R2(Bas_m2_brms)

####################
# Checking residuals from Bas_m2_brms
#################### 
resid.pred <- as.data.frame(predict(Bas_m2_brms))
ind.log.mr <-as.data.frame(log(bassiana.data$O2_min)) %>% 
  rename(log.mr = `log(bassiana.data$O2_min)`)

df <- bind_cols(ind.log.mr, resid.pred) %>% 
  mutate(diff = (log.mr - Estimate))
hist(df$diff)
ggplot(df, aes(x=diff)) + 
  geom_density()

# prediction of overall model using spread draws
Bas_m2_brms %>%
      spread_draws(b_Intercept, b_zlogMass) %>%
      mutate(zlogmass = list(seq(-0.4, 0.4, length.out = 100))) %>% #the observed value range of zlogmass
      unnest(zlogmass) %>%
      mutate(pred = exp(b_zlogMass + b_zlogMass*zlogmass)/(1+exp(b_zlogMass + b_zlogMass*zlogmass))) %>%
      group_by(zlogmass) %>%
      summarise(pred_m = mean(pred, na.rm = TRUE),
                pred_low = quantile(pred, prob = 0.025),
                pred_high = quantile(pred, prob = 0.975)) %>%
      ggplot(aes(x = zlogmass, y = pred_m)) +
      geom_line() +
      geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha=0.2)
ggsave(filename ="figures/bassiana.mod2.predicition.pdf",  height = 5, width = 7)

####################
# Model comparison
####################

    loo_compare(Bas_m1_brms, Bas_m2_brms)

####################  
# extract posteriors for best model + Plotting 
####################
post_Bas_m2 <- posterior_samples(Bas_m2_brms, pars = "^b")
dimnames(post_Bas_m2)

# extracting posteriors 
####### CHECK XX FEMALE
    XXf <- as.array(post_Bas_m2[,"b_zlogMass"])
    XXm <- as.array(post_Bas_m2[,"b_zlogMass"] + post_Bas_m2[,"b_sexXX_Male:zlogMass"])
    XYm <- as.array(post_Bas_m2[,"b_zlogMass"] + post_Bas_m2[,"b_sexXY_Male:zlogMass"])
   
    bass.dat <- cbind(XXf, XXm, XYm )
      
    # plotting posteriors lot
    mcmc_intervals(bass.dat, 
                   pars = c("XXf", "XXm", "XYm"), 
                   prob = 0.95, 
                   prob_outer = 0.99, 
                   point_est = "mean") 
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
# Predictions from best model for plots
####################
    #XY female
    newdata <- data.frame(
      sex = "XX_Female",
      zlogMass = seq(-0.3856391, 0.3710889, length.out = 100),
      ztime = 0)
    prXX_female <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXX_female <- prXX_female %>% 
      mutate(sex = "XX_Female")
    
    #XY male 
    newdata <- data.frame(
      sex = "XY_Male",
      zlogMass = seq(-0.3856391, 0.3710889, length.out = 100),
      ztime = 0)
    prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXY_male <- prXY_male %>% 
      mutate(sex = "XY_Male")
    
    #XX male
    newdata <- data.frame(
      sex = "XX_Male",
      zlogMass = seq(-0.3856391, 0.3710889, length.out = 100),
      ztime = 0)
    prXX_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXX_male <- prXX_male %>% 
      mutate(sex = "XX_Male")
    
    # setting up predicition data into one df 
    bass.mod.dat <- rbind(prXX_female, prXX_male, prXY_male) %>% 
      group_by(sex) %>% 
      mutate(sd = sd(Estimate),
             se= sd/sqrt(n()))
    # adding est error col with se
    bass.mod.dat$upper <- bass.mod.dat$Estimate + bass.mod.dat$se
    bass.mod.dat$lower <- bass.mod.dat$Estimate - bass.mod.dat$se    
    
############## ############### ############## 
#### df for summarizing raw datapoints #####
############## ############### ##############  
bass.raw.summary <- bassiana.data %>% 
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

    
####################
# Plots
####################
    # 1st ORDER: XYM
    mycolors <- c("#333333", "#990000", "#3399FF")
    ggplot(data = bass.raw.summary, aes(zlogMass, MR, group = sex, color= sex )) +
      geom_point(alpha =.6)+
      geom_errorbar(aes(ymin = MR-MR.se, ymax = MR+MR.se)) + 
      geom_errorbarh(aes(xmin = zlogMass-mass.se, xmax = zlogMass+mass.se))+
      geom_smooth(data = bass.mod.dat, aes(x=zlogMass, y=Estimate))+
      geom_ribbon(data = bass.mod.dat, aes(y = NULL, ymin = lower, ymax = upper, fill = sex), alpha = .5)+
      scale_fill_manual(values = mycolors, guide = FALSE) +
      scale_color_manual(values = mycolors, guide = FALSE) +
      theme_bw() +
      theme(axis.text = element_text(size=12)) +
      theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
      labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
    
    ### save plot ##
    ggsave(filename ="figures/bassiana..mod2.regression.pdf",  height = 5, width = 7)
      
      
    # Regression Plot accounting for log metabolic rate and log mass across sex
    ggplot(bassiana.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
      geom_smooth(method="lm") +
      geom_point(size=3, alpha = 0.2) +
      theme_bw() + 
      xlab("Log Mass") +
      ylab("Log Metabolic Rate") +
      ggtitle("Bassiana MR") 

    ### BRMS try
    bassiana.data %>%
      group_by(sex) %>%
      modelr::data_grid(zlogMass = zlogMass,ztime =ztime, day = day, id = id) %>%
      add_fitted_draws(Bas_m2_brms) %>%
      ggplot(aes(x = zlogMass, y = log(O2_min), color = ordered(sex))) +
      stat_lineribbon(aes(y = .value)) +
      geom_point(data = bassiana.data) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_brewer(palette = "Set2")
    dimnames(post_Bas_m2)
    
    
    
    
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

##############
## Model 1
##############
Pog_m1_brms <- brm(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                   family = gaussian(), data = pogona.data, iter= 2000, warmup = 1000, thin = 1, control = list(adapt_delta=0.95), cores = 4)
Pog_m1_brms <- add_criterion(Pog_m1_brms, c("loo", "waic"), moment_match = TRUE)
saveRDS(Pog_m1_brms, "./models/Pog_m1_brms")
summary(Pog_m1_brms)
bayes_R2(Pog_m1_brms)

##############
## Model 2
##############
        mod <- bf(log(O2_min) ~ sex*zlogMass + ztime + (1 + ztime | id) + (1 | day),  
                        sigma ~ zlogMass + ztime)
Pog_m2_brms <- brm(mod, family = gaussian(), data = pogona.data, iter= 2000, warmup = 1000, thin = 1, control = list(adapt_delta=0.95), cores = 4, save_all_pars = TRUE)
Pog_m2_brms <- add_criterion(Pog_m2_brms, c("loo", "waic"), moment_match = TRUE, reloo = TRUE)
saveRDS(Pog_m2_brms, "./models/Pog_m2_brms")

# import file
Pog_m2_brms <- readRDS(file="models/Pog_m2_brms")
summary(Pog_m2_brms)


# Model checks
plot(Pog_m2_brms)
summary(Pog_m2_brms)

#R2 of full model
bayes_R2(Pog_m2_brms)

######## ########### ######## 
### checking mod2 residuals ### 
######## ########### ######## 
resid.pred <- as.data.frame(predict(Pog_m2_brms))
ind.log.mr <-as.data.frame(log(pogona.data$O2_min)) %>% 
  rename(log.mr = `log(pogona.data$O2_min)`)

df <- bind_cols(ind.log.mr, resid.pred) %>% 
  mutate(diff = (log.mr - Estimate))
hist(df$diff)
ggplot(df, aes(x=diff)) + 
  geom_density()

# Prediction of model2 
Pog_m2_brms %>%
  spread_draws(b_Intercept, b_zlogMass) %>%
  mutate(zlogmass = list(seq(-0.4, 0.4, length.out = 100))) %>% #the observed value range of zlogmass
  unnest(zlogmass) %>%
  mutate(pred = exp(b_zlogMass + b_zlogMass*zlogmass)/(1+exp(b_zlogMass + b_zlogMass*zlogmass))) %>%
  group_by(zlogmass) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low = quantile(pred, prob = 0.025),
            pred_high = quantile(pred, prob = 0.975)) %>%
  ggplot(aes(x = zlogmass, y = pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha=0.2)
ggsave(filename ="figures/pogona.mod2.predicition.pdf",  height = 5, width = 7)


####################
# Model comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)

####################  
# extract posteriors for best model + Plotting
####################
post_Pog_m2_brms <- posterior_samples(Pog_m2_brms, "^b")
dimnames(post_Pog_m2_brms)

# Slope contrast
######### CHECK ZWF slope!!!!
ZZf <- as.array(post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZf:zlogMass"])
ZZm <- as.array(post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZm:zlogMass"])
ZWf <- as.array(post_Pog_m2_brms[,"b_zlogMass"])
pog.dat <- cbind(ZZf, ZZm, ZWf)

# plot
mcmc_intervals(pog.dat, 
               pars = c("ZWf", "ZZm", "ZZf"),
               prob = 0.95, 
               prob_outer = 0.99, 
               point_est = "mean")

mcmc_areas(
  pog.dat, 
  pars = c("ZWf", "ZZf", "ZZm"),
  prob = 0.95, 
  prob_outer = 0.99, 
  point_est = "mean")+ 
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope") 
ggsave(filename ="figures/pogona.mod2.posterior.pdf",  height = 5, width = 7)



# slope contrast phenotype
contrast_ZZ <- as.mcmc(ZZf - ZZm)
mean(contrast_ZZ)
HPDinterval(contrast_ZZ)
# slope contrast phenotype
contrast_ZWf_ZZf <- as.mcmc(ZWf - ZZf)
mean(contrast_ZWf_ZZf)
HPDinterval(contrast_ZWf_ZZf)



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
  mutate(sex = "ZWf")

#ZZ male 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.5029249, 1.23776, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZm <- prZZm %>% 
  mutate(sex = "ZZm")

#ZZ females
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.5029249, 1.23776, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZf <- prZZf %>% 
  mutate(sex = "ZZf")

# setting up prediction data into one df and calculating SE
mod.pog.dat <- rbind(prZWf, prZZm, prZZf) %>% 
  group_by(sex) %>% 
  mutate(sd = sd(Estimate),
       se= sd/sqrt(n()))

# adding est error col
mod.pog.dat$upper <- mod.pog.dat$Estimate + mod.pog.dat$se
mod.pog.dat$lower <- mod.pog.dat$Estimate - mod.pog.dat$se  


####################
# Plots
####################
mycolors <- c("#333333", "#990000", "#3399FF")
ggplot(data = pogona.data, aes(zlogMass, log(O2_min), group = sex, color= sex)) +
  geom_point(alpha = .6) +
  geom_line(data = mod.pog.dat, aes(x=zlogMass, y=Estimate))+
  geom_ribbon(data = mod.pog.dat, aes(y = NULL, ymin = lower, ymax = upper, fill = sex), alpha = .5)+
  scale_fill_manual(values = mycolors, guide = FALSE) +
  scale_color_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
  labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
ggsave(filename ="figures/pogona.mod2.regression.pdf",  height = 5, width = 7)
