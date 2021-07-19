################################################### 
# Analysis of sexreversal status bassiana experiment
################################################### 

# Packages
pacman::p_load("lme4", "tidyverse", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "emmeans", "latex2exp", "DHARMa", tidybayes)

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
    bayes_R2(post_Bas_m1)

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
    bayes_R2(Bas_m2_brms)

    # Model checks
    plot(Bas_m2_brms)
    summary(Bas_m2_brms)


####################
# Model comparison
####################
    loo_compare(Bas_m1_brms, Bas_m2_brms)

####################
# Predictions from best model for plots
####################
Bas_m2_brms <- readRDS(file = "models/Bas_m2_brms")
    #XY female
    newdata <- data.frame(
      sex = "XX_Female",
      zlogMass = seq(-0.4, 0.4, length.out = 100),
      ztime = 0)
    prXY_female <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXX_female <- prXX_female %>% 
      mutate(sex = "XX_female")
    
    #XY male 
    newdata <- data.frame(
      sex = "XY_Male",
      zlogMass = seq(-0.4, 0.4, length.out = 100),
      ztime = 0)
    prXY_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXY_male <- prXY_male %>% 
      mutate(sex = "XY_male")
    
    #XX male
    newdata <- data.frame(
      sex = "XX_Male",
      zlogMass = seq(-0.4, 0.4, length.out = 100),
      ztime = 0)
    prXX_male <- data.frame(cbind(predict(Bas_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
    prXX_male <- prXX_male %>% 
      mutate(sex = "XX_male")
    
    # setting up predicition data into one df 
    mod.dat <- rbind(prXX_female, prXX_male, prXY_male) %>% 
      group_by(sex)
    # adding est error col
    mod.dat$upper <- mod.dat$Estimate + mod.dat$Est.Error
    mod.dat$lower <- mod.dat$Estimate - mod.dat$Est.Error    
    
    
    
####################
# Plots
####################
    # ORDER: xxfemale, xxmale (SR), xymale
    mycolors <- c("#333333", "#990000", "#3399FF", "#333333", "#990000", "#3399FF")
    
    ### plot from predictions
     ggplot(data = bassiana.data, aes(zlogMass, log(O2_min), group = sex, color= sex)) +
        geom_point()  +
        geom_rug(sides = "b", size = 1) +
        geom_ribbon(data = mod.dat, aes(y = NULL, ymin = lower, ymax = upper, fill = sex), alpha = .5)+
        geom_smooth(data = mod.dat, aes(y = Estimate), size = 1) +
        scale_color_manual(values = mycolors, guide = FALSE) + 
        scale_fill_manual(values = mycolors, guide = FALSE) +
        theme_bw() +
        theme(axis.text = element_text(size=12)) +
        theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
        labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
      ### save plot ##
     ggsave(filename = "bassiana.pdf",  height = 5, width = 7)
      
      
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
bayes_R2(Pog_m2_brms)

## Posteriors
post_Pog_m2_brms <- posterior_samples(Pog_m2_brms, "^b")

# Slope contrast
b_sexZZf <- post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZf:zlogMass"]
b_sexZZm <- post_Pog_m2_brms[,"b_zlogMass"] + post_Pog_m2_brms[,"b_sexZZm:zlogMass"]

contrast_ZZ <- as.mcmc(b_sexZZf - b_sexZZm)
mean(contrast_ZZ)
HPDinterval(contrast_ZZ)

####################
# Model comparison
####################
loo_compare(Pog_m2_brms,Pog_m1_brms)


####################
# Predictions from best model for plots
####################
Pog_m2_brms <- readRDS(file="models/Pog_m2_brms")
#ZW female
newdata <- data.frame(
  sex = "ZWf",
  zlogMass = seq(-0.6, 1.3, length.out = 100),
  ztime = 0)
prZWf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZWf <- prZWf %>% 
  mutate(sex = "ZWf")

#ZZ male 
newdata <- data.frame(
  sex = "ZZm",
  zlogMass = seq(-0.6, 1.3, length.out = 100),
  ztime = 0)
prZZm <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZm <- prZZm %>% 
  mutate(sex = "ZZm")

#ZZ female
newdata <- data.frame(
  sex = "ZZf",
  zlogMass = seq(-0.6, 1.3, length.out = 100),
  ztime = 0)
prZZf <- data.frame(cbind(predict(Pog_m2_brms, newdata = newdata, re_formula = NA, summary = TRUE), zlogMass=newdata$zlogMass))
prZZf <- prZZf %>% 
  mutate(sex = "ZZf")

# setting up predicition data into one df 
mod.pog.dat <- rbind(prZWf, prZZm, prZZf) %>% 
  group_by(sex)
# adding est error col
mod.pog.dat$upper <- mod.pog.dat$Estimate + mod.pog.dat$Est.Error
mod.pog.dat$lower <- mod.pog.dat$Estimate - mod.pog.dat$Est.Error  


####################
# Plots
####################
# ORDER: ZWf, ZZm, ZZf(SR)
mycolors <- c("#333333", "#3399FF", "#990000", "#333333", "#3399FF", "#990000")

### plot from predictions
ggplot(data = pogona.data, aes(zlogMass, log(O2_min), group = sex, color= sex)) +
  geom_point()  +
  geom_rug(sides = "b", size = 1) +
  geom_ribbon(data = mod.pog.dat, aes(y = NULL, ymin = lower, ymax = upper, fill = sex), alpha = .5)+
  geom_smooth(data = mod.pog.dat, aes(y = Estimate), size = 1) +
  scale_color_manual(values = mycolors, guide = FALSE) + 
  scale_fill_manual(values = mycolors, guide = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))+
labs(y = TeX("log Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "log Mass (g)") 
### save plot ##
ggsave(filename = "pogona.pdf", height = 5, width = 7)


# Regression Plot accounting for log metabolic rate and log mass across sex
ggplot(pogona.data, aes(log(mass_g), log(O2_min), shape=sex, colour=sex, fill=sex)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Log Mass") +
  ylab("Log Metabolic Rate") +
  ggtitle("Pogona") 



