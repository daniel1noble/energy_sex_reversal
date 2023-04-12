# Reviewer correlations
pacman::p_load("dplyr", "tidyverse","brms")

# O2 - Pogona
pogona.data <- read.csv("~/Dropbox/energy_sex_reversal/final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv") %>% 
  rename(day =date.dd.mm.yy.,
         time = marker_sample,
         id = bd_liz_id, 
         O2_min = Final_MR_O2_min,
         sex = Geno.pheno,
         mass_g = mass) %>% 
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
pogona_final <- merge(x = pogona.data, y = pogona_growth, by = "id", all = TRUE)
bassiana_final <- merge(x = bassiana.data, y = bassiana_growth, by = "id", all = TRUE)



######################
# analysis for survival 
######################
######
# pogona
######
pog_surv_brms <- brm(status ~ log(mean_02), 
                     family = bernoulli(link = "logit"), 
                     data = pogona_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 4, cores = 8)
# assumption plot checks
plot(pog_surv_brms)
# summary & r2
summary(pog_surv_brms)
bayes_R2(pog_surv_brms)
# overall O2 on survival - high o2 high survival; dead = 0 alive = 1
plot(conditional_effects(pog_surv_brms, "mean_02"), ask = FALSE) 

######
# bassiana
######
bass_surv_brms <- brm(status ~ log(mean_02),  
                     family = bernoulli(link = "logit"), 
                     data = bassiana_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 4, cores = 8)
# assumption plot checks
plot(bass_surv_brms)
# summary & r2
summary(bass_surv_brms)
bayes_R2(bass_surv_brms)
# overall O2 on survival
plot(conditional_effects(bass_surv_brms, "mean_02"), ask = FALSE)



######################
# O2 on growth 
######################
######
# pogona
######
pog_o2_growth <- brm(log(mean_02) ~Growth.rate.mass., family = "gaussian", 
                     data = pogona_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 5, cores = 8) # warning for NA's are animals that died
# assumption plot checks
plot(pog_o2_growth)
# summary & r2
summary(pog_o2_growth)
bayes_R2(pog_o2_growth)
# overall O2 on growth - slight positive relationship with high MR and high GR
plot(conditional_effects(pog_o2_growth, "Growth.rate.mass."), ask = FALSE)

######
# Bassiana
######
bass_o2_growth <- brm(log(mean_02) ~ Growth.rate.mass., family = "gaussian", 
                     data = bassiana_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 5, cores = 8) # warning for NA's are animals that died
# assumption plot checks
plot(bass_o2_growth)
# summary & r2
summary(bass_o2_growth)
bayes_R2(bass_o2_growth)
# overall O2 on grwoth
plot(conditional_effects(bass_o2_growth, "Growth.rate.mass."), ask = FALSE)

