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
pogona_final <- merge(x = pogona.data, y = pogona_growth, by = "id", all = TRUE) %>% 
  rename(sex = sex.x)
bassiana_final <- merge(x = bassiana.data, y = bassiana_growth, by = "id", all = TRUE) %>% 
  rename(sex = sex.x)



######################
# analysis for survival 
######################
######
# pogona
######
pog_surv_brms <- brm(status ~ log(mean_02)+ Mass.1.g + sex, 
                     family = bernoulli(link = "logit"), 
                     data = pogona_final, 
                     iter= 5000, warmup = 1000, 
                     thin = 4, cores = 8)
# assumption plot checks
plot(pog_surv_brms)
# summary & r2
summary(pog_surv_brms)
bayes_R2(pog_surv_brms)

# extract 
pog_post <- posterior_samples(pog_surv_brms)
a<- pog_post[,2]
a <- as.array(a)
pmcmc(a)


# overall O2 on survival - high o2 high survival; dead = 0 alive = 1
plot(conditional_effects(pog_surv_brms, "mean_02"), ask = FALSE) 

######
# bassiana
######
bass_surv_brms <- brm(status ~ log(mean_02)+ Mass.1.g + sex,  
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




