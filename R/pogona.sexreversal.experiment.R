############################################################
# Code for processing Pogona SEX REVERSAL 20 June 2021
# extracts the change in CO2/O2
############################################################
# WD and packages
setwd("~/OneDrive - Australian National University/Respirometry")
library(plyr)
library(tidyverse)
library(metabR)
library(latex2exp)

########### MR Functions	 
# CO2 
MRCO2 <- function (v, time, bp, t, data) {
  VolCO2 = v * (data$change)
  k <- t + 273.15
  CorVolCO2 = (VolCO2 * bp * 273.15)/(101.325 * k)
  mCO2 <- (CorVolCO2)/time
  return(mCO2)
}
# O2 
MRO2 <- function (v, time, bp, t, data) {
  VolO2 = v * (data$change)
  k <- t + 273.15
  CorVolO2 = (VolO2 * bp * 273.15)/(101.325 * k)
  mO2 <- (CorVolO2)/time
  return(mO2)
}


## Find files
files <- paste0("./data/pogona/", list.files("./data/pogona/"))

# EXP data import 
data <- lapply(files, function(x) read.exp(x))
names(data) <- files
tau <- c(0.5)

# Plot the data; extract the peaks automatically. 
# Save the plots to pdf in figures output folder. 
# Need to check all auto-detection. 
plotting <- function(data, files, channel, col = "red", method = "spline", ...){
  quartz(width = 24.493392, height=  7.251101)
  gas <- plot_resp(data, channel = channel, col = col, main = files, method = method, ...)
  quartz.save(file = paste0(gsub("./data/pogona/", "./figures/pogona_resp_figures/", files), "_", channel, ".pdf"), type = "pdf")
  dat_gas <- resp_data(gas)
  return(dat_gas)
}

CO2_change <- mapply(function(x,y,z) plotting(x, y, channel = "CO2", tau = z, method = "spline"), x = data, y = files, z = tau, SIMPLIFY = FALSE)
O2_change <- mapply(function(x,y) plotting(x, y, channel = "O2", col = "blue"), x = data, y = files, SIMPLIFY = FALSE)

########## renaming to help with merging  data here in a few steps
file_names <- gsub("./resp data/", "", files)
file_names <- gsub(".exp", "", file_names)

## create new df with co2 and 02 data
datCO2 <- plyr::ldply(CO2_change)
datO2 <- plyr::ldply(O2_change)

### CO2 DATA 
# subset from data file name
datCO2$.id <- gsub("./data/pogona/", "", datCO2$.id)
datCO2$match <-  as.character(interaction(datCO2$.id, datCO2$marker))

### O2 DATA
# subset from data file name
datO2$.id <- gsub("./data/pogona/", "", datO2$.id)
datO2$match <-  as.character(interaction(datO2$.id, datO2$marker))

### METADATA 
# Bring in the meta-data file with weight etc
dat <- read.csv("final.pogona.metadata/pogona.metabolic.sexreversal.csv", stringsAsFactors = FALSE)
dat$marker <- as.character(tolower(dat$marker))
dat$match <- as.character(interaction(dat$file, dat$marker))

########### Merge MR data(CO2 & 02)and metadata together
final_CO2dat <- merge(dat, datCO2, by = "match")
final_O2dat <- merge(dat, datO2, by = "match")

######### cleaning CO2 & O2 data for MR function
# calculating mass for each df
final_CO2data<- final_CO2dat %>% 
  dplyr::group_by(date.dd.mm.yy., bd_liz_id) %>%
  dplyr::mutate(mass = mean(start_mass_g + end_mass_g))

# calculating mass for each df
final_O2data<- final_O2dat %>% 
  group_by(date.dd.mm.yy., bd_liz_id) %>%
  dplyr::mutate(mass = mean(start_mass_g + end_mass_g))


# calculating volume  by camber v 251.33 - mass of animal
final_CO2data$v <- final_CO2data$camber_vol_cm3 - final_CO2data$mass
final_O2data$v <- final_O2data$camber_vol_cm3 - final_O2data$mass

# final df's with MR for both CO2 and O2
# time needs to be converted from min to seconds
final_CO2data$MR_CO2_min <- MRCO2(v = final_CO2data$v, 
                                   time = (final_CO2data$time_in_chamber_min),
                                   bp = final_CO2data$bp, 
                                   t = final_CO2data$temp, 
                                   data = final_CO2data)

final_O2data$MR_O2_min <- MRO2(v = final_O2data$v, 
                                time = (final_O2data$time_in_chamber_min),
                                bp = final_O2data$bp, 
                                t = final_O2data$temp, 
                                data = final_O2data)

######### arranging final O2 data frame
O2_df <- final_O2data %>% 
  dplyr::group_by(date.dd.mm.yy., bd_liz_id) %>% 
  dplyr::arrange(marker_sample, .by_group = TRUE) %>% 
  dplyr::select(MR_O2_min,
                marker_sample, 
                FMS,
                Geno.pheno,
                Date.Hatched,
                Phenotypic.Sex,
                Genotypic.Sex,
                mass,
                svl_mm)

# removing 2 data points (2 hours) for each experiment
O2_1hr_removed <- O2_df %>%
  group_by(date.dd.mm.yy., bd_liz_id) %>%
  slice(-1)
O2_2hr_removed <- O2_1hr_removed %>%
  group_by(date.dd.mm.yy., bd_liz_id) %>%
  slice(-1)

# creating col with BMR (O2) now that first 2hrs of have been removed
O2_BMR <- O2_2hr_removed %>% 
  group_by(bd_liz_id, date.dd.mm.yy.) %>% 
  dplyr::mutate(Final_MR_O2_min = MR_O2_min/60, 
    BMR_O2 = min(Final_MR_O2_min)) %>% 
  filter(!is.na(Final_MR_O2_min))
  

# quick O2 Figure
fig <- ggplot(O2_BMR, aes(x = Geno.pheno, y = Final_MR_O2_min)) + 
  geom_point() +
  geom_violin() + 
  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
  geom_jitter(fill = "grey", alpha = 0.3) + 
  labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,O^2}{min}\\right)$"), x = "Treatment") + 
  theme_bw()
fig

write.csv(O2_BMR, file = "final.analysis.data/Pogona.finalO2.sexreversal.analysis.data.clean.csv")

######### arranging final CO2 data frame
CO2_df <- final_CO2data %>% 
  dplyr::group_by(date.dd.mm.yy., bd_liz_id) %>% 
  dplyr::arrange(marker_sample, .by_group = TRUE) %>% 
  dplyr::select(MR_CO2_min,
                marker_sample, 
                FMS,
                Geno.pheno,
                Date.Hatched,
                Phenotypic.Sex,
                Genotypic.Sex,
                mass,
                svl_mm)


# removing 2 data points (2 hours) for each experiment
# 1 hour
CO2_1hr_removed <- CO2_df %>%
  group_by(date.dd.mm.yy., bd_liz_id) %>%
  slice(-1)
# 2nd hour
CO2_2hr_removed <- CO2_1hr_removed %>%
  group_by(date.dd.mm.yy., bd_liz_id) %>%
  slice(-1)

# creating col with BMR (O2) now that first data point is removed
CO2_BMR <- CO2_2hr_removed %>% 
  group_by(date.dd.mm.yy., bd_liz_id) %>% 
  dplyr::mutate(Final_MR_CO2_min = MR_CO2_min/60, 
                BMR_CO2 = min(Final_MR_CO2_min))%>% 
  filter(!is.na(Final_MR_CO2_min))

fig <- ggplot(CO2_BMR, aes(x = Geno.pheno, y = Final_MR_CO2_min)) + 
  geom_point() +
  geom_violin() + 
  geom_jitter(fill = "grey", alpha = 0.3) + 
  labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,CO^2}{min}\\right)$"), x = "Treatment") + 
  theme_bw()
fig

write.csv(CO2_BMR, file = "final.analysis.data/Pogona.finalCO2.sexreversal.analysis.data.clean.csv")

