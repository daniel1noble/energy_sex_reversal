############################################################
# Code DILUTION EXPERIMENT 20 MAY 2021
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
	files <- paste0("./data/dilutionexperiment/", list.files("./data/dilutionexperiment/"))

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
		quartz.save(file = paste0(gsub("./data/dilutionexperiment/", "./figures/dilutionexperiment_figures/", files), "_", channel, ".pdf"), type = "pdf")
		dat_gas <- resp_data(gas)
		return(dat_gas)
	}

	CO2_change <- mapply(function(x,y,z) plotting(x, y, channel = "CO2", tau = z, method = "spline"), x = data, y = files, z = tau, SIMPLIFY = FALSE)
	O2_change <- mapply(function(x,y) plotting(x, y, channel = "O2", col = "blue"), x = data, y = files, SIMPLIFY = FALSE)
	
## renaming to help with merging  data
	 file_names <- gsub("./resp data/", "", files)
	 file_names <- gsub(".exp", "", file_names)

# create new file with co2 and 02 data
datCO2 <- plyr::ldply(CO2_change)
datO2 <- plyr::ldply(O2_change)

########### CO2 
# subset from data file name
datCO2$.id <- gsub("./data/dilutionexperiment/", "", datCO2$.id)
datCO2$match <-  as.character(interaction(datCO2$.id, datCO2$marker))

########### O2 
# subset from data file name
datO2$.id <- gsub("./data/dilutionexperiment/", "", datO2$.id)
datO2$match <-  as.character(interaction(datO2$.id, datO2$marker))

########### METADATA 
# Bring in the meta-data file with weight etc.
meta_data <- read.csv("final.bassiana.metadata/bassiana.metabolic.dilution.experiment.csv", stringsAsFactors = FALSE)
meta_data$marker <- as.character(tolower(meta_data$marker))
meta_data$match <- as.character(interaction(meta_data$file, meta_data$marker))

########### Merge data together
final_CO2data <- merge(meta_data, datCO2, by = "match")
final_O2data <- merge(meta_data, datO2, by = "match")


######### cleaning CO2 & O2 data for MR function
# calculating volume  by camber v 23.56cm3 - mass of animal
final_CO2data$v <- with(final_CO2data, camber_vol_cm3 - end_mass_g)
final_O2data$v <- with(final_O2data, camber_vol_cm3 - end_mass_g)

# final df's with MR for both CO2 and O2
final_CO2data$MRCO2_final <- MRCO2(v = final_CO2data$v, 
                             time = final_CO2data$time_in_chamber_min, 
                             bp = final_CO2data$bp, 
                             t = final_CO2data$temp, 
                             data = final_CO2data)

final_O2data$MRO2_final <- MRO2(v = final_O2data$v, 
                          time = final_O2data$time_in_chamber_min, 
                          bp = final_O2data$bp, 
                          t = final_O2data$temp, 
                          data = final_O2data)


# summarizing by individual across experiments
df <- final_CO2data %>%
  dplyr::group_by(bd_liz_id, experiment,experiment_group, marker_sample, day) %>% 
  dplyr::summarise(MR_CO2_min = unique(MRCO2_final),
                   experiment = unique(experiment),
                   experiment_group = unique(experiment_group),
                   FMS = unique(FMS))
  

# removing first data point or hour for each experiment
# just check to be sure slice worked
final_df <- df %>%
  group_by(bd_liz_id, experiment, experiment_group) %>%
  slice(-2)

#check if the first row for each individual and experiment was remeoved
view(df)
view(final_df)

# had an na on one animal? from group_2 3x Anhydrous???
final_na_drop <- final_df %>% 
  drop_na()

# quick plot
fig <- ggplot(final_na_drop, aes(x = experiment, y = MR_CO2_min)) + 
geom_point() +
geom_violin() + 
geom_boxplot(width=0.1, color="red", alpha=0.2) + 
geom_jitter(fill = "grey", alpha = 0.3) + 
facet_grid(~experiment_group)+
labs(y = TeX("Metabolic Rate $\\left(\\frac{mL\\,CO^2}{min}\\right)$"), x = "Treatment") +
  scale_x_discrete(labels = c('1x Anhydrous','1x Depleted','3x Anhydrous', '3x Depleted', "Control")) + 
theme_bw()
fig

write.csv(final_na_drop, file = "final.dilution.analysis.data.clean.csv")





 