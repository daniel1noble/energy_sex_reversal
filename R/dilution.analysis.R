
################################################### 
# Analysis of anhydrous / dilution experiment
################################################### 
setwd("~/OneDrive - Australian National University/Respirometry")

# Packages
	pacman::p_load(lme4, MASS, tidyverse, brms, MCMCglmm, quantreg, lmerTest, emmeans, latex2exp, multcomp)

# Load data
	data <- read.csv("./final.dilution.analysis.data.clean.csv")

# Rename marker sample to time
	colnames(data)[5] <- "time"

# New variables
final_dat <- data %>% 
			mutate(vol = factor(substr(experiment, 1, 1)),
				   drierite = gsub("[3]_", "", experiment),
				   scl_time = scale(time),
				        FMS_2 = interaction(FMS, date),
				        experiment = factor(experiment, levels = c("Control", "1_anhydrous", "1_depleted", "3_anhydrous", "3_depleted")))

# quick plot
	fig <- ggplot(final_dat, aes(x = experiment, y = MR_CO2_min)) + 
	  geom_point() +
	  geom_violin() + 
	  geom_boxplot(width=0.1, color="red", alpha=0.2) + 
	  geom_jitter(fill = "grey", alpha = 0.3) + 
	  #facet_grid(~experiment_group) + 
	  theme_bw()
	fig

# Models 

model.lmer <- lmer(scale(MR_CO2_min) ~ experiment + scl_time + experiment:scl_time + (1 + scl_time | bd_liz_id) + (1 | date), data = final_dat)
summary(model.lmer)

AIC(model.lmer)
anova(model.lmer)

# OK, but within each volume / depleted pair, do we get signififcant differences. This is what we predcit if volume or anhydrous effect signals. 
emmeans(model.lmer, pairwise ~ experiment)

## We don;t see mean differences, because the slopes cancel each other out. in controla nd depleted they are negative, in the sense that animals are habituating and metabolism decreasing through time, but because anhydrous is grabbing CO2 we don't get the strong initial signal and as a consequence see the opposite effect on time. These would then cancel out potentially.

# check assumptions for final model
anova(model.lmer)
hist(residuals(model.lmer))
plot(residuals(model.lmer))

