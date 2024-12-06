# Created by Selina Cheng
# Modeling how long oysters are open during "open" events

# Steps:
# Need to subset out data where there are NAs for predictors
# Stepwise backwards model selection for different models
# Compare models using AIC
# Report pseudo R squared value

# ----- Set up ------
rm(list = ls())
gc()

library(tidyverse)
library(lubridate)
library(GGally)
library(performance)
library(DHARMa)

source_dir <- "~/Library/CloudStorage/OneDrive-USNH/Oyster Biosensor/0_InProgressExperiments/2023 - 2024 Continuous system on oysters and mussels/2c_data for analysis"

oyster_dat_event <- fread(list.files(source_dir, pattern = "event", full.names=T))

# ---- # Amount of time open per time open is also Poisson negative binomial glm ---------
# oyster_dat_event

# Check for NAs in predictors
which(is.na(oyster_dat_event$mean_temp_c))
which(is.na(oyster_dat_event$mean_do))
which(is.na(oyster_dat_event$mean_sal))
which(is.na(oyster_dat_event$mean_ph))

# # Are predictors correlated with one another?
corrplot <- oyster_dat_event %>%
  dplyr::select(mean_temp_c, time_above_threshold:mean_ph)

# Create correlation plot
ggcorr(corrplot, palette = "RdBu", label = TRUE)
ggpairs(corrplot, lower=list(continuous='smooth'), axisLabels = 'show')

# Again, DO and temp are highly negatively correlated
# DO and salinity are highly negatively correlated
# Salinity and temp are positively correlated
# Salinity and ph are a little bit negatively correlated?

# Check VIF of pairs
r2 <- round(summary(lm(mean_sal ~ mean_do, data = oyster_dat_event))$r.squared, 1)
vif <- round(1 / (1 - r2), 1)
ggplot(oyster_dat_event, aes(x = mean_do, y = mean_sal)) + 
  geom_point() + 
  theme_classic() +
  ggtitle(paste("VIF = ", vif, ". R2 = ", 100*r2, "%", sep = ""))

# Drop DO? VIF of DO and temp is 5

# -------- Compare different temp variables to original temp data ------
oyster_dat_event <- oyster_dat_event %>%
  filter(grepl("open", event_id))

# Lag 1 day
# Subset data for no NAs in predictors
lag1_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag1),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag1_dat)

mod2 <- glm.nb(formula = length_minute ~ mean_lag1, data = lag1_dat)

AIC(mod1)
AIC(mod2)
# 1 day of lag is actually a better indicator here?  

# Lag 2 days
lag2_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag2),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag2_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag2, data = lag2_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag1, data = lag2_dat)

AIC(mod1)
AIC(mod2)
AIC(mod3)
# 1 day and 2 days of lag are similar AIC

# Lag 3 days
lag3_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag3),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag3_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag3, data = lag3_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag2, data = lag3_dat)
mod4 <- glm.nb(formula = length_minute ~ mean_lag1, data = lag3_dat)

AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)
# 3 day lag is best here 

# Lag 4 days
lag4_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag4),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag4_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag4_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag3, data = lag4_dat)
mod4 <- glm.nb(formula = length_minute ~ mean_lag2, data = lag4_dat)
mod5 <- glm.nb(formula = length_minute ~ mean_lag1, data = lag4_dat)

AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)
AIC(mod5)
# 4 day lag is best

# Lag 5 days
lag5_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag5),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag5_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag5, data = lag5_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag5_dat)
mod4 <- glm.nb(formula = length_minute ~ mean_lag3, data = lag5_dat)
mod5 <- glm.nb(formula = length_minute ~ mean_lag2, data = lag5_dat)
mod6 <- glm.nb(formula = length_minute ~ mean_lag1, data = lag5_dat)

AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)
AIC(mod5)
AIC(mod6)
# 4 and 5 day lag  are similar, 4 day lag still better (barely) (AIC lower)

# Lag 6 days
lag6_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag6),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag6_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag6, data = lag6_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag6_dat)
mod4 <- glm.nb(formula = length_minute ~ mean_lag5, data = lag6_dat)


AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)
# 4 day lag still performs better (AIC lower)

# Lag 7 days
lag7_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag7),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag7_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_lag7, data = lag7_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag7_dat)

AIC(mod1)
AIC(mod2)
AIC(mod3)
# 4 day lag still performs better (AIC lower)

# number of continuous minutes above threshold?
mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = oyster_dat_event)
mod2 <- glm.nb(formula = length_minute ~ time_above_threshold, data = oyster_dat_event)

AIC(mod1)
AIC(mod2)
# time_above_threshold does a good job (AIC lower)

# What if we compare temp vs date?
lag4_dat <- oyster_dat_event[!is.na(oyster_dat_event$mean_lag4),]

mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = lag4_dat)
mod2 <- glm.nb(formula = length_minute ~ start_time, data = lag4_dat)
mod3 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag4_dat)


AIC(mod1)
AIC(mod2)
AIC(mod3)
# start time is a better predictor...but adding time is a confusing/confounding variable

# What if we compare temp vs do?
mod1 <- glm.nb(formula = length_minute ~ mean_temp_c, data = oyster_dat_event)
mod2 <- glm.nb(formula = length_minute ~ mean_do, data = oyster_dat_event)

AIC(mod1)
AIC(mod2)
# And DO is a better predictor than temperature. huh

# Compare lag3, DO, start_time, time above threshold
mod1 <- glm.nb(formula = length_minute ~ mean_lag4, data = lag4_dat)
mod2 <- glm.nb(formula = length_minute ~ mean_do, data = lag4_dat)
mod3 <- glm.nb(formula = length_minute ~ time_above_threshold, data = lag4_dat)
mod4 <- glm.nb(formula = length_minute ~ start_time, data = lag4_dat)

AIC(mod1)
(mod1$null.deviance - mod1$deviance) / mod1$null.deviance
AIC(mod2)
(mod2$null.deviance - mod2$deviance) / mod2$null.deviance
AIC(mod3)
(mod3$null.deviance - mod3$deviance) / mod3$null.deviance
AIC(mod4)
(mod4$null.deviance - mod4$deviance) / mod4$null.deviance

# Hold on a minute...
check <- oyster_dat_event[oyster_dat_event$time_above_threshold == oyster_dat_event$length_minute,]
# Ok lol time above threshold is basically the same as oyster_dat_event, which is why it's such a "good predictor". Remove from model.

# Mean DO seems like a better predictor than temperature?

# ------------------- Model selection -----------------
mod1 <- glm.nb(formula = length_minute ~ mean_do + mean_sal + mean_ph, data = oyster_dat_event)

summary(mod1)
AIC(mod1)
# Calculate McFadden's pseudo R^2
1-(mod1$deviance/mod1$null.deviance)

# Salinity is not significant, drop?
mod2 <- glm.nb(formula = length_minute ~ mean_do + mean_ph, data = oyster_dat_event)
summary(mod2)
AIC(mod2)
# Calculate McFadden's pseudo R^2
1-(mod2$deviance/mod2$null.deviance)

# What if I drop ph?
mod3 <- glm.nb(formula = length_minute ~ mean_do, data = oyster_dat_event)
summary(mod3)
AIC(mod3)
# Calculate McFadden's pseudo R^2
1-(mod3$deviance/mod3$null.deviance)

# What if I substitute temp for do?
mod4a <- glm.nb(formula = length_minute ~ mean_lag4 + mean_ph, data = lag4_dat)
summary(mod4a)
AIC(mod4a)
# Calculate McFadden's pseudo R^2
1-(mod4a$deviance/mod4a$null.deviance)

mod4b <- glm.nb(formula = length_minute ~ mean_do + mean_ph, data = lag4_dat)
summary(mod4b)
AIC(mod4b)
# Calculate McFadden's pseudo R^2
1-(mod4b$deviance/mod4b$null.deviance)

# Model with DO is still lower AIC

# ------------------

# Mod 2 is best performing model based on AIC (lowest AIC)
summary(mod2)
AIC(mod2)
# Calculate McFadden's pseudo R^2
1-(mod2$deviance/mod2$null.deviance)

# Look at residuals and tests
simulationOutput <- simulateResiduals(mod2, n=1000)
plot(simulationOutput)  
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

# Prelim figure
ggplot(data = oyster_dat_event, aes(x = mean_do, y = length_minute))+
  geom_point()+
  stat_smooth(method = "glm.nb")






