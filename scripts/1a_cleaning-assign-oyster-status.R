# Created by Selina Cheng
# Last modified 6 Nov 2024
# Clean JEL data 
# Goals: Assign status of "open" or "closed" 
# I finagled a bunch of methods due to messy data...so it is not very clean
rm(list = ls())
gc()

# ---- Functions ----
# Create mean, SD, range functions 
mean_rm <- function(x){mean(x, na.rm = T)}
sd_rm <- function(x){sd(x, na.rm = T)}
range_rm <- function(x){
  range <- range(x, na.rm = T)
  return(diff(range))
  }

# Create mean function that trims upper and lower 25% 
mean_trim <- function(x){mean(x, trim = 0.25, na.rm = T)}

# -------- SET UP ----------
# load libraries
library(data.table)
library(tidyverse)
library(lubridate)
library(zoo)

# Set source directory
oyster_dir <- "~/Library/CloudStorage/OneDrive-USNH/Oyster Biosensor/0_InProgressExperiments/2023 - 2024 Continuous system on oysters and mussels/2a_processed sensor data"

# --------- LOAD DATA ----------
# Read in "rough" cleaned oyster data (code not included in project)
oyster_dat <- fread(list.files(oyster_dir, pattern = "rough", full.names=T))

# -------------- First round to nearest minute -------------------
# Round 12s data to 1min average oyster gape data
oyster_dat <- oyster_dat %>%
  # Create new column for minute flooring
  mutate(minute_floor_est = floor_date(timestamp_est, "minute"))

# OK..... I give up on the mussel data. Filter for just oysters
# I think mussel and oyster data have to be cleaned differently
oyster_dat <- oyster_dat %>%
  dplyr::filter(species == "oyster")

# Create a 1-min rounded dataset
start <- Sys.time()
oyster_dat_1min <- oyster_dat %>%
  group_by(minute_floor_est, unique_id, species, died) %>%
  summarise(mean_value_clean = mean_rm(value_clean),
            mean_value_raw = mean_rm(value_raw))
end <- Sys.time()
end-start

# Save data
# oyster_dat_1min %>%
#   mutate(minute_floor_est = format(minute_floor_est, format = "%Y-%m-%d %H:%M:%S")) %>%
#   write.table(file = file.path(oyster_dir, paste0("JEL-continuous-oyster_1min.csv")),
#               append = F, sep = ",", na = "NA", dec = ".", row.names = F, col.names = T)

# --------- Assign open/closed status -------------
# Read in 1min oyster data (from previous section)
oyster_dat_1min <- fread(list.files(oyster_dir, pattern = "1min", full.names=T))

# Choose just one oyster
oyster6_dat <- oyster_dat_1min %>%
  dplyr::filter(unique_id == "ID6-Y6")

# Use data.table library for faster processing
setDT(oyster6_dat)

# Functions to get the 5th and 95th quantile of the data
fifth_quantile <- function(x){
  y <- quantile(x, 0.05, na.rm = T)
  return(y)
}

ninetyfifth_quantile <- function(x){
  y <- quantile(x, 0.95, na.rm = T)
  return(y)
}

# Return midpoint ("midrange" of data)
middle <- function(x){
  midpt <- (diff(range(x, na.rm = T))*0.5)+min(x, na.rm = T)
  return(midpt)
  
  # if(length(midpt) > 1){
  #   warning("2 points returned")
  # } else{
  # }
}

# This function is misnamed...returns a point at 70% of max.
range10 <- function(x){
  midpt <- (diff(range(x, na.rm = T))*0.7)+min(x, na.rm = T)
  return(midpt)
  
  # if(length(midpt) > 1){
  #   warning("2 points returned")
  # } else{
  # }
}

# Create a line that roughly runs through the middle of the entire dataset (difficult due to changing baselines)
# Create random bins of data, find the rough midpoints of each bin and connect them.
# Create row ID
oyster6_dat$index <- c(1:nrow(oyster6_dat))

# Create separate bins for the winter (no changing baselines = fewer bins) and creates bins for the spring (changing baselines = need more bins to find accurate midpoints)
bins_winter <- seq(from = 1, to = oyster6_dat$index[which(date(oyster6_dat$minute_floor_est) == "2024-04-20")[1]], 
                   length.out = 6)
bins_spring <- seq(from = oyster6_dat$index[which(date(oyster6_dat$minute_floor_est) == "2024-05-05")[1]], 
                   to = nrow(oyster6_dat), length.out = 51)

bins_winter <- round(bins_winter)
bins_spring <- round(bins_spring)
bins <- c(bins_winter, bins_spring)

# Create empty vectors to fill with for loop
indices <- vector("numeric", length = length(bins)-1)
midpoints <- vector("numeric", length = length(bins)-1)

# For each bin, find midpoint of the interval. Save to empty vectors for midpoint and index
for(n in 2:length(bins)){
  midpoint <- middle(oyster6_dat$mean_value_clean[(bins[n-1]):bins[n]])
  midpoints[n-1] <- midpoint
  indices[n-1] <- round(((bins[n]-bins[n-1])/2)+bins[n-1])
}

# Add beginning and end midpoints for this "midpoint throughline" I'm making
indices <- c(1, indices)
indices <- c(indices, 393792)
midpoints <- c(530, midpoints)
midpoints <- c(midpoints, 525)

# Turn the empty vectors into a dataframe, add timestamps
midpoint_dat <- data.frame(index= indices, midpoints=midpoints)
midpoint_dat$minute_floor_est <- oyster6_dat$minute_floor_est[midpoint_dat$index]

# Remove a few points that are messing up the line
midpoint_dat <- midpoint_dat %>%
  dplyr::filter(month(minute_floor_est) != 1, month(minute_floor_est) != 4, month(minute_floor_est) != 3)

# Add some more points to the dataframe that I think would be helpful...
midpoint_dat <- midpoint_dat %>%
  add_row(index = oyster6_dat$index[which(date(oyster6_dat$minute_floor_est) == "2024-05-20")[1]],
          midpoints = 510, minute_floor_est = ymd_hms("2024-05-20 00:00:00")) %>%
  add_row(index = oyster6_dat$index[which(date(oyster6_dat$minute_floor_est) == "2024-04-20")[1]],
          midpoints = 526, minute_floor_est = ymd_hms("2024-04-20 00:00:00")) %>%
  add_row(index = oyster6_dat$index[which(date(oyster6_dat$minute_floor_est) == "2024-03-05")[1]],
          midpoints = 526, minute_floor_est = ymd_hms("2024-03-05 00:00:00"))
  
# Join midpoint dat with oyster data
oyster6_dat <- left_join(oyster6_dat, midpoint_dat[,-3], by = "index")

# Now linear interpolate between each midpoint to create a continuous line
oyster6_dat$midpoints <- na.approx(oyster6_dat$midpoints, na.rm = F)

# Find SD in each 1 hr moving window
start <- Sys.time()
oyster6_dat[, sd_day := frollapply(mean_value_clean, 1440, sd_rm, fill = NA, align = c("center"))]
end <- Sys.time()
end-start
# Fill in leading and trailing NAs for SD
oyster6_dat$sd_day[1:719] <- sd_rm(oyster6_dat$mean_value_clean[1:719])
oyster6_dat$sd_day[393073:393792] <- sd_rm(oyster6_dat$mean_value_clean[393073:393792])

# Create "70% of max" line for each rolling 36 hr window = "ten_pct_range" (misnomer...)
oyster6_dat[, ten_pct_range := frollapply(mean_value_clean, 2160, range10, fill = NA, align = c("center"))]
# Fill in leading and trailing NAs 
oyster6_dat$ten_pct_range[1:1079] <- range10(oyster6_dat$mean_value_clean[1:1079])
oyster6_dat$ten_pct_range[392713:393792] <- range10(oyster6_dat$mean_value_clean[392713:393792])

# Create narrow window 95th quant and 5th quant
# Moving window for 6 hours, find 5th quantile
oyster6_dat[, fifth_quant_narrow := frollapply(mean_value_clean, 360, fifth_quantile, fill = NA, align = c("center"))]
# Moving window for 6 hours, find 95th quantile
oyster6_dat[, ninetyfifth_quantile_narrow := frollapply(mean_value_clean, 360, ninetyfifth_quantile, fill = NA, align = c("center"))]
# Fill in leading and trailing NAs
oyster6_dat$fifth_quant_narrow[1:179] <- quantile(oyster6_dat$mean_value_clean[1:179], 0.05, na.rm = T)
oyster6_dat$fifth_quant_narrow[393613:393792] <- quantile(oyster6_dat$mean_value_clean[393613:393792], 0.05, na.rm = T)
oyster6_dat$ninetyfifth_quantile_narrow[1:179] <- quantile(oyster6_dat$mean_value_clean[1:179], 0.95, na.rm = T)
oyster6_dat$ninetyfifth_quantile_narrow[393613:393792] <- quantile(oyster6_dat$mean_value_clean[393613:393792], 0.95, na.rm = T)

# Assess variability over the day
# Variability is "low" if SD of the interval is less than 40% quantile of SD across the whole dataset
# OR variability is "low" if 95th quantile is very close to the 5th quantile (within 2)
oyster6_dat[, low_var := ifelse(sd_day < quantile(sd_day, 0.4, na.rm = T) | (ninetyfifth_quantile_narrow - fifth_quant_narrow) <= 2, "low", "high")]

# First, by default if mean value is below the "70% of max" line, assign open. otherwise, closed
oyster6_dat[, status := ifelse(mean_value_clean < ten_pct_range, "open", "closed")]
                               
# Then, if variability is low, look at if value is above or below midpoints line. If below, open. If above, closed
oyster6_dat[, status := ifelse(low_var == "low", ifelse(mean_value_clean > midpoints, "closed", "open"), status)]
oyster6_dat[, status_numeric := ifelse(status == "open", 1, 0)]

# Save data
oyster6_dat %>%
  mutate(minute_floor_est = format(minute_floor_est, format = "%Y-%m-%d %H:%M:%S")) %>%
  write.table(file = file.path(oyster_dir, paste0("JEL-continuous-oyster6.csv")),
              append = F, sep = ",", na = "NA", dec = ".", row.names = F, col.names = T)

