# Created by Selina Cheng
# Last modified 19 Nov 2024
# Detect erroneous status assignments and reassign status if erroneous
rm(list = ls())
gc()

# ---------- SET UP ------------
# load libraries
library(data.table)
library(tidyverse)
library(lubridate)
library(zoo)

# Set dirs
source_dir <- "~/Library/CloudStorage/OneDrive-USNH/Oyster Biosensor/0_InProgressExperiments/2023 - 2024 Continuous system on oysters and mussels/2a_processed sensor data"

output_dir <- "~/Library/CloudStorage/OneDrive-USNH/Oyster Biosensor/0_InProgressExperiments/2023 - 2024 Continuous system on oysters and mussels/2b_derived sensor data"

# ---------- READ DATA -----------
# Read in rough cleaned oyster data
oyster_dat <- fread(list.files(source_dir, pattern = "6", full.names=T))

# Remove NAs from data to get clean "runs" 
oyster_dat <- oyster_dat %>%
  dplyr::filter(!is.na(mean_value_clean))

# Get lengths and values of runs of equal values in a vector
# Get runs of status
runs <- rle(oyster_dat$status)
runs <- data.frame(status = runs$values, length = runs$lengths)

# Create diff (detect when it goes open to closed)
oyster_dat$change <- c(diff(oyster_dat$status_numeric),0)
# Are mean values above midpoints?
oyster_dat$above_midpoints <- oyster_dat$mean_value_clean > oyster_dat$midpoints

# Diffs represent status of (n+1) being different from (n)
# A diff of 1 means an oyster changed from closed to open
# A diff of -1 means an oyster changed from open to closed
# When change != 0, it coincides with the end of a status (not the beginning of a new status)

# Add where the change occurs so we know where to look
runs$change_rows <- c(which(oyster_dat$change != 0), NA)
runs$above_mid <- c(oyster_dat$above_midpoints[which(oyster_dat$change != 0)], T)

# Which rows correspond to May through September? -- this is where there is a lot of short runs of data due to spawning behavior? Use this to differentiate how following ifelse statement works
which(month(oyster_dat$minute_floor_est)==5)[1]

# Flag is TRUE if the runs are particularly short relative to what we think they should be
# If runs are after May 1 and are above midpoint line, only set flag to TRUE if runs are between 8 and 20 rows long
# Otherwise, runs less than 20 are TRUE
runs$flag <- ifelse(runs$change_rows-runs$length+1 >= 204439 & runs$above_mid, 
                    ifelse(runs$length < 20 & runs$length > 8, T, F), 
                    ifelse(runs$length < 20, T, F))

# NOW do runs of flag = T to identify when there are multiple switches back and forth of short "open" or "closed" stints. Oysters should not switch back and forth for such short amounts of time so much
runs_flag <- rle(runs$flag)
# Expand lengths
flag_lengths <- rep(runs_flag$lengths, times = runs_flag$lengths)
runs <- cbind(runs, flag_lengths)
# # When diff !=0 it is the start of a new group of T or F flags
runs$diff <- c(1, diff(runs$flag))
runs$diff <- ifelse(runs$diff != 0, "start", "group")

# Create new oyster status var (preserve old one)
oyster_dat$new_status <- oyster_dat$status

# Find where flag run lengths are greater than one (switching back and forth)
runs_fix <- runs %>% dplyr::filter(flag, flag_lengths > 1)
# Find starts of oscillation events
events <- which(runs_fix$diff == "start")

# If runs of true are greater than 3 (oscillating 3 times), change status
start <- Sys.time()
for(n in events){
  # Get first row of issues
    first_row <- runs_fix$change_rows[n]
  # get length of issues
    length <- runs_fix$flag_lengths[n]
  # Get last row of issues
    last_row <- runs_fix$change_rows[n+length-1]
    
    # Get length of first event
    length_dat <- runs_fix$length[n]
  # From first row to last row, change status to the one before the first row
    oyster_dat$new_status[(first_row-length_dat+1):last_row] <- oyster_dat$status[(first_row-length_dat)]
}
end <- Sys.time()
end-start

# The runs code should correct anywhere that there are a bunch of erroneous status assignments 
# Save runs data
oyster_dat %>%
  mutate(minute_floor_est = format(minute_floor_est, format = "%Y-%m-%d %H:%M:%S")) %>%
  write.table(file = file.path(output_dir, paste0("JEL-continuous-oyster6_fix.csv")),
              append = F, sep = ",", na = "NA", dec = ".", row.names = F, col.names = T)









