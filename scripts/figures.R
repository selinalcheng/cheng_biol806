# Creating figures for BIOL 806 project (base figures were modified in Microsoft Powerpoint)

# ------ Set up --------
library(tidyverse)
library(lubridate)
library(data.table)

source_dir <- "~/Library/CloudStorage/OneDrive-USNH/Oyster Biosensor/0_InProgressExperiments/2023 - 2024 Continuous system on oysters and mussels/2c_data for analysis"

oyster_dat_1min <- fread(list.files(source_dir, pattern = "1min", full.names = T))

# ----- Theme ---------
selina_theme <- function(text_size){
  theme_classic()+
    theme(panel.grid.major.y = element_line(colour = "#ececec"),
          panel.grid.major.x = element_line(color = "#ececec"),
          text = element_text(size = (text_size + 1)),
          axis.text.x = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = text_size, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = text_size),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks = element_line(linewidth = 0.5, color = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.title.x = element_text(vjust= -1),
          legend.title = element_text(size = (text_size +1)),
          legend.text = element_text(size = text_size))
}

# -------- Create figures? ------------------
# Create a figure for the whole timeseries
fig1a <- ggplot(data = oyster_dat_1min, aes(x = minute_floor_est, y = mean_value_clean, color = new_status))+
  geom_point()+
  labs(x="Time", y = "Mean voltage (mV)", color = "Status")+
  scale_x_datetime(date_breaks = "1 month", 
                   date_labels="%b-%Y") +
  scale_color_manual(labels = c("Closed", "Open"), values = c("#f8766d", "#00bfc4")) +
  annotate("rect", xmin = as.POSIXct('06/01/2024', format = "%m/%d/%Y"), 
           xmax = as.POSIXct('06/07/2024', format = "%m/%d/%Y"),
           ymin = 490, ymax = 540,
           alpha = .1,fill = "blue")+
  selina_theme(14)

ggsave("figures/figure1_wholeseries.png", plot = fig1a, width = 12, height = 6, dpi = 300)

# Create a figure zoomed in on a week of data
fig1b <- ggplot(data = oyster_dat_1min, aes(x = minute_floor_est, y = mean_value_clean, color = new_status))+
  geom_point()+
  labs(x="Time", y = "Mean voltage (mV)", color = "Status")+
  scale_x_datetime(date_breaks = "1 day", 
                   limits = as.POSIXct(c('06/01/2024', '06/07/2024'), format="%m/%d/%Y"),
                   date_labels="%d-%b-%Y") +
  scale_color_manual(labels = c("Closed", "Open"), values = c("#f8766d", "#00bfc4")) +
  scale_y_continuous(limits = c(515, 535))+
  selina_theme(14)

ggsave("figures/figure1_zoom.png", plot = fig1b, width = 12, height = 6, dpi = 300)



# Create midpoints plot
# ggplot()+
#   geom_point(data = oyster6_dat, aes(x = minute_floor_est, y = mean_value_clean),
#              color= "blue")+
#   geom_point(data = oyster6_dat, aes(x = minute_floor_est, y= midpoints),
#              color = "red")+
#   geom_point(data = midpoint_dat, aes(x = minute_floor_est, y = midpoints), size = 3)