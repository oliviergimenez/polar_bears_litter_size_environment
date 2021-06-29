#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#                                 Figures                                      #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(lubridate)
library(raster)
library(grDevices)
Sys.setenv(LANG = "en")

# + 1. Temperature -------------------------------------------------------------

T_data <- read.csv("06_processed_data/temperature_data/temperature_processed_monthly") %>%
  mutate(date = as.Date(date))

# fun_color_range <- colorRampPalette(c("#56b1f7", "#132b43")) # colors from the default ggplot gradient (blue)
fun_color_range <- colorRampPalette(c("#85c7f9", "#0b1928"))
my_colors <- fun_color_range(29)

ggplot(data = T_data, aes(x = month, y = T_monthly, color = as.factor(year))) +
  geom_line(size = 0.3) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_color_manual(values = my_colors,
                     breaks = c("1990", "2000", "2010", "2018")) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("Janvier", "Avril", "Août", "Octobre")) +
  labs(x = "",
       y = "Température",
       color = "Année")

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/temperature_1990-2019.svg",
       width = 8, height = 4, unit = "cm")


# + 2. Sea ice concentration ----------------------------------------------------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData") 

years <- seq(from = 1990, to = 2020, by = 1)

daily_sea_ice_df <- data.frame(SI = unlist(daily_sea_ice_D_3),
                               date = seq(as.Date("1990-01-01"), as.Date("2020-12-31"), by="days")) %>%
  mutate(year = year(date),
         month = month(date))


monthly_sea_ice_df <- daily_sea_ice_df %>%
  group_by(month, year) %>%
  summarize(monthly_SI = mean(SI))

fun_color_range <- colorRampPalette(c("#85c7f9", "#0b1928"))
my_colors <- fun_color_range(31)

ggplot(data = monthly_sea_ice_df, aes(x = month, y = monthly_SI, color = as.factor(year))) +
  geom_line(size = 0.3) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_color_manual(values = my_colors,
                     breaks = c("1990", "2000", "2010", "2020")) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("Janvier", "Avril", "Août", "Octobre")) +
  labs(x = "",
       y = "Concentration de banquise",
       color = "Année")
# title = "Figure 3. Monthly mean temperature by year")

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/SI_concentration_1990-2019.svg",
       width = 8, height = 4, unit = "cm")



# + 3. Sea ice raster 1992 VS 2018 ---------------------------------------------

library(raster)
library(sp)
# Load the raster for year 1992 and 2018
load("04_raw_data/sea_ice/1992/raster_stack/sea_ice_raster_stack_1992.RData")
sea_ice_1992_50 <- sea_ice_raster_stack[[50]]

load("04_raw_data/sea_ice/2018/raster_stack/sea_ice_raster_stack_2018.RData")
sea_ice_2018_50 <- sea_ice_raster_stack[[50]]

#• Crop
crop_square <- extent(0, 2000000, -1500000, 500000)
sea_ice_1992_50_svalbard <- raster::crop(sea_ice_1992_50, crop_square)
sea_ice_2018_50_svalbard <- raster::crop(sea_ice_2018_50, crop_square)

writeRaster(sea_ice_2018_50_svalbard, "D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/2018_day_50.tif")


# Convert to df
df_sea_ice_1992 <- as.data.frame(sea_ice_1992_50_svalbard, 
                                 xy = TRUE) %>%
  mutate(year = 1992) %>%
  rename(value = sea_ice_1992_02_19)

df_sea_ice_2018 <- as.data.frame(sea_ice_2018_50_svalbard, 
                                 xy = TRUE) %>%
  mutate(year = 2018) %>%
  rename(value = sea_ice_2018_02_19)

f_sea_ice <- rbind(df_sea_ice_1992, df_sea_ice_2018)

# Plot
ggplot() +
  geom_raster(data = df_sea_ice, aes(x = x, y = y, fill = value)) +
  # geom_raster(data = df_sea_ice_1992, aes(x = x, y = y, fill = sea_ice_1992_02_19)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "",
       y = "",
       fill = "Couverture \nde banquise") +
  scale_fill_continuous(low = "#132b43", high = "white", 
                        guide = "colorbar", na.value = "#4d4d4d", # "#cccccc",
                        breaks = c(0, 50, 100),
                        labels = c("0%", "50%", "100%")) +
  facet_wrap(.~ year)


ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/raster_1992_VS_2018.svg",
       width = 16, height = 8, unit = "cm")


# + 4. Sea ice raster 1992 -----------------------------------------------------

library(raster)
library(sp)
# Load the raster for year 1992 and 2018
load("04_raw_data/sea_ice/1992/raster_stack/sea_ice_raster_stack_1992.RData")
sea_ice_1992_50 <- sea_ice_raster_stack[[50]]


# Crop
crop_square <- extent(0, 2000000, -1500000, 500000)
sea_ice_1992_50_svalbard <- raster::crop(sea_ice_1992_50, crop_square)


# Convert to df
df_sea_ice_1992 <- as.data.frame(sea_ice_1992_50_svalbard, 
                                 xy = TRUE) %>%
  rename(value = sea_ice_1992_02_19)


# Plot
ggplot() +
  geom_raster(data = df_sea_ice_1992, aes(x = x, y = y, fill = value)) +
  # geom_raster(data = df_sea_ice_1992, aes(x = x, y = y, fill = sea_ice_1992_02_19)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "",
       y = "",
       fill = "Couverture \nde banquise") +
  scale_fill_continuous(low = "#132b43", high = "white", 
                        guide = "colorbar", na.value = "#4d4d4d", # "#cccccc",
                        breaks = c(0, 50, 100),
                        labels = c("0%", "50%", "100%"))

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/raster_1992.svg",
       width = 8, height = 6, unit = "cm")

# + 5. Arctic oscillation ------------------------------------------------------

AO_data <- read.csv("06_processed_data/AO_data/AO_data_1990-2019_processed.csv") %>%
  mutate(date = as.Date(date))



ggplot(data = AO_data, aes(x = date, y = index, fill = as.factor(positive))) +
  geom_col() +
  theme_bw() +
  theme(legend.position = "",
        axis.title.x = element_text(size = 7, vjust = 5),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  scale_fill_manual(limits = c("1", "0"),
                    values = c("#132b43", "#56b1f7")) +
  labs(x = "Année",
       y = "Indice d'oscillation Arctique")

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/AO_1990-2020.svg",
       width = 8, height = 4, unit = "cm")



# + 6. Daily sea ice lineplot year 1992 ----------------------------------------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated.RData")
years <- seq(from = 1990, to = 2020, by = 1)

daily_sea_ice_D_3_df <- c()
for (k in 1:length(daily_sea_ice_D_3)) {
  daily_SI <- daily_sea_ice_D_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_D_3_df <- rbind(daily_sea_ice_D_3_df, df_k)
}
daily_sea_ice_D_3_df <- as.data.frame(daily_sea_ice_D_3_df) %>%
  filter(year == 1992)

ggplot(daily_sea_ice_D_3_df, aes(x = day_nbr, y = daily_SI)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/daily_sea_ice_lineplot.svg",
       width = 9, height = 6, unit = "cm")


# + 7. Daily sea ice lineplot all years ----------------------------------------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1990, to = 2020, by = 1)

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

daily_sea_ice_D_3_df <- c()
for (k in 1:length(daily_sea_ice_D_3)) {
  daily_SI <- daily_sea_ice_D_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_D_3_df <- rbind(daily_sea_ice_D_3_df, df_k)
}
daily_sea_ice_D_3_df <- as.data.frame(daily_sea_ice_D_3_df) %>%
  filter(year <2020)

ggplot(daily_sea_ice_D_3_df, aes(x = day_nbr, y = daily_SI, group = year)) +
  geom_line() +
  geom_hline(yintercept = get_halfway_sea_ice_concentration(daily_sea_ice_D_3),
             linetype="dashed", color = "red") +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/facetted_daily_sea_ice.svg",
       width = 9, height = 7.5)


# + 8. Date of sea ice retreat & advance ---------------------------------------

sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")

sea_ice_data %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     limits = c("day_retreat", "day_advance"),
                     labels = c("sea ice retreat", "sea ice advance")) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("D:/ETUDES/4A 2020-2021/Entretien PLR/graphes/sea_ice_retreat_advance.svg",
       width = 5, height = 3.2)

