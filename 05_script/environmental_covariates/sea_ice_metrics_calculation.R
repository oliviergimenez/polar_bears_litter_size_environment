#==============================================================================#
#                                                                              #
#              Calculation of sea ice metrics used in the model                #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rgeos)
library(ff)
library(broom)
library(grDevices)
Sys.setenv(LANG = "en")


# A. Whole Barents Sea =========================================================

# ~ 1. Extract the daily sea ice concentration for the area for each year ------

barents_sea_boundary_polar <- readOGR("06_processed_data/subpopulation_boundary",
                                      layer = "barents_sea_subpop_polar")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice_A <- list()

# Less than 2-3min for 1991 -> 2020
start <- Sys.time()
for (i in 1:length(years)) {
  # Load the sea ice stack for year i
  sea_ice_stack_i <- raster::stack(x = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                                              "sea_ice_raster_stack_", years[i], ".tif"))
  load(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Load the ff matrix
              "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  daily_sea_ice_year_i <- rep(NA, times = nlayers(sea_ice_stack_i)) # Initialize the daily sea ice concentration vector for year i
  
  raster_ID <- raster(sea_ice_stack_i[[1]])        # create raster with same number of cells as the sea ice raster but empty
  raster_ID[] <- 1:ncell(sea_ice_stack_i[[1]])     # fill each cell with it's ID (run 'plot(raster_ID)' to better understand)
  cells_to_extract_ID <- extract(raster_ID,        # crop the raster containing the ID in its cells
                                 barents_sea_boundary_polar,
                                 weights = TRUE,          
                                 normalizeWeights = TRUE,   
                                 sp = FALSE)[[1]]
  all.days.extracted <- mat[cells_to_extract_ID[, 1], ] # Keep only the cells that correspond to these ID in the ff
                                                        # matrix
  NA.list <- which(is.na(all.days.extracted[, 1]))      
  all.days.extracted.clean <- all.days.extracted[-NA.list, ]  # Remove the ID of the cells that are NA (because with
                                                              # Weight = TRUE and normalizeWeights = TRUE, to get the
                                                              # mean SI concentration for each day, I can sum the
                                                              # "cell values*weight" if the sum of weights = 1. Need
                                                              # to remove the NA cells and re-adjust the weights accordingly
                                                              # so that their sum equals 1
  correction.cnst <- sum(cells_to_extract_ID[-NA.list, 2])
  for (j in 1:ncol(all.days.extracted)) {
    daily_sea_ice_year_i[j] <- sum(all.days.extracted.clean[, j] * (cells_to_extract_ID[-NA.list, 2]/correction.cnst))
    
  }
  
  daily_sea_ice_A[[i]] <- daily_sea_ice_year_i
  print(paste("done", years[i]))
  print(Sys.time())
}
end <- Sys.time()
end - start

save(daily_sea_ice_A, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_barents_sea.RData")


# The code above allows to go much faster than by cropping as usual. 




# ~ 2. Interpolate missing values ----------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_A_2 <- interpolate(daily_sea_ice_A, years)

save(daily_sea_ice_A_2,
     file = "06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated.RData")


# Plot

daily_sea_ice_A_df <- c()
for (k in 2:length(daily_sea_ice_A_2)) {
  daily_SI <- daily_sea_ice_A_2[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_A_df <- rbind(daily_sea_ice_A_df, df_k)
}
daily_sea_ice_A_df <- as.data.frame(daily_sea_ice_A_df)

# seq.months <- c(rep(1, times = 31), rep(2, times = 28),
#                 rep(3, times = 31), rep(4, times = 30),
#                 rep(5, times = 31), rep(6, times = 30),
#                 rep(7, times = 31), rep(8, times = 31),
#                 rep(9, times = 30), rep(10, times = 31),
#                 rep(11, times = 30), rep(12, times = 31))
# seq.months.leap <- c(rep(1, times = 31), rep(2, times = 29),
#                      rep(3, times = 31), rep(4, times = 30),
#                      rep(5, times = 31), rep(6, times = 30),
#                      rep(7, times = 31), rep(8, times = 31),
#                      rep(9, times = 30), rep(10, times = 31),
#                      rep(11, times = 30), rep(12, times = 31))
# 
# seq.days <- c(seq(1, 31, 1), seq(1, 28, 1), seq(1, 31, 1), seq(1, 30, 1),
#               seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1), seq(1, 31, 1),
#               seq(1, 30, 1), seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1))
# seq.days.leap <- c(seq(1, 31, 1), seq(1, 29, 1), seq(1, 31, 1), seq(1, 30, 1),
#                    seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1), seq(1, 31, 1),
#                    seq(1, 30, 1), seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1))
# 
# col.months <- c()
# col.days <- c()
# for (k in 1:length(unique(daily_sea_ice_A_df$year))) {
#   x <- nrow(daily_sea_ice_A_df[daily_sea_ice_A_df$year == unique(daily_sea_ice_A_df$year)[k], ])
#   if (x == 366) {
#     col.months <- c(col.months, seq.months.leap)
#     col.days <- c(col.days, seq.days.leap)
#     
#   } else {
#     col.months <- c(col.months, seq.months)
#     col.days <- c(col.days, seq.days)
#   }
# }
# 
# library(lubridate)
# daily_sea_ice_A_df <- daily_sea_ice_A_df %>%
#   mutate(month = col.months,
#          day = col.days,
#          date = as.Date(paste0(year, "-", month, "-", day, format = "%Y-%m-%d")))


ggplot(daily_sea_ice_A_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/daily_sea_ice_barents_sea.png",
       width = 10, height = 7.5)

ggplot(daily_sea_ice_A_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_barents_sea.png",
       width = 10, height = 7.5)

rm(daily_sea_ice_A_df)


# ~ 3. Correct wrong values ----------------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_A_3 <- correct_SI_mistakes(daily_sea_ice_A_2, years)

save(daily_sea_ice_A_3, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated_corrected.RData")


# Plot corrected daily SI 
daily_sea_ice_A_3_df <- c()
for (k in 2:length(daily_sea_ice_A_3)) {
  daily_SI <- daily_sea_ice_A_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_A_3_df <- rbind(daily_sea_ice_A_3_df, df_k)
}
daily_sea_ice_A_3_df <- as.data.frame(daily_sea_ice_A_3_df)

ggplot(daily_sea_ice_A_3_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_A_barents_sea.png",
       width = 10, height = 7.5)

# Check year 2009 and 2015
daily_sea_ice_A_3_df %>%
  filter(year %in% c(2009, 2015)) %>%
  ggplot() +
  geom_line(aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)



# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load data
load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run function
SI_retreat_advance_date_A <- get_SI_retreat_advance_date(daily_sea_ice_A_3, 
                                                         years)
# Save
write_csv(SI_retreat_advance_date_A, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_A.csv")

# plot
SI_retreat_advance_date_A %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_barents_sea.png",
       width = 6, height = 4)



# ~ 5. Plot approach A ---------------------------------------------------------

barents_sea_boundary_polar <- readOGR("06_processed_data/subpopulation_boundary",
                                      layer = "barents_sea_subpop_polar")

svalbard <- readOGR(dsn = "06_processed_data/svalbard_map",
                    layer = "svalbard_map_polar")


plot(barents_sea_boundary_polar)
plot(svalbard,
     add = TRUE)





# B. Buffer around Svalbard ====================================================

# ~ 1. Extract the daily sea ice concentration in buffer around islands --------

# ~~ a. Create the buffer around Svalbard --------------------------------------

svalbard.planar <- readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map_planar")
buffer.size <- 500000
svalbard.buffers.planar <- gBuffer(svalbard.planar, # generate the buffer around Svalbard
                                   byid = TRUE, 
                                   id = NULL, 
                                   width = buffer.size, 
                                   quadsegs = 5, 
                                   capStyle = "ROUND",
                                   joinStyle = "ROUND",
                                   mitreLimit = 1.0)
plot(svalbard.buffers.planar)
plot(svalbard.planar,
     add = TRUE)

sea_ice_stack <- raster::stack(x = "04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.tif")
svalbard.buffers <- spTransform(svalbard.buffers.planar, 
                                CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates



# ~~ b. Extract the daily mean sea ice concentration ---------------------------

years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice_B <- list()
start <- Sys.time()
for (i in 1:length(years)) {
  # Load the sea ice stack for year i
  sea_ice_stack_i <- raster::stack(x = paste0("04_raw_data/sea_ice/hamburg/", years[i], 
                                              "/raster_stack/", # Folder name
                                              "sea_ice_raster_stack_", years[i], ".tif"))
  load(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
              "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  daily_sea_ice_year_i <- rep(NA, times = nlayers(sea_ice_stack_i)) # Initialize the daily sea ice concentration vector for year i
  
  raster_ID <- raster(sea_ice_stack_i[[1]])
  raster_ID[] <- 1:ncell(sea_ice_stack_i[[1]])
  cells_to_extract_ID <- extract(raster_ID, 
                                 svalbard.buffers,
                                 weights = TRUE,          
                                 normalizeWeights = TRUE,   
                                 sp = FALSE)[[1]]
  
  
  all.days.extracted <- mat[cells_to_extract_ID[, 1], ]
  NA.list <- which(is.na(all.days.extracted[, 1]))
  all.days.extracted.clean <- all.days.extracted[-NA.list, ]
  
  correction.cnst <- sum(cells_to_extract_ID[-NA.list, 2])
  for (j in 1:ncol(all.days.extracted)) {
    daily_sea_ice_year_i[j] <- sum(all.days.extracted.clean[, j] * 
                                     (cells_to_extract_ID[-NA.list, 2]/correction.cnst))
    
  }
  
  daily_sea_ice_B[[i]] <- daily_sea_ice_year_i
  print(paste("done", years[i]))
  print(Sys.time())
}
end <- Sys.time()
end - start

save(daily_sea_ice_B, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer.RData")



# ~~ c. Plot the daily mean sea ice concentration -----------------------------+

load("06_processed_data/sea_ice_data/daily_sea_ice_100km_buffer.RData")

daily_sea_ice_B_df <- c()
for (k in 2:length(daily_sea_ice_B)) {
  daily_SI <- daily_sea_ice_B[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_B_df <- rbind(daily_sea_ice_B_df, df_k)
}
daily_sea_ice_B_df <- as.data.frame(daily_sea_ice_B_df)


seq.months <- c(rep(1, times = 31), rep(2, times = 28),
                rep(3, times = 31), rep(4, times = 30),
                rep(5, times = 31), rep(6, times = 30),
                rep(7, times = 31), rep(8, times = 31),
                rep(9, times = 30), rep(10, times = 31),
                rep(11, times = 30), rep(12, times = 31))
seq.months.leap <- c(rep(1, times = 31), rep(2, times = 29),
                     rep(3, times = 31), rep(4, times = 30),
                     rep(5, times = 31), rep(6, times = 30),
                     rep(7, times = 31), rep(8, times = 31),
                     rep(9, times = 30), rep(10, times = 31),
                     rep(11, times = 30), rep(12, times = 31))

seq.days <- c(seq(1, 31, 1), seq(1, 28, 1), seq(1, 31, 1), seq(1, 30, 1),
              seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1), seq(1, 31, 1),
              seq(1, 30, 1), seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1))
seq.days.leap <- c(seq(1, 31, 1), seq(1, 29, 1), seq(1, 31, 1), seq(1, 30, 1),
                   seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1), seq(1, 31, 1),
                   seq(1, 30, 1), seq(1, 31, 1), seq(1, 30, 1), seq(1, 31, 1))

col.months <- c()
col.days <- c()
for (k in 1:length(unique(daily_sea_ice_B_df$year))) {
  x <- nrow(daily_sea_ice_B_df[daily_sea_ice_B_df$year == unique(daily_sea_ice_B_df$year)[k], ])
  if (x == 366) {
    col.months <- c(col.months, seq.months.leap)
    col.days <- c(col.days, seq.days.leap)
    
  } else {
    col.months <- c(col.months, seq.months)
    col.days <- c(col.days, seq.days)
  }
}

library(lubridate)
daily_sea_ice_B_df <- daily_sea_ice_B_df %>%
  mutate(month = col.months,
         day = col.days,
         date = as.Date(paste0(year, "-", month, "-", day, format = "%Y-%m-%d")))


ggplot(daily_sea_ice_B_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration")



# ~ 1. Interpolate missing values ----------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_B_2 <- interpolate(daily_sea_ice_B, years)

save(daily_sea_ice_B_2,
     file = "06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer_interpolated.RData")


# Plot
daily_sea_ice_B_df <- c()
for (k in 2:length(daily_sea_ice_B_2)) {
  daily_SI <- daily_sea_ice_B_2[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_B_df <- rbind(daily_sea_ice_B_df, df_k)
}
daily_sea_ice_B_df <- as.data.frame(daily_sea_ice_B_df)


ggplot(daily_sea_ice_B_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/daily_sea_ice_500km_buffer.png",
       width = 10, height = 7.5)

ggplot(daily_sea_ice_B_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_500km_buffer.png.png",
       width = 10, height = 7.5)

rm(daily_sea_ice_B_df)



# ~ 3. Correct wrong values ----------------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_B_3 <- correct_SI_mistakes(daily_sea_ice_B_2, years)

save(daily_sea_ice_B_3, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer_interpolated_corrected.RData")


# Plot corrected daily SI 
daily_sea_ice_B_3_df <- c()
for (k in 2:length(daily_sea_ice_B_3)) {
  daily_SI <- daily_sea_ice_B_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_B_3_df <- rbind(daily_sea_ice_B_3_df, df_k)
}
daily_sea_ice_B_3_df <- as.data.frame(daily_sea_ice_B_3_df)

ggplot(daily_sea_ice_B_3_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_B_500km_buffer.png",
       width = 10, height = 7.5)



# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
SI_retreat_advance_date_B <- get_SI_retreat_advance_date(daily_sea_ice_B_3, 
                                                         years)

write_csv(SI_retreat_advance_date_B, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_B.csv")

# Plot
SI_retreat_advance_date_B %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_500km_buffer.png",
       width = 6, height = 4)



# ~ 5. Plot approach B ---------------------------------------------------------


svalbard <- readOGR(dsn = "06_processed_data/svalbard_map",
                    layer = "svalbard_map_polar")


plot(svalbard.buffers)
plot(svalbard,
     add = TRUE)










# C. Barents sea where depth <300m =============================================

# ~ 1. Extract the daily sea ice concentration for the area for each year ------

# ~~ a. Get ID of cells with <300m depth ---------------------------------------
barents_sea_boundary_polar <- readOGR("06_processed_data/subpopulation_boundary",
                                      layer = "barents_sea_subpop_polar")

bathymetry_barents_sea <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped_Barents_sea_final.tif")

# As the sea ice raster stacks are not cropped following the Barents Sea boundary,
# I need to make artificially extend the bathymetry raster to the same size, so that
# the cell IDs that I will retrieve will be the right ones for the sea ice raster

crop.square <- extent(-1000000, 2500000, -2000000, 1000000) # Same extent as that of the sea ice rasterstack
bath_extended <- extend(x = bathymetry_barents_sea, 
                        y = crop.square, 
                        value = -5000)
# Extract the cell with depth >= -300 or with NA
shallow_cells <- Which(bath_extended >= -300, cells = TRUE)
NA_cells <- Which(is.na(bath_extended), cells = TRUE)



# ~~ b. Extract the daily mean sea ice concentration ---------------------------

years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice_C <- list()

start <- Sys.time()
for (i in 1:length(years)) {
  # Load the sea ice stack for year i
  sea_ice_stack_i <- raster::stack(x = paste0("04_raw_data/sea_ice/hamburg/", years[i], 
                                              "/raster_stack/", # Folder name
                                              "sea_ice_raster_stack_", years[i], ".tif"))
  load(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
              "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  daily_sea_ice_year_i <- rep(NA, times = nlayers(sea_ice_stack_i)) # Initialize the daily sea ice concentration vector for year i
  
  raster_ID <- raster(sea_ice_stack_i[[1]])
  raster_ID[] <- 1:ncell(sea_ice_stack_i[[1]])
  cells_to_extract_ID <- as.data.frame(extract(raster_ID, 
                                       barents_sea_boundary_polar,
                                       weights = TRUE,          
                                       normalizeWeights = TRUE,   
                                       sp = FALSE)[[1]]) %>%
    filter(value %in% c(shallow_cells, NA_cells))
  
  all.days.extracted <- mat[cells_to_extract_ID[, 1], ]
  NA.list <- which(is.na(all.days.extracted[, 1]))
  all.days.extracted.clean <- all.days.extracted[-NA.list, ]
  
  correction.cnst <- sum(cells_to_extract_ID[-NA.list, 2])
  for (j in 1:ncol(all.days.extracted)) {
    daily_sea_ice_year_i[j] <- sum(all.days.extracted.clean[, j] * 
                                     (cells_to_extract_ID[-NA.list, 2]/correction.cnst))
  }
  
  daily_sea_ice_C[[i]] <- daily_sea_ice_year_i
  print(paste("done", years[i]))
  print(Sys.time())
}
end <- Sys.time()
end - start

save(daily_sea_ice_C, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea.RData")

rm(barents_sea_boundary_polar, bathymetry_barents_sea, crop.square, bath_extended, 
   shallow_cells, NA_cells, correction.cnst, daily_sea_ice_year_i, raster_ID, mat, 
   sea_ice_stack_i, cells_to_extract_ID, NA.list, all.days.extracted, 
   all.days.extracted.clean, i, j)



daily_sea_ice_C_df <- c()
for (k in 2:length(daily_sea_ice_C)) {
  daily_SI <- daily_sea_ice_C[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_C_df <- rbind(daily_sea_ice_C_df, df_k)
}
daily_sea_ice_C_df <- as.data.frame(daily_sea_ice_C_df)

ggplot(daily_sea_ice_C_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

rm(df_k, daily_sea_ice_C_df, daily_SI, day_nbr, year, k)


# ~ 2. Interpolate missing values ----------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_C_2 <- interpolate(daily_sea_ice_C, years)

save(daily_sea_ice_C_2,
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea_interpolated.RData")



# ~ 3. Correct wrong values ----------------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_C_3 <- correct_SI_mistakes(daily_sea_ice_C_2, years)

save(daily_sea_ice_C_3, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea_interpolated_corrected.RData")


# Plot corrected daily SI 
daily_sea_ice_C_3_df <- c()
for (k in 2:length(daily_sea_ice_C_3)) {
  daily_SI <- daily_sea_ice_C_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_C_3_df <- rbind(daily_sea_ice_C_3_df, df_k)
}
daily_sea_ice_C_3_df <- as.data.frame(daily_sea_ice_C_3_df)

ggplot(daily_sea_ice_C_3_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_C_300m_depth_barents_sea.png",
       width = 10, height = 7.5)




# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run function
SI_retreat_advance_date_C <- get_SI_retreat_advance_date(daily_sea_ice_C_3, 
                                                         years)
# Save
write_csv(SI_retreat_advance_date_C, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_C.csv")

# plot
SI_retreat_advance_date_C %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_depth_300m_barents_sea.png",
       width = 6, height = 4)


# ~ 5. Plot approach C ---------------------------------------------------------

bathymetry_barents_sea <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped_Barents_sea_final.tif")

barents_sea_boundary_polar <- readOGR("06_processed_data/subpopulation_boundary",
                                      layer = "barents_sea_subpop_polar")

svalbard <- readOGR(dsn = "06_processed_data/svalbard_map",
                    layer = "svalbard_map_polar")

cuts <- c(0, -300, -1000, -2000, -3000, -6000) #set breaks
pal <- colorRampPalette(c("black", "grey99"))
plot(bathymetry_barents_sea,
     breaks = cuts, 
     col = pal(7))
plot(barents_sea_boundary_polar,
     add = TRUE)
plot(svalbard,
     add = TRUE)




# D. Buffer around Svalbard where depth <300m ==================================

# ~ 1. Extract the daily sea ice concentration for the area for each year ------

# ~~ a. Create the buffer around Svalbard --------------------------------------

svalbard.planar <- readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map_planar")
buffer.size <- 500000
svalbard.buffers.planar <- gBuffer(svalbard.planar, # generate the buffer around Svalbard
                                   byid = TRUE, 
                                   id = NULL, 
                                   width = buffer.size, 
                                   quadsegs = 5, 
                                   capStyle = "ROUND",
                                   joinStyle = "ROUND",
                                   mitreLimit = 1.0)

plot(svalbard.buffers.planar)
plot(svalbard.planar,
     add = TRUE)

sea_ice_stack <- raster::stack(x = "04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.tif")
svalbard.buffers <- spTransform(svalbard.buffers.planar, 
                                CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates
svalbard <- spTransform(svalbard.planar, 
                        CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates

# ~~ b. Get ID of cells with <300m depth ---------------------------------------

bathymetry_barents_sea <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped_Barents_sea_final.tif")

# As the sea ice raster stacks are not cropped following the Barents Sea boundary,
# I need to make artificially extend the bathymetry raster to the same size, so that
# the cell IDs that I will retrieve will be the right ones for the sea ice raster

crop.square <- extent(-1000000, 2500000, -2000000, 1000000) # Same extent as that of the sea ice rasterstack
bath_extended <- extend(x = bathymetry_barents_sea, 
                        y = crop.square, 
                        value = -5000)
# Extract the cell with depth >= -300 or with NA
shallow_cells <- Which(bath_extended >= -300, cells = TRUE)
NA_cells <- Which(is.na(bath_extended), cells = TRUE)


# ~~ c. Extract the daily mean sea ice concentration ---------------------------

years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice_D <- list()

start <- Sys.time()
for (i in 1:length(years)) {
  # Load the sea ice stack for year i
  sea_ice_stack_i <- raster::stack(x = paste0("04_raw_data/sea_ice/hamburg/", years[i], 
                                              "/raster_stack/", # Folder name
                                              "sea_ice_raster_stack_", years[i], ".tif"))
  load(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
              "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  daily_sea_ice_year_i <- rep(NA, times = nlayers(sea_ice_stack_i)) # Initialize the daily sea ice concentration vector for year i
  
  raster_ID <- raster(sea_ice_stack_i[[1]])
  raster_ID[] <- 1:ncell(sea_ice_stack_i[[1]])
  cells_to_extract_ID <- as.data.frame(extract(raster_ID, 
                                               svalbard.buffers,
                                               weights = TRUE,          
                                               normalizeWeights = TRUE,   
                                               sp = FALSE)[[1]]) %>%
    filter(value %in% c(shallow_cells, NA_cells))
  
  all.days.extracted <- mat[cells_to_extract_ID[, 1], ]
  NA.list <- which(is.na(all.days.extracted[, 1]))
  all.days.extracted.clean <- all.days.extracted[-NA.list, ]
  
  correction.cnst <- sum(cells_to_extract_ID[-NA.list, 2])
  for (j in 1:ncol(all.days.extracted)) {
    daily_sea_ice_year_i[j] <- sum(all.days.extracted.clean[, j] * 
                                     (cells_to_extract_ID[-NA.list, 2]/correction.cnst))
  }
  
  daily_sea_ice_D[[i]] <- daily_sea_ice_year_i
  print(paste("done", years[i]))
  print(Sys.time())
}
end <- Sys.time()
end - start

save(daily_sea_ice_D, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer.RData")

rm(bathymetry_barents_sea, crop.square, bath_extended, buffer.size,
   svalbard.buffers.planar, svalbard.buffers, svalbard.planar,
   shallow_cells, NA_cells, correction.cnst, daily_sea_ice_year_i, raster_ID, mat, 
   sea_ice_stack_i, cells_to_extract_ID, NA.list, all.days.extracted, 
   all.days.extracted.clean, i, j)



daily_sea_ice_D_df <- c()
for (k in 2:length(daily_sea_ice_D)) {
  daily_SI <- daily_sea_ice_D[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_D_df <- rbind(daily_sea_ice_D_df, df_k)
}
daily_sea_ice_D_df <- as.data.frame(daily_sea_ice_D_df)

ggplot(daily_sea_ice_D_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

rm(df_k, daily_sea_ice_D_df, daily_SI, day_nbr, year, k)


# ~ 2. Interpolate missing values ----------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_D_2 <- interpolate(daily_sea_ice_D, years)

save(daily_sea_ice_D_2,
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer_interpolated.RData")



# ~ 3. Correct wrong values ----------------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_D_3 <- correct_SI_mistakes(daily_sea_ice_D_2, years)

save(daily_sea_ice_D_3, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer_interpolated_corrected.RData")


# Plot corrected daily SI 
daily_sea_ice_D_3_df <- c()
for (k in 2:length(daily_sea_ice_D_3)) {
  daily_SI <- daily_sea_ice_D_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_D_3_df <- rbind(daily_sea_ice_D_3_df, df_k)
}
daily_sea_ice_D_3_df <- as.data.frame(daily_sea_ice_D_3_df)

ggplot(daily_sea_ice_D_3_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_D_300m_depth_500km_buffer.png",
       width = 10, height = 7.5)




# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load data
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_300m_500km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run function
SI_retreat_advance_date_D <- get_SI_retreat_advance_date(daily_sea_ice_D_3, 
                                                         years)
# Save
write_csv(SI_retreat_advance_date_D, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")

# plot
SI_retreat_advance_date_D %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_depth_300m_500km_buffer.png",
       width = 6, height = 4)


# ~ 5. Plot approach D ---------------------------------------------------------

bathymetry_barents_sea <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped_Barents_sea_final.tif")


svalbard.planar <- readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map_planar")
buffer.size <- 500000
svalbard.buffers.planar <- gBuffer(svalbard.planar, # generate the buffer around Svalbard
                                   byid = TRUE, 
                                   id = NULL, 
                                   width = buffer.size, 
                                   quadsegs = 5, 
                                   capStyle = "ROUND",
                                   joinStyle = "ROUND",
                                   mitreLimit = 1.0)

sea_ice_stack <- raster::stack(x = "04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.tif")
svalbard.buffers <- spTransform(svalbard.buffers.planar, 
                                CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates
svalbard <- spTransform(svalbard.planar, 
                        CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates



cuts <- c(0, -300, -1000, -2000, -3000, -6000) #set breaks
pal <- colorRampPalette(c("black", "grey99"))
plot(bathymetry_barents_sea,
     breaks = cuts, 
     col = pal(7))
plot(svalbard.buffers,
     add = TRUE)
plot(svalbard,
     add = TRUE)


bathymetry_df <- as.data.frame(bathymetry_barents_sea, 
                               xy = TRUE) %>%
  rename(depth = ICBAO_bathymetry_polar_cropped_Barents_sea_final)
svalbard_df <- tidy(svalbard) 

buffers_df <- tidy(svalbard.buffers) 

ggplot() +
  geom_raster(data = bathymetry_df, aes(x = x, y = y, fill = depth)) +
  
  geom_polygon(data = svalbard_df, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 1) +
  geom_polygon(data = buffers_df, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 1) +
  scale_fill_continuous(low = "#132b43", high = "lightblue",
                        guide = "colorbar", na.value = "white",
                        breaks = c(0, 50, 100),
                        labels = c("0%", "50%", "100%")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  labs(x = "",
       y = "",
       fill = "Sea ice \ncover")


# E. Buffer around Svalbard + pelagic bears ====================================

# ~ 1. Extract the daily sea ice concentration for the area for each year ------

# ~~ a. Create the buffer around Svalbard --------------------------------------

svalbard.planar <- readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map_planar")
buffer.size <- 500000
svalbard.buffers.planar <- gBuffer(svalbard.planar, # generate the buffer around Svalbard
                                   byid = TRUE, 
                                   id = NULL, 
                                   width = buffer.size, 
                                   quadsegs = 5, 
                                   capStyle = "ROUND",
                                   joinStyle = "ROUND",
                                   mitreLimit = 1.0)

# Reproject the buffer into polar coordinates
sea_ice_stack <- raster::stack(x = "04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.tif")
svalbard.buffers <- spTransform(svalbard.buffers.planar, 
                                CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates
svalbard <- spTransform(svalbard.planar, 
                        CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates


# ~~ b. Add the pelagic bear area ----------------------------------------------

pelagic_bears_area_polar <- readOGR(dsn = "06_processed_data/sea_ice_data",
                              layer = "pelagic_bears_area_polar")

study_area_E <- aggregate(rbind(svalbard.buffers, pelagic_bears_area_polar))


# ~~ c. Extract the daily mean sea ice concentration ---------------------------

years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice_E <- list()
start <- Sys.time()
for (i in 1:length(years)) {
  # Load the sea ice stack for year i
  sea_ice_stack_i <- raster::stack(x = paste0("04_raw_data/sea_ice/hamburg/", years[i], 
                                              "/raster_stack/", # Folder name
                                              "sea_ice_raster_stack_", years[i], ".tif"))
  load(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
              "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  daily_sea_ice_year_i <- rep(NA, times = nlayers(sea_ice_stack_i)) # Initialize the daily sea ice concentration vector for year i
  
  raster_ID <- raster(sea_ice_stack_i[[1]])
  raster_ID[] <- 1:ncell(sea_ice_stack_i[[1]])
  cells_to_extract_ID <- raster::extract(raster_ID, 
                                         study_area_E,
                                         weights = TRUE,          
                                         normalizeWeights = TRUE,   
                                         sp = FALSE)[[1]]
  
  
  all.days.extracted <- mat[cells_to_extract_ID[, 1], ]
  NA.list <- which(is.na(all.days.extracted[, 1]))
  all.days.extracted.clean <- all.days.extracted[-NA.list, ]
  
  correction.cnst <- sum(cells_to_extract_ID[-NA.list, 2])
  for (j in 1:ncol(all.days.extracted)) {
    daily_sea_ice_year_i[j] <- sum(all.days.extracted.clean[, j] * 
                                     (cells_to_extract_ID[-NA.list, 2]/correction.cnst))
    
  }
  
  daily_sea_ice_E[[i]] <- daily_sea_ice_year_i
  print(paste("done", years[i]))
}
end <- Sys.time()
end - start

save(daily_sea_ice_E, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area.RData")


# ~ 2. Interpolate missing values ----------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_E_2 <- interpolate(daily_sea_ice_E, years)

save(daily_sea_ice_E_2,
     file = "06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated.RData")

# ~ 3. Correct wrong values ----------------------------------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load the data
load("06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run the function
daily_sea_ice_E_3 <- correct_SI_mistakes(daily_sea_ice_E_2, years)


save(daily_sea_ice_E_3, 
     file = "06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated_corrected.RData")


# Plot corrected daily SI 
daily_sea_ice_E_3_df <- c()
for (k in 2:length(daily_sea_ice_E_3)) {
  daily_SI <- daily_sea_ice_E_3[[k]]
  day_nbr <- seq(from = 1, to = length(daily_SI), by = 1)
  year <- rep(years[k], times = length(daily_SI))
  df_k <- cbind(year, day_nbr, daily_SI)
  daily_sea_ice_E_3_df <- rbind(daily_sea_ice_E_3_df, df_k)
}
daily_sea_ice_E_3_df <- as.data.frame(daily_sea_ice_E_3_df)

ggplot(daily_sea_ice_E_3_df, aes(x = day_nbr, y = daily_SI, group = year, color = year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of the year",
       y = "Daily sea ice concentration") +
  facet_wrap(~year)

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_sea_ice_E_buffer_pelagic_area.png",
       width = 10, height = 7.5)




# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load data
load("06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run function
SI_retreat_advance_date_E <- get_SI_retreat_advance_date(daily_sea_ice_E_3, 
                                                         years)
# Save
write_csv(SI_retreat_advance_date_E, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_E.csv")

# plot
SI_retreat_advance_date_E %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_E_buffer_pelagic_area.png",
       width = 6, height = 4)


# ~ 5. Plot approach E ---------------------------------------------------------

svalbard.planar <- readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map_planar")
buffer.size <- 500000
svalbard.buffers.planar <- gBuffer(svalbard.planar, # generate the buffer around Svalbard
                                   byid = TRUE, 
                                   id = NULL, 
                                   width = buffer.size, 
                                   quadsegs = 5, 
                                   capStyle = "ROUND",
                                   joinStyle = "ROUND",
                                   mitreLimit = 1.0)

sea_ice_stack <- raster::stack(x = "04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.tif")
svalbard.buffers <- spTransform(svalbard.buffers.planar, 
                                CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates
svalbard <- spTransform(svalbard.planar, 
                        CRSobj = crs(sea_ice_stack)) # Re project on a sphere, in polar coordinates


pelagic_bears_area_polar <- readOGR(dsn = "06_processed_data/sea_ice_data",
                                    layer = "pelagic_bears_area_polar")

study_area_E <- aggregate(rbind(svalbard.buffers, pelagic_bears_area_polar))

plot(study_area_E)
plot(svalbard,
     add = TRUE)

study_area_E_df <- tidy(study_area_E) 
svalbard_df <- tidy(svalbard) 


ggplot() +
  geom_polygon(data = svalbard_df, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 1) +
  geom_polygon(data = study_area_E_df, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 1) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  labs(x = "",
       y = "",
       fill = "Sea ice \ncover")

ggsave("10_meetings/2021-05-06 Meeting with Sarah/area E.svg")
# F. Buffer around Svalabrd with a 60% sea ice concentration threshold ============

# ~ 4. Calculate sea ice retreat and advance dates -----------------------------

source("05_script/environmental_covariates/sea_ice_metrics_calculation_functions.R")

# Load data
load("06_processed_data/sea_ice_data/daily_sea_ice_500km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

# Run function
SI_retreat_advance_date_F <- get_SI_retreat_advance_date_any_threshold(daily_sea_ice_B_3, 
                                                                       years,
                                                                       threshold = 60)


# Save
write_csv(SI_retreat_advance_date_F, 
          "06_processed_data/sea_ice_data/retreat_advance_ice_free_days_F.csv")

# plot
SI_retreat_advance_date_C %>%
  pivot_longer(cols = c("day_retreat", "day_advance")) %>%
  ggplot() +
  geom_point(aes(x = year, y = value, group = name, color = name)) +
  theme_bw() +
  labs(x = "Year", 
       y = "Day of the year",
       color = "")

ggsave("06_processed_data/sea_ice_data/graphs/retreat_advance_day_depth_300m_barents_sea.png",
       width = 6, height = 4)


# ~ 5. Plot approach B ---------------------------------------------------------

bathymetry_barents_sea <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped_Barents_sea_final.tif")

barents_sea_boundary_polar <- readOGR("06_processed_data/subpopulation_boundary",
                                      layer = "barents_sea_subpop_polar")

svalbard <- readOGR(dsn = "06_processed_data/svalbard_map",
                    layer = "svalbard_map_polar")

cuts <- c(0, -300, -1000, -2000, -3000, -6000) #set breaks
pal <- colorRampPalette(c("black", "grey99"))
plot(bathymetry_barents_sea,
     breaks = cuts, 
     col = pal(7))
plot(barents_sea_boundary_polar,
     add = TRUE)
plot(svalbard,
     add = TRUE)




# For walid
costline <- readOGR(dsn = "D:/temp_walid",
                           layer = "ne_10m_coastline")
crop.square <- extent(90, 170, -25, 30)
costline_south_east_asia <- crop(costline, crop.square)

plot(costline_south_east_asia)

long <- c()
lat <- c()
id <- c()
for (k in 1:length(costline_south_east_asia@lines)) {
  long <- c(long, costline_south_east_asia@lines[[k]]@Lines[[1]]@coords[, 1])
  lat <- c(lat, costline_south_east_asia@lines[[k]]@Lines[[1]]@coords[, 2])
  id <- c(id, 
          rep(costline_south_east_asia@lines[[k]]@ID,
          times = length(costline_south_east_asia@lines[[k]]@Lines[[1]]@coords[, 1])))
}
coords <- data.frame(id = id,
                     long = long,
                     lat = lat)
write_csv(coords, "D:/temp_walid/coords_costline_southeast_asia.csv")
