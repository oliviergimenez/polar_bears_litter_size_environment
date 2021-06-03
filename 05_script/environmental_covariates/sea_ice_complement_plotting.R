#==============================================================================#
#                                                                              #
#             Facetted plots of consecutive ice-free/covered days              #
#                             (+ other plots)                                  #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rgeos)
library(ff)
library(viridis)
Sys.setenv(LANG = "en")

# O. Functions =================================================================

# get_ice_free_days_in_row <- function(daily_sea_ice, days_in_a_row, area_considered) {
#   years <- seq(from = 1991, to = 2020, by = 1)
#   
#   means_for_max <- c()
#   start_days_for_max <- c()
#   years_for_max <- c()
#   df_max <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
#     
#     # Maximum sea ice extent
#     means_for_max_i <- c()
#     start_days_for_max_i <- c()
#     for (j in 1:100) {
#       means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
#       start_days_for_max_i <- c(start_days_for_max_i, days[j])
#     }
#     year_for_max_i <- rep(years[i], times = length(means_for_max_i))
#     
#     df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
#     max_year_i <- df_i[max(df_i$means_for_max), ]
#     df_max <- rbind(df_max, max_year_i)
#   }
#   start_day_max <- round(mean(df_max$start_days_for_max_i))
#   
#   # Find best window for min sea ice extent
#   means_for_min <- c()
#   start_days_for_min <- c()
#   years_for_min <- c()
#   df_min <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
#     
#     # Maximum sea ice extent
#     means_for_min_i <- c()
#     start_days_for_min_i <- c()
#     for (j in 200:280) {
#       means_for_min_i <- c(means_for_min_i,
#                            mean(daily_sea_ice[[i]][j:(j+30)]))
#       start_days_for_min_i <- c(start_days_for_min_i, 
#                                 days[j])
#     }
#     year_for_min_i <- rep(years[i], times = length(means_for_min_i))
#     
#     df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
#     min_year_i <- df_i[min(df_i$means_for_min), ]
#     df_min <- rbind(df_min, min_year_i)
#   }
#   start_day_min <- round(mean(df_min$start_days_for_min_i))
#   
#   
#   # Calculate the average yearly max and min sea ice extent
#   means_max <- c()
#   means_min <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     means_max <- c(means_max, 
#                    mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
#     means_min <- c(means_min, 
#                    mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
#   }
#   max_sea_ice <- mean(means_max)
#   min_sea_ice <- mean(means_min)
#   
#   # Calculate the halfway value
#   halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
#   print(halfway_sea_ice)
#   
#   # Retrive the day number of retreat/advance for each year
#   # RETREAT 
#   retreat_days <- c()
#   year <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     year <- c(year, years[i])
#     five_day_SI <- c()
#     for (j in 1:250) {
#       five_day_SI <- c( five_day_SI,
#                         sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
#       
#     }
#     j <- 250
#     while(five_day_SI[j] == 5) {
#       j <- j - 1
#     }
#     retreat_days <- c(retreat_days, j+1)
#   }
#   
#   # To plot (facetted plot)
#   n_day_SI_test <- c()
#   year <- c()
#   day_start <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     n_day_SI <- c()
#     year <- c(year, rep(years[i], times = 250))
#     day_start <- c(day_start,seq(1, 250, by = 1))
#     for (j in 1:250) {
#       n_day_SI <- c( n_day_SI,
#                         sum(daily_sea_ice[[i]][seq(from = j, to = (j + days_in_a_row - 1), 
#                                                    by = 1)] < halfway_sea_ice) )
#       
#     }
#     n_day_SI_test <- c(n_day_SI_test, n_day_SI)
#     
#   }
#   
#   df_plot_SI_retreat <- data.frame(year = year,
#                                    day_start = day_start,
#                                    value = n_day_SI_test)
#   
#   ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
#     geom_line() + 
#     theme_bw() +
#     facet_wrap(year) +
#     labs(x = "Start day of 5 days in a row",
#          y = "Number of ice-free days in a row")
#   
#   ggsave(filename = paste0("06_processed_data/sea_ice_data/graphs/",
#                            days_in_a_row, "_ice_free_days_in_a_row_",
#                            area_considered, ".png"),
#          width = 10, height = 7.5)
#   
# }


# get_ice_covered_days_in_row <- function(daily_sea_ice, days_in_a_row, area_considered) {
#   years <- seq(from = 1991, to = 2020, by = 1)
#   means_for_max <- c()
#   start_days_for_max <- c()
#   years_for_max <- c()
#   df_max <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
#     
#     # Maximum sea ice extent
#     means_for_max_i <- c()
#     start_days_for_max_i <- c()
#     for (j in 1:100) {
#       means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
#       start_days_for_max_i <- c(start_days_for_max_i, days[j])
#     }
#     year_for_max_i <- rep(years[i], times = length(means_for_max_i))
#     
#     df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
#     max_year_i <- df_i[max(df_i$means_for_max), ]
#     df_max <- rbind(df_max, max_year_i)
#   }
#   start_day_max <- round(mean(df_max$start_days_for_max_i))
#   
#   # Find best window for min sea ice extent
#   means_for_min <- c()
#   start_days_for_min <- c()
#   years_for_min <- c()
#   df_min <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
#     
#     # Maximum sea ice extent
#     means_for_min_i <- c()
#     start_days_for_min_i <- c()
#     for (j in 200:280) {
#       means_for_min_i <- c(means_for_min_i,
#                            mean(daily_sea_ice[[i]][j:(j+30)]))
#       start_days_for_min_i <- c(start_days_for_min_i, 
#                                 days[j])
#     }
#     year_for_min_i <- rep(years[i], times = length(means_for_min_i))
#     
#     df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
#     min_year_i <- df_i[min(df_i$means_for_min), ]
#     df_min <- rbind(df_min, min_year_i)
#   }
#   start_day_min <- round(mean(df_min$start_days_for_min_i))
#   
#   
#   # Calculate the average yearly max and min sea ice extent
#   means_max <- c()
#   means_min <- c()
#   for (i in 2:length(daily_sea_ice)) {
#     means_max <- c(means_max, 
#                    mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
#     means_min <- c(means_min, 
#                    mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
#   }
#   max_sea_ice <- mean(means_max)
#   min_sea_ice <- mean(means_min)
#   
#   # Calculate the halfway value
#   halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
#   print(halfway_sea_ice)
#   
#   # Retrive the day number of retreat/advance for each year
#   # ADVANCE
#   advance_days <- c()
#   year <- c()
#   for (i in 2:( length(daily_sea_ice) - 1)) {
#     daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
#                          daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
#     year <- c(year, years[i])
#     five_day_SI <- c()
#     for (j in 1:length(daily_sea_ice_i)) {
#       five_day_SI <- c( five_day_SI,
#                         sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
#       
#     }
#     j <- 1
#     while(five_day_SI[j] < 5) {
#       j <- j + 1
#     }
#     advance_days <- c(advance_days, j + 249)
#   }
# 
#   # To plot (facetted plot)
#   five_day_SI_test <- c()
#   year <- c()
#   day_start <- c()
#   for (i in 2:( length(daily_sea_ice) - 1 )) {
#     daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
#                          daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
#     five_day_SI <- c()
#     year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
#     day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
#     for (j in 1:length(daily_sea_ice_i)) {
#       five_day_SI <- c( five_day_SI,
#                         sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
#       
#     }
#     five_day_SI_test <- c(five_day_SI_test, five_day_SI)
#     
#   }
#   
#   df_plot_SI_advance <- data.frame(year = year,
#                                    day_start = day_start,
#                                    value = five_day_SI_test)
#   
#   ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
#     geom_line() + 
#     theme_bw() +
#     facet_wrap(year) +
#     labs(x = "Start day of 5 days in a row",
#          y = "Number of ice-covered days in a row")
#   
#   ggsave(filename = paste0("06_processed_data/sea_ice_data/graphs/",
#                            days_in_a_row, "_ice_covered_days_in_a_row_",
#                            area_considered, ".png"),
#          width = 10, height = 7.5)
#   
# }


# Load the "get_haldway_sea_ice_concentration" function
source("05_script/environmental_covariates/sea_ice_metrics_calc;ulation_functions.R")


get_ice_free_days_in_row <- function(daily_sea_ice, days_in_a_row, area_considered) {
  years <- seq(from = 1991, to = 2020, by = 1)
  
  # Calculate the halfway sea ice concentration
  halfway_sea_ice <- get_halfway_sea_ice_concentration(daily_sea_ice)
  
  # Determine whether each day satisfies the number-of-ice-free-days-in-a-row condition
  n_day_SI_all <- c()
  duration_all <- c()
  day_start_all <- c()
  year_all <- c()
  for (k in 1:length(days_in_a_row)) { # A loop for each number of days in a row: 3, 5, and 7
    n_day_SI_k <- c()
    year_k <- c()
    day_start_k <- c()
    day_start_i <- c()
    for (i in 2:length(daily_sea_ice)) { # A loop for each year
      n_day_SI_i <- c()
      year_k <- c(year_k, rep(years[i], times = 250))
      for (j in 1:250) {                 # A loop for each day of the year for which to test whether it satisfies the condition
        n_day_SI_j <- 0
        if (daily_sea_ice[[i]][j] < halfway_sea_ice) { # If the first day of the 3 days interval to be tested has a sea ice concentration lower than the threshold
          l <- 1
          cdt <- TRUE
          while (l <= (days_in_a_row[k]) && cdt == TRUE) {
            cdt <- daily_sea_ice[[i]][j + l] < halfway_sea_ice
            if (isTRUE(cdt)) {
              n_day_SI_j <- n_day_SI_j + 1
            }
            l <- l + 1
          }
        }
        n_day_SI_i <- c(n_day_SI_i, n_day_SI_j)
      }
      n_day_SI_k <- c(n_day_SI_k, n_day_SI_i)
      n_day_SI_k[n_day_SI_k < days_in_a_row[k]] <- 0
      day_start_k <- c(day_start_k, day_start_i)
    }
    n_day_SI_all <- c(n_day_SI_all, n_day_SI_k)
    duration_all <- c(duration_all, rep(days_in_a_row[k], times = length(n_day_SI_k)))
    day_start_all <- c(day_start_all, day_start_k)
    year_all <- c(year_all, year_k)
  }
  
  df_plot_SI_retreat <- data.frame(year = year_all,
                                   day_start = rep(x = seq(1, 250, 1), 
                                                   times = (length(days_in_a_row) * (length(years) - 1))),
                                   value = n_day_SI_all,
                                   days_in_a_row = duration_all)
  
  # Plot 
  ggplot(data = df_plot_SI_retreat, 
         aes(x = day_start, y = value, group = as.factor(days_in_a_row), color = as.factor(days_in_a_row))) +
    geom_line() +
    theme_bw() +
    scale_color_viridis(discrete = TRUE,
                        labels = c("3 day", "5 days", "7 days")) +
    labs(x = "day",
         y = "number of ice-free days in a row",
         color = "threshold") +
    facet_wrap(~ year)
  
  # Save
  ggsave(filename = paste0("06_processed_data/sea_ice_data/graphs/", "3_5_7_ice_free_days_in_a_row_",
                           area_considered, ".png"),
         width = 10, height = 7.5)
  
}



get_ice_covered_days_in_row <- function(daily_sea_ice, days_in_a_row, area_considered) {
  years <- seq(from = 1991, to = 2020, by = 1)
  
  # Calculate the halfway sea ice concentration
  halfway_sea_ice <- get_halfway_sea_ice_concentration(daily_sea_ice)
  
  # Determine whether each day satisfies the number-of-ice-covered-days-in-a-row condition
  n_day_SI_all <- c()
  duration_all <- c()
  day_start_all <- c()
  year_all <- c()
  for (k in 1:length(days_in_a_row)) {    # A loop for each number of days in a row: 3, 5, and 7
    n_day_SI_k <- c()
    year_k <- c()
    day_start_k <- c()
    for (i in 2:(length(daily_sea_ice) -1 )) {   # A loop for each year
      
      # Here we fuse days from year i to days from year i+1 because sea ice may advance in early year i+1
      daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                           daily_sea_ice[[i+1]][seq(from = 1, to = 140, by = 1)])
      n_day_SI_i <- c()
      year_k <- c(year_k, rep(years[i], times = 250))
      day_start_i <- seq(250, 249 + 250, by = 1)
      for (j in 1:250) {                          # A loop for each day of the year for which to test whether it satisfies the condition
        n_day_SI_j <- 0
        if (daily_sea_ice_i[j] > halfway_sea_ice) { # If the first day of the 3 (or 5 or 7) day interval to be tested
          l <- 1                                    # has a sea ice concentration > sea ice threshold
          cdt <- TRUE
          while (l <= (days_in_a_row[k]) && cdt == TRUE) {
            cdt <- daily_sea_ice_i[j + l] > halfway_sea_ice
            if (isTRUE(cdt)) {
              n_day_SI_j <- n_day_SI_j + 1
            }
            l <- l + 1
          }
        }
        
        n_day_SI_i <- c(n_day_SI_i, n_day_SI_j)
      }
      n_day_SI_k <- c(n_day_SI_k, n_day_SI_i)
      n_day_SI_k[n_day_SI_k < days_in_a_row[k]] <- 0
      day_start_k <- c(day_start_k, day_start_i)
    }
    n_day_SI_all <- c(n_day_SI_all, n_day_SI_k)
    duration_all <- c(duration_all, rep(days_in_a_row[k], times = length(n_day_SI_k)))
    day_start_all <- c(day_start_all, day_start_k)
    year_all <- c(year_all, year_k)
  }
  
  df_plot_SI_advance <- data.frame(year = year_all,
                                   day_start = day_start_all,
                                   #rep(x = seq(1, 250, 1), 
                                   #             times = (length(days_in_a_row) * (length(years) - 1))),
                                   value = n_day_SI_all,
                                   days_in_a_row = duration_all)
  
  # Plot
  ggplot(data = df_plot_SI_advance, 
         aes(x = day_start, y = value, group = as.factor(days_in_a_row), color = as.factor(days_in_a_row))) +
    geom_line() +
    theme_bw() +
    scale_color_viridis(discrete = TRUE,
                        labels = c("3 day", "5 days", "7 days")) +
    labs(x = "day",
         y = "number of ice-covered days in a row",
         color = "threshold") +
    facet_wrap(~ year)
  
  # Save
  ggsave(filename = paste0("06_processed_data/sea_ice_data/graphs/", "3_5_7_ice_covered_days_in_a_row_",
                           area_considered, ".png"),
         width = 10, height = 7.5)
  
}






# Area considered: A (Barents Sea)
load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated_corrected.RData")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                         days_in_a_row = 3, area_considered = "A")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                         days_in_a_row = 5, area_considered = "A")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                         days_in_a_row = 7, area_considered = "A")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                            days_in_a_row = 3, area_considered = "A")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                            days_in_a_row = 5, area_considered = "A")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_A_3,
                            days_in_a_row = 7, area_considered = "A")

# Area considered: B (100km buffer around Svalbard)
load("06_processed_data/sea_ice_data/daily_sea_ice_100km_buffer_interpolated_corrected.RData")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                         days_in_a_row = 3, area_considered = "B")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                         days_in_a_row = 5, area_considered = "B")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                         days_in_a_row = 7, area_considered = "B")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                            days_in_a_row = 3, area_considered = "B")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                            days_in_a_row = 5, area_considered = "B")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_B_3,
                            days_in_a_row = 7, area_considered = "B")

# Area considered: C (<600m depth Svalabrd)
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_barents_sea_interpolated_corrected.RData")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                         days_in_a_row = 3, area_considered = "C")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                         days_in_a_row = 5, area_considered = "C")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                         days_in_a_row = 7, area_considered = "C")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                            days_in_a_row = 3, area_considered = "C")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                            days_in_a_row = 5, area_considered = "C")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_C_3,
                            days_in_a_row = 7, area_considered = "C")


# Area considered: D (<600m depth in 100km buffer around Svalbard)
load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                         days_in_a_row = 3, area_considered = "D")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                         days_in_a_row = 5, area_considered = "D")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                         days_in_a_row = 7, area_considered = "D")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                            days_in_a_row = 3, area_considered = "D")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                            days_in_a_row = 5, area_considered = "D")
get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_D_3,
                            days_in_a_row = 7, area_considered = "D")

# Area considered: E (100km buffer around Svalbard + pelagic bear area)
load("06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated_corrected.RData")
get_ice_free_days_in_row(daily_sea_ice = daily_sea_ice_E_3,
                         days_in_a_row = c(3, 5, 7), area_considered = "E")

get_ice_covered_days_in_row(daily_sea_ice = daily_sea_ice_E_3,
                            days_in_a_row = c(3, 5, 7), area_considered = "E")


# A. Whole Barents Sea =========================================================

load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_A_3 


# + 1. Facetted plots with consecutive ice-free/covered days -------------------

# Find best window for max sea ice extent
means_for_max <- c()
start_days_for_max <- c()
years_for_max <- c()
df_max <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_max_i <- c()
  start_days_for_max_i <- c()
  for (j in 1:100) {
    means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+29)]))
    start_days_for_max_i <- c(start_days_for_max_i, days[j])
  }
  year_for_max_i <- rep(years[i], times = length(means_for_max_i))
  
  df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
  max_year_i <- df_i[max(df_i$means_for_max), ]
  df_max <- rbind(df_max, max_year_i)
}
start_day_max <- round(mean(df_max$start_days_for_max_i))

# Find best window for min sea ice extent
means_for_min <- c()
start_days_for_min <- c()
years_for_min <- c()
df_min <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_min_i <- c()
  start_days_for_min_i <- c()
  for (j in 200:280) {
    means_for_min_i <- c(means_for_min_i,
                         mean(daily_sea_ice[[i]][j:(j+29)]))
    start_days_for_min_i <- c(start_days_for_min_i, 
                              days[j])
  }
  year_for_min_i <- rep(years[i], times = length(means_for_min_i))
  
  df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
  min_year_i <- df_i[min(df_i$means_for_min), ]
  df_min <- rbind(df_min, min_year_i)
}
start_day_min <- round(mean(df_min$start_days_for_min_i))


# Calculate the average yearly max and min sea ice extent
means_max <- c()
means_min <- c()
for (i in 2:length(daily_sea_ice)) {
  means_max <- c(means_max, 
                 mean(daily_sea_ice[[i]][start_day_max:(start_day_max+29)]))
  means_min <- c(means_min, 
                 mean(daily_sea_ice[[i]][start_day_min:(start_day_min+29)]))
}
max_sea_ice <- mean(means_max)
min_sea_ice <- mean(means_min)

# Calculate the halfway value
halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
print(halfway_sea_ice)

# Retrive the day number of retreat/advance for each year
# RETREAT 
retreat_days <- c()
year <- c()
for (i in 2:length(daily_sea_ice)) {
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  j <- 250
  while(five_day_SI[j] == 5) {
    j <- j - 1
  }
  retreat_days <- c(retreat_days, j+1)
}
retreat <- data.frame(year = year,
                      day_retreat = retreat_days)
# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:length(daily_sea_ice)) {
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = 250))
  day_start <- c(day_start,seq(1, 250, by = 1))
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)

}

df_plot_SI_retreat <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-free days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_free_days_in_a_row_A_barents_sea.png",
       width = 10, height = 7.5)

# ADVANCE
advance_days <- c()
year <- c()
for (i in 2:( length(daily_sea_ice) - 1)) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  j <- 1
  while(five_day_SI[j] < 5) {
    j <- j + 1
  }
  advance_days <- c(advance_days, j + 249)
}
advance <- data.frame(year = year,
                      day_advance = advance_days)

ggplot(data = advance, aes(x = year, y = day_advance)) +
  geom_point() +
  geom_smooth(method = lm)


# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:( length(daily_sea_ice) - 1 )) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
  day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_advance <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-covered days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_covered_days_in_a_row_A_barents_sea.png",
       width = 10, height = 7.5)





# + 2. Calculate the 30-year median/mean sea ice concentration for each day ----

load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

daily_sea_ice_df <- c()
for (i in 2:length(daily_sea_ice)) {
  daily_SI_year_i <- daily_sea_ice[[i]]
  if (length(daily_SI_year_i) == 365) {
    daily_SI_year_i <- c(daily_SI_year_i[1:59], NA, daily_SI_year_i[60:365])
  }
  daily_sea_ice_df <- cbind(daily_sea_ice_df, daily_SI_year_i)
}
daily_sea_ice_df <- as.data.frame(daily_sea_ice_df) 
colnames(daily_sea_ice_df) <- years[-1]

mean_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                    FUN = mean, na.rm = TRUE)
median_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                      FUN = median, na.rm = TRUE)

daily_sea_ice_df <- daily_sea_ice_df %>%
  mutate(mean_daily = mean_daily,
         median_daily = median_daily,
         day_nbr = seq(1, 366, by = 1))


write_csv(daily_sea_ice_df, "06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_df_30-year_mean.csv")
daily_sea_ice_df <- read_csv("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_df_30-year_mean.csv")

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = mean_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/mean_daily_sea_ice_A_barents_sea.png",
       width = 6, height = 4)

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = median_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/median_daily_sea_ice_A_barents_sea.png",
       width = 6, height = 4)




# + 3. Plot the 2-day variation in SI concentration to spot mistakes -----------


load("06_processed_data/sea_ice_data/daily_sea_ice_barents_sea_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)


delta <- c()
year <- c()
day_nbr <- c()
for (i in 1:length(daily_sea_ice_A_2)) {
  year_i <- rep(years[i], times = length(daily_sea_ice_A_2[[i]]) -1)
  day_nbr_i <- seq(1, length(daily_sea_ice_A_2[[i]]) -1)
  delta_i <- c()
  for (j in 1:( length(daily_sea_ice_A_2[[i]]) - 1 )) {
    #delta_i <- c(delta_i, abs(daily_sea_ice_A_2[[i]][j] - daily_sea_ice_A_2[[i]][j + 1]))
    delta_i <- c(delta_i, daily_sea_ice_A_2[[i]][j] - daily_sea_ice_A_2[[i]][j + 1])
    
  }
  year <- c(year, year_i)
  delta <- c(delta, delta_i)
  day_nbr <- c(day_nbr, day_nbr_i)
}
df_delta <- data.frame(day_nbr = day_nbr,
                       year = year,
                       delta = delta)


ggplot(df_delta, aes(x = day_nbr, y = delta)) +
  geom_line() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day",
       y = "2-day difference in sea ice concentration")



ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_differences_SI_barents_sea.png",
       width = 10, height = 7.5)




# B. Buffer around Svalbard ====================================================

load("06_processed_data/sea_ice_data/daily_sea_ice_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_B_3 


# + 1. Facetted plots with consecutive ice-free/covered days -------------------

# Find best window for max sea ice extent
means_for_max <- c()
start_days_for_max <- c()
years_for_max <- c()
df_max <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_max_i <- c()
  start_days_for_max_i <- c()
  for (j in 1:100) {
    means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_max_i <- c(start_days_for_max_i, days[j])
  }
  year_for_max_i <- rep(years[i], times = length(means_for_max_i))
  
  df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
  max_year_i <- df_i[max(df_i$means_for_max), ]
  df_max <- rbind(df_max, max_year_i)
}
start_day_max <- round(mean(df_max$start_days_for_max_i))

# Find best window for min sea ice extent
means_for_min <- c()
start_days_for_min <- c()
years_for_min <- c()
df_min <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_min_i <- c()
  start_days_for_min_i <- c()
  for (j in 200:280) {
    means_for_min_i <- c(means_for_min_i,
                         mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_min_i <- c(start_days_for_min_i, 
                              days[j])
  }
  year_for_min_i <- rep(years[i], times = length(means_for_min_i))
  
  df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
  min_year_i <- df_i[min(df_i$means_for_min), ]
  df_min <- rbind(df_min, min_year_i)
}
start_day_min <- round(mean(df_min$start_days_for_min_i))


# Calculate the average yearly max and min sea ice extent
means_max <- c()
means_min <- c()
for (i in 2:length(daily_sea_ice)) {
  means_max <- c(means_max, 
                 mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
  means_min <- c(means_min, 
                 mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
}
max_sea_ice <- mean(means_max)
min_sea_ice <- mean(means_min)

# Calculate the halfway value
halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
print(halfway_sea_ice)

# Retrive the day number of retreat/advance for each year
# RETREAT 
retreat_days <- c()
year <- c()
for (i in 2:length(daily_sea_ice)) {
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  j <- 250
  while(five_day_SI[j] == 5) {
    j <- j - 1
  }
  retreat_days <- c(retreat_days, j+1)
}
retreat <- data.frame(year = year,
                      day_retreat = retreat_days)
# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:length(daily_sea_ice)) {
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = 250))
  day_start <- c(day_start,seq(1, 250, by = 1))
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_retreat <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-free days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_free_days_in_a_row_B_100km_buffer.png",
       width = 10, height = 7.5)

# ADVANCE
advance_days <- c()
year <- c()
for (i in 2:( length(daily_sea_ice) - 1)) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  j <- 1
  while(five_day_SI[j] < 5) {
    j <- j + 1
  }
  advance_days <- c(advance_days, j + 249)
}
advance <- data.frame(year = year,
                      day_advance = advance_days)

ggplot(data = advance, aes(x = year, y = day_advance)) +
  geom_point() +
  geom_smooth(method = lm)


# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:( length(daily_sea_ice) - 1 )) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
  day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_advance <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-covered days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_covered_days_in_a_row_B_100km_buffer.png",
       width = 10, height = 7.5)





# + 2. Calculate the 30-year median/mean sea ice concentration for each day ----

load("06_processed_data/sea_ice_data/daily_sea_ice_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)

daily_sea_ice_df <- c()
for (i in 2:length(daily_sea_ice)) {
  daily_SI_year_i <- daily_sea_ice[[i]]
  if (length(daily_SI_year_i) == 365) {
    daily_SI_year_i <- c(daily_SI_year_i[1:59], NA, daily_SI_year_i[60:365])
  }
  daily_sea_ice_df <- cbind(daily_sea_ice_df, daily_SI_year_i)
}
daily_sea_ice_df <- as.data.frame(daily_sea_ice_df) 
colnames(daily_sea_ice_df) <- years[-1]

mean_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                    FUN = mean, na.rm = TRUE)
median_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                      FUN = median, na.rm = TRUE)

daily_sea_ice_df <- daily_sea_ice_df %>%
  mutate(mean_daily = mean_daily,
         median_daily = median_daily,
         day_nbr = seq(1, 366, by = 1))


write_csv(daily_sea_ice_df, "06_processed_data/sea_ice_data/daily_sea_ice_50km_buffer_df_30-year_mean.csv")
daily_sea_ice_df <- read_csv("06_processed_data/sea_ice_data/daily_sea_ice_50km_buffer_df_30-year_mean.csv")

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = mean_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/mean_daily_sea_ice_B_100km_buffer.png",
       width = 6, height = 4)

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = median_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/median_daily_sea_ice_B_100km_buffer.png",
       width = 6, height = 4)




# + 3. Plot the 2-day variation in SI concentration to spot mistakes -----------


load("06_processed_data/sea_ice_data/daily_sea_ice_100km_buffer_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)


delta <- c()
year <- c()
day_nbr <- c()
for (i in 1:length(daily_sea_ice_B_2)) {
  year_i <- rep(years[i], times = length(daily_sea_ice_B_2[[i]]) -1)
  day_nbr_i <- seq(1, length(daily_sea_ice_B_2[[i]]) -1)
  delta_i <- c()
  for (j in 1:( length(daily_sea_ice_B_2[[i]]) - 1 )) {
    #delta_i <- c(delta_i, abs(daily_sea_ice_B_2[[i]][j] - daily_sea_ice_B_2[[i]][j + 1]))
    delta_i <- c(delta_i, daily_sea_ice_B_2[[i]][j] - daily_sea_ice_B_2[[i]][j + 1])
    
  }
  year <- c(year, year_i)
  delta <- c(delta, delta_i)
  day_nbr <- c(day_nbr, day_nbr_i)
}
df_delta <- data.frame(day_nbr = day_nbr,
                       year = year,
                       delta = delta)


ggplot(df_delta, aes(x = day_nbr, y = delta)) +
  geom_line() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day",
       y = "2-day difference in sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_differences_SI_100_km_buffer.png",
       width = 10, height = 7.5)



# C. <600m depth Barents Sea ===================================================

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_C_3 


# + 1. Facetted plots with consecutive ice-free/covered days -------------------

# Find best window for max sea ice extent
means_for_max <- c()
start_days_for_max <- c()
years_for_max <- c()
df_max <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_max_i <- c()
  start_days_for_max_i <- c()
  for (j in 1:100) {
    means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_max_i <- c(start_days_for_max_i, days[j])
  }
  year_for_max_i <- rep(years[i], times = length(means_for_max_i))
  
  df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
  max_year_i <- df_i[max(df_i$means_for_max), ]
  df_max <- rbind(df_max, max_year_i)
}
start_day_max <- round(mean(df_max$start_days_for_max_i))

# Find best window for min sea ice extent
means_for_min <- c()
start_days_for_min <- c()
years_for_min <- c()
df_min <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_min_i <- c()
  start_days_for_min_i <- c()
  for (j in 200:280) {
    means_for_min_i <- c(means_for_min_i,
                         mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_min_i <- c(start_days_for_min_i, 
                              days[j])
  }
  year_for_min_i <- rep(years[i], times = length(means_for_min_i))
  
  df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
  min_year_i <- df_i[min(df_i$means_for_min), ]
  df_min <- rbind(df_min, min_year_i)
}
start_day_min <- round(mean(df_min$start_days_for_min_i))


# Calculate the average yearly max and min sea ice extent
means_max <- c()
means_min <- c()
for (i in 2:length(daily_sea_ice)) {
  means_max <- c(means_max, 
                 mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
  means_min <- c(means_min, 
                 mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
}
max_sea_ice <- mean(means_max)
min_sea_ice <- mean(means_min)

# Calculate the halfway value
halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
print(halfway_sea_ice)

# Retrive the day number of retreat/advance for each year
# RETREAT 
retreat_days <- c()
year <- c()
for (i in 2:length(daily_sea_ice)) {
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  j <- 250
  while(five_day_SI[j] == 5) {
    j <- j - 1
  }
  retreat_days <- c(retreat_days, j+1)
}
retreat <- data.frame(year = year,
                      day_retreat = retreat_days)
# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:length(daily_sea_ice)) {
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = 250))
  day_start <- c(day_start,seq(1, 250, by = 1))
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_retreat <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-free days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_free_days_in_a_row_C_600m_depth_barents_sea.png",
       width = 10, height = 7.5)

# ADVANCE
advance_days <- c()
year <- c()
for (i in 2:( length(daily_sea_ice) - 1)) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  j <- 1
  while(five_day_SI[j] < 5) {
    j <- j + 1
  }
  advance_days <- c(advance_days, j + 249)
}
advance <- data.frame(year = year,
                      day_advance = advance_days)

ggplot(data = advance, aes(x = year, y = day_advance)) +
  geom_point() +
  geom_smooth(method = lm)


# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:( length(daily_sea_ice) - 1 )) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
  day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_advance <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-covered days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_covered_days_in_a_row_C_600m_depth_barents_sea.png",
       width = 10, height = 7.5)





# + 2. Calculate the 30-year median/mean sea ice concentration for each day ----

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_barents_sea_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_C_3 

daily_sea_ice_df <- c()
for (i in 2:length(daily_sea_ice)) {
  daily_SI_year_i <- daily_sea_ice[[i]]
  if (length(daily_SI_year_i) == 365) {
    daily_SI_year_i <- c(daily_SI_year_i[1:59], NA, daily_SI_year_i[60:365])
  }
  daily_sea_ice_df <- cbind(daily_sea_ice_df, daily_SI_year_i)
}
daily_sea_ice_df <- as.data.frame(daily_sea_ice_df) 
colnames(daily_sea_ice_df) <- years[-1]

mean_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                    FUN = mean, na.rm = TRUE)
median_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                      FUN = median, na.rm = TRUE)

daily_sea_ice_df <- daily_sea_ice_df %>%
  mutate(mean_daily = mean_daily,
         median_daily = median_daily,
         day_nbr = seq(1, 366, by = 1))


write_csv(daily_sea_ice_df, "06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_barents_sea_df_30-year_mean.csv")
daily_sea_ice_df <- read_csv("06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_barents_sea_df_30-year_mean.csv")

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = mean_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/mean_daily_sea_ice_C_600m_depth_barents_sea.png",
       width = 6, height = 4)

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = median_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/median_daily_sea_ice_C_600m_depth_barents_sea.png",
       width = 6, height = 4)




# + 3. Plot the 2-day variation in SI concentration to spot mistakes -----------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_barents_sea_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_C_2 


delta <- c()
year <- c()
day_nbr <- c()
for (i in 1:length(daily_sea_ice)) {
  year_i <- rep(years[i], times = length(daily_sea_ice[[i]]) -1)
  day_nbr_i <- seq(1, length(daily_sea_ice[[i]]) -1)
  delta_i <- c()
  for (j in 1:( length(daily_sea_ice[[i]]) - 1 )) {
    #delta_i <- c(delta_i, abs(daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1]))
    delta_i <- c(delta_i, daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1])
    
  }
  year <- c(year, year_i)
  delta <- c(delta, delta_i)
  day_nbr <- c(day_nbr, day_nbr_i)
}
df_delta <- data.frame(day_nbr = day_nbr,
                       year = year,
                       delta = delta)


ggplot(df_delta, aes(x = day_nbr, y = delta)) +
  geom_line() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day",
       y = "2-day difference in sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_differences_SI_600m_depth_barents_sea.png",
       width = 10, height = 7.5)






# D. <600m depth buffer around Svalbard ========================================


load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_D_3 


# + 1. Facetted plots with consecutive ice-free/covered days -------------------

# Find best window for max sea ice extent
means_for_max <- c()
start_days_for_max <- c()
years_for_max <- c()
df_max <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_max_i <- c()
  start_days_for_max_i <- c()
  for (j in 1:100) {
    means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_max_i <- c(start_days_for_max_i, days[j])
  }
  year_for_max_i <- rep(years[i], times = length(means_for_max_i))
  
  df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
  max_year_i <- df_i[max(df_i$means_for_max), ]
  df_max <- rbind(df_max, max_year_i)
}
start_day_max <- round(mean(df_max$start_days_for_max_i))

# Find best window for min sea ice extent
means_for_min <- c()
start_days_for_min <- c()
years_for_min <- c()
df_min <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_min_i <- c()
  start_days_for_min_i <- c()
  for (j in 200:280) {
    means_for_min_i <- c(means_for_min_i,
                         mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_min_i <- c(start_days_for_min_i, 
                              days[j])
  }
  year_for_min_i <- rep(years[i], times = length(means_for_min_i))
  
  df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
  min_year_i <- df_i[min(df_i$means_for_min), ]
  df_min <- rbind(df_min, min_year_i)
}
start_day_min <- round(mean(df_min$start_days_for_min_i))


# Calculate the average yearly max and min sea ice extent
means_max <- c()
means_min <- c()
for (i in 2:length(daily_sea_ice)) {
  means_max <- c(means_max, 
                 mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
  means_min <- c(means_min, 
                 mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
}
max_sea_ice <- mean(means_max)
min_sea_ice <- mean(means_min)

# Calculate the halfway value
halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
print(halfway_sea_ice)

# Retrive the day number of retreat/advance for each year
# RETREAT 
retreat_days <- c()
year <- c()
for (i in 2:length(daily_sea_ice)) {
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  j <- 250
  while(five_day_SI[j] == 5) {
    j <- j - 1
  }
  retreat_days <- c(retreat_days, j+1)
}
retreat <- data.frame(year = year,
                      day_retreat = retreat_days)
# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:length(daily_sea_ice)) {
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = 250))
  day_start <- c(day_start,seq(1, 250, by = 1))
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_retreat <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-free days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_free_days_in_a_row_D_600m_depth_100km_buffer.png",
       width = 10, height = 7.5)

# ADVANCE
advance_days <- c()
year <- c()
for (i in 2:( length(daily_sea_ice) - 1)) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  j <- 1
  while(five_day_SI[j] < 5) {
    j <- j + 1
  }
  advance_days <- c(advance_days, j + 249)
}
advance <- data.frame(year = year,
                      day_advance = advance_days)

ggplot(data = advance, aes(x = year, y = day_advance)) +
  geom_point() +
  geom_smooth(method = lm)


# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:( length(daily_sea_ice) - 1 )) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
  day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_advance <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-covered days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_covered_days_in_a_row_D_600m_depth_100km_buffer.png",
       width = 10, height = 7.5)





# + 2. Calculate the 30-year median/mean sea ice concentration for each day ----

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_D_3 

daily_sea_ice_df <- c()
for (i in 2:length(daily_sea_ice)) {
  daily_SI_year_i <- daily_sea_ice[[i]]
  if (length(daily_SI_year_i) == 365) {
    daily_SI_year_i <- c(daily_SI_year_i[1:59], NA, daily_SI_year_i[60:365])
  }
  daily_sea_ice_df <- cbind(daily_sea_ice_df, daily_SI_year_i)
}
daily_sea_ice_df <- as.data.frame(daily_sea_ice_df) 
colnames(daily_sea_ice_df) <- years[-1]

mean_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                    FUN = mean, na.rm = TRUE)
median_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                      FUN = median, na.rm = TRUE)

daily_sea_ice_df <- daily_sea_ice_df %>%
  mutate(mean_daily = mean_daily,
         median_daily = median_daily,
         day_nbr = seq(1, 366, by = 1))


write_csv(daily_sea_ice_df, "06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_100km_buffer_df_30-year_mean.csv")
daily_sea_ice_df <- read_csv("06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_100km_buffer_df_30-year_mean.csv")

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = mean_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/mean_daily_sea_ice_D_600m_depth_100km_buffer.png",
       width = 6, height = 4)

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = median_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/median_daily_sea_ice_D_600m_depth_100km_buffer.png",
       width = 6, height = 4)




# + 3. Plot the 2-day variation in SI concentration to spot mistakes -----------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_C_2 


delta <- c()
year <- c()
day_nbr <- c()
for (i in 1:length(daily_sea_ice)) {
  year_i <- rep(years[i], times = length(daily_sea_ice[[i]]) -1)
  day_nbr_i <- seq(1, length(daily_sea_ice[[i]]) -1)
  delta_i <- c()
  for (j in 1:( length(daily_sea_ice[[i]]) - 1 )) {
    #delta_i <- c(delta_i, abs(daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1]))
    delta_i <- c(delta_i, daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1])
    
  }
  year <- c(year, year_i)
  delta <- c(delta, delta_i)
  day_nbr <- c(day_nbr, day_nbr_i)
}
df_delta <- data.frame(day_nbr = day_nbr,
                       year = year,
                       delta = delta)


ggplot(df_delta, aes(x = day_nbr, y = delta)) +
  geom_line() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day",
       y = "2-day difference in sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_differences_D.png",
       width = 10, height = 7.5)



# E. Buffer around Svalbard + pelagic bears ====================================


load("06_processed_data/sea_ice_data/daily_sea_ice_E_buffer_pelagic_area_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_E_3 


# + 1. Facetted plots with consecutive ice-free/covered days -------------------

# Find best window for max sea ice extent
means_for_max <- c()
start_days_for_max <- c()
years_for_max <- c()
df_max <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_max_i <- c()
  start_days_for_max_i <- c()
  for (j in 1:100) {
    means_for_max_i <- c(means_for_max_i, mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_max_i <- c(start_days_for_max_i, days[j])
  }
  year_for_max_i <- rep(years[i], times = length(means_for_max_i))
  
  df_i <- as.data.frame(cbind(year_for_max_i, start_days_for_max_i, means_for_max_i))
  max_year_i <- df_i[max(df_i$means_for_max), ]
  df_max <- rbind(df_max, max_year_i)
}
start_day_max <- round(mean(df_max$start_days_for_max_i))

# Find best window for min sea ice extent
means_for_min <- c()
start_days_for_min <- c()
years_for_min <- c()
df_min <- c()
for (i in 2:length(daily_sea_ice)) {
  days <- seq(from = 1, to = length(daily_sea_ice[[i]]), by = 1)
  
  # Maximum sea ice extent
  means_for_min_i <- c()
  start_days_for_min_i <- c()
  for (j in 200:280) {
    means_for_min_i <- c(means_for_min_i,
                         mean(daily_sea_ice[[i]][j:(j+30)]))
    start_days_for_min_i <- c(start_days_for_min_i, 
                              days[j])
  }
  year_for_min_i <- rep(years[i], times = length(means_for_min_i))
  
  df_i <- as.data.frame(cbind(year_for_min_i, start_days_for_min_i, means_for_min_i))
  min_year_i <- df_i[min(df_i$means_for_min), ]
  df_min <- rbind(df_min, min_year_i)
}
start_day_min <- round(mean(df_min$start_days_for_min_i))


# Calculate the average yearly max and min sea ice extent
means_max <- c()
means_min <- c()
for (i in 2:length(daily_sea_ice)) {
  means_max <- c(means_max, 
                 mean(daily_sea_ice[[i]][start_day_max:(start_day_max+30)]))
  means_min <- c(means_min, 
                 mean(daily_sea_ice[[i]][start_day_min:(start_day_min+30)]))
}
max_sea_ice <- mean(means_max)
min_sea_ice <- mean(means_min)

# Calculate the halfway value
halfway_sea_ice <- median(c(max_sea_ice, min_sea_ice))
print(halfway_sea_ice)

# Retrive the day number of retreat/advance for each year
# RETREAT 
retreat_days <- c()
year <- c()
for (i in 2:length(daily_sea_ice)) {
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  j <- 250
  while(five_day_SI[j] == 5) {
    j <- j - 1
  }
  retreat_days <- c(retreat_days, j+1)
}
retreat <- data.frame(year = year,
                      day_retreat = retreat_days)
# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:length(daily_sea_ice)) {
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = 250))
  day_start <- c(day_start,seq(1, 250, by = 1))
  for (j in 1:250) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice[[i]][seq(from = j, to = j+4, by = 1)] < halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_retreat <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_retreat, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-free days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_free_days_in_a_row_E_100km_buffer_pelagic_area.png",
       width = 10, height = 7.5)

# ADVANCE
advance_days <- c()
year <- c()
for (i in 2:( length(daily_sea_ice) - 1)) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  year <- c(year, years[i])
  five_day_SI <- c()
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  j <- 1
  while(five_day_SI[j] < 5) {
    j <- j + 1
  }
  advance_days <- c(advance_days, j + 249)
}
advance <- data.frame(year = year,
                      day_advance = advance_days)

ggplot(data = advance, aes(x = year, y = day_advance)) +
  geom_point() +
  geom_smooth(method = lm)


# To plot (facetted plot)
five_day_SI_test <- c()
year <- c()
day_start <- c()
for (i in 2:( length(daily_sea_ice) - 1 )) {
  daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                       daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
  five_day_SI <- c()
  year <- c(year, rep(years[i], times = length(daily_sea_ice_i)))
  day_start <- c(day_start, seq(250, 249 + length(daily_sea_ice_i), by = 1))
  for (j in 1:length(daily_sea_ice_i)) {
    five_day_SI <- c( five_day_SI,
                      sum(daily_sea_ice_i[seq(from = j, to = j+4, by = 1)] > halfway_sea_ice) )
    
  }
  five_day_SI_test <- c(five_day_SI_test, five_day_SI)
  
}

df_plot_SI_advance <- data.frame(year = year,
                                 day_start = day_start,
                                 value = five_day_SI_test)

ggplot(df_plot_SI_advance, aes(x = day_start, y = value)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(year) +
  labs(x = "Start day of 5 days in a row",
       y = "Number of ice-covered days in a row")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_five_ice_covered_days_in_a_row_E_100km_buffer_pelagic_area.png",
       width = 10, height = 7.5)





# + 2. Calculate the 30-year median/mean sea ice concentration for each day ----

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_100km_buffer_interpolated_corrected.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_E_3 

daily_sea_ice_df <- c()
for (i in 2:length(daily_sea_ice)) {
  daily_SI_year_i <- daily_sea_ice[[i]]
  if (length(daily_SI_year_i) == 365) {
    daily_SI_year_i <- c(daily_SI_year_i[1:59], NA, daily_SI_year_i[60:365])
  }
  daily_sea_ice_df <- cbind(daily_sea_ice_df, daily_SI_year_i)
}
daily_sea_ice_df <- as.data.frame(daily_sea_ice_df) 
colnames(daily_sea_ice_df) <- years[-1]

mean_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                    FUN = mean, na.rm = TRUE)
median_daily <- apply(X = daily_sea_ice_df, MARGIN = 1, 
                      FUN = median, na.rm = TRUE)

daily_sea_ice_df <- daily_sea_ice_df %>%
  mutate(mean_daily = mean_daily,
         median_daily = median_daily,
         day_nbr = seq(1, 366, by = 1))


write_csv(daily_sea_ice_df, "06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_100km_buffer_df_30-year_mean.csv")
daily_sea_ice_df <- read_csv("06_processed_data/sea_ice_data/daily_sea_ice_600m_depth_100km_buffer_df_30-year_mean.csv")

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = mean_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/mean_daily_sea_ice_E_100km_buffer_pelagic_area.png",
       width = 6, height = 4)

ggplot(daily_sea_ice_df, aes(x = day_nbr, y = median_daily)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day",
       y = "Sea ice concentration")
ggsave("06_processed_data/sea_ice_data/graphs/median_daily_sea_ice_E_100km_buffer_pelagic_area.png",
       width = 6, height = 4)




# + 3. Plot the 2-day variation in SI concentration to spot mistakes -----------

load("06_processed_data/sea_ice_data/daily_sea_ice_depth_600m_barents_sea_interpolated.RData")
years <- seq(from = 1991, to = 2020, by = 1)
daily_sea_ice <- daily_sea_ice_C_2 


delta <- c()
year <- c()
day_nbr <- c()
for (i in 1:length(daily_sea_ice)) {
  year_i <- rep(years[i], times = length(daily_sea_ice[[i]]) -1)
  day_nbr_i <- seq(1, length(daily_sea_ice[[i]]) -1)
  delta_i <- c()
  for (j in 1:( length(daily_sea_ice[[i]]) - 1 )) {
    #delta_i <- c(delta_i, abs(daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1]))
    delta_i <- c(delta_i, daily_sea_ice[[i]][j] - daily_sea_ice[[i]][j + 1])
    
  }
  year <- c(year, year_i)
  delta <- c(delta, delta_i)
  day_nbr <- c(day_nbr, day_nbr_i)
}
df_delta <- data.frame(day_nbr = day_nbr,
                       year = year,
                       delta = delta)


ggplot(df_delta, aes(x = day_nbr, y = delta)) +
  geom_line() +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day",
       y = "2-day difference in sea ice concentration")

ggsave("06_processed_data/sea_ice_data/graphs/facetted_daily_differences_SI_600m_depth_barents_sea.png",
       width = 10, height = 7.5)

