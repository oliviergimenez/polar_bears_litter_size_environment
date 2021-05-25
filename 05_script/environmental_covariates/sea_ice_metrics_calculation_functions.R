#==============================================================================#
#                                                                              #
#      Functions for the calculation of sea ice metrics used in the model      #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rgeos)
library(ff)
Sys.setenv(LANG = "en")


interpolate <- function(daily_sea_ice, years) {
  index_j_NA <- c()
  index <- c()
  years.index <- c()
  for (i in 2:length(years)) {
    x <- which(is.na(daily_sea_ice[[i]]))
    if (length(x) > 0) {
      year <- rep(years[i], times = length(x))
      years.index <- c(years.index,
                       year)
      index_i <- rep(i, times = length(x))
      index <- c(index, index_i)
      index_j_NA <- c(index_j_NA,
                      x)
    }
  }
  missing_values <- data.frame(year = years.index,
                               index = index,
                               days_missing = index_j_NA)
  # Replace the missing values by the median/mean between the previous and following available point
  years_missing <- unique(missing_values$year)
  index_missing <- unique(missing_values$index)
  for (i in 1:length(unique(years.index))) {
    missing_values_i <- missing_values %>%
      filter(year == years_missing[i])
    
    if (nrow(missing_values_i) == 1) { # If a single value is missing
      daily_sea_ice[[ index_missing[i] ]][missing_values_i$days_missing] <-
        median(daily_sea_ice[[ index_missing[i] ]][missing_values_i$days_missing - 1],
               daily_sea_ice[[ index_missing[i] ]][missing_values_i$days_missing + 1])
      
    } else {
      sequencial_days <- c()
      for (j in 1:(nrow(missing_values_i)- 1)) {
        sequencial_days <- c(sequencial_days, missing_values_i$days_missing[j + 1] == missing_values_i$days_missing[j] + 1)
        
      }
      if (sum(sequencial_days) == (nrow(missing_values_i)- 1)) { # If all the missing days are one after the other
        daily_sea_ice[[index_missing[i]]][missing_values_i$days_missing] <-
          median(daily_sea_ice[[ index_missing[i] ]][missing_values_i$days_missing[1] - 1],
                 daily_sea_ice[[ index_missing[i] ]][missing_values_i$days_missing[nrow(missing_values_i)] + 1])
      } else {
        
        print(paste0("do it manually for ", years_missing[i]))
        
      }
    }
  }
  return(daily_sea_ice)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
correct_SI_mistakes <- function(daily_sea_ice, years) {
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
  wrong_measures <- df_delta %>%
    filter(abs(delta) > 10)
  
  # Get the year number
  i_years <- c()
  for (k in 1:nrow(wrong_measures)) {
    i_years <- c(i_years, which(years == wrong_measures$year[k]))
  }
  wrong_measures <- wrong_measures %>%
    mutate(year_i = i_years)
  
  daily_sea_ice_corrected <- daily_sea_ice
  index <- seq(from = 2, to = 14, by = 2)
  for (k in index) {
    daily_sea_ice_corrected[[wrong_measures$year_i[k]]][wrong_measures$day_nbr[k]] <- 
      median(x = c(daily_sea_ice_corrected[[wrong_measures$year_i[k]]][wrong_measures$day_nbr[k] - 1],
                   daily_sea_ice_corrected[[wrong_measures$year_i[k]]][wrong_measures$day_nbr[k] + 1]))
    
  }
  return(daily_sea_ice_corrected)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
correct_SI_mistakes_trial_update <- function(daily_sea_ice, years) {
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
  wrong_measures <- df_delta %>%
    filter(abs(delta) > 15) %>%
    mutate(days = 1)
  
  # How many different mistakes there are
  years_mistakes <- distinct(wrong_measures$year)
  ID_mistakes <- c()
  for (i in 1:length(years_mistakes)) {
    mistake_year_i <- wrong_measures %>%
      filter(year = years_mistakes[i])
    if (nrow(mistake_year_i) == 2) {
      ID_mistakes <- c(ID_mistakes, i, i)
    } else {
      remaining_rows <- nrow(mistake_year_i)
      while(remaining_rows > 0) {
        j <- 1
        while(mistake_year_i$day[j + 1] - mistake_year_i$day[j + 1] == 1 && j <= nrow(mistake_year_i)) {
          j <- j + 1
        }
        ID_mistakes <- c(ID_mistakes, rep(i, times = j))
        mistake_year_i <- mistake_year_i[-c(1:j),]
        remaining_rows <- nrow(mistake_year_i)
        
      }
    }
    
  }
  wrong_measures <- wrong_measures %>%
    mutate(ID_mistakes = wrong_measures)
  # Get the year number
  i_years <- c()
  for (k in 1:nrow(wrong_measures)) {
    i_years <- c(i_years, which(years == wrong_measures$year[k]))
  }
  wrong_measures <- wrong_measures %>%
    mutate(year_i = i_years)
  
  daily_sea_ice_corrected <- daily_sea_ice
  index <- seq(from = 2, to = nrow(wrong_measures), by = 2)
  for (k in index) {
    daily_sea_ice_corrected[[ wrong_measures$year_i[k] ]][wrong_measures$day_nbr[k]] <- 
      median(x = c(daily_sea_ice_corrected[[wrong_measures$year_i[k]]][wrong_measures$day_nbr[k] - 1],
                   daily_sea_ice_corrected[[wrong_measures$year_i[k]]][wrong_measures$day_nbr[k] + 1]))
    
  }
  return(daily_sea_ice_corrected)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
get_halfway_sea_ice_concentration <- function(daily_sea_ice) {
  years <- seq(from = 1991, to = 2020, by = 1)
  
  # Find best window for max sea ice extent: Which 30-day interval has the highest
  # sea ice concentration average over the ~30 years?
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
  
  # Find best window for min sea ice extent: Which 30-day interval has the lowest
  # sea ice concentration average over the ~30 years?
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
  return(halfway_sea_ice)
}





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
get_SI_retreat_advance_date <- function(daily_sea_ice, years, days_in_a_row) {
  years <- seq(from = 1991, to = 2020, by = 1)
  
  # Calculate the halfway sea ice concentration
  halfway_sea_ice <- get_halfway_sea_ice_concentration(daily_sea_ice)
  
  ### Retrieve the day number of the retreat date for each date
  # First: calculate the number of days out of days_in_a_row with ice < halfway
  retreat_days <- c()
  year <- c()
  for (i in 2:length(daily_sea_ice)) {
    year <- c(year, years[i])
    days_SI <- c()
    for (j in 1:250) {
      days_SI <- c( days_SI,
                    sum(daily_sea_ice[[i]][seq(from = j, to = j + days_in_a_row - 1, by = 1)] < halfway_sea_ice) )
      
    }
    # Second: starting from a date with many ice-free days in a row, go back in time 
    # day by day until the number of ice free days out of days_in_a_row is lower than 
    # the maximum (i.e. lower than days_in_a_row).
    j <- 250
    while(days_SI[j] == days_in_a_row) {
      j <- j - 1
    }
    retreat_days <- c(retreat_days, j + 1)
  }
  retreat <- data.frame(year = year,
                        day_retreat = retreat_days)
  
  
  ### Retrieve the day number of the retreat date for each date
  # First: calculate the number of days out of days_in_a_row with ice > halfway
  advance_days <- c()
  year <- c()
  for (i in 2:( length(daily_sea_ice) - 1)) {
    # Here we fuse days from year i to days from year i+1 because sea ice advance 
    # sometimes occur in early year i + 1
    daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                         daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
    year <- c(year, years[i])
    days_SI <- c()
    for (j in 1:length(daily_sea_ice_i)) {
      days_SI <- c( days_SI,
                    sum(daily_sea_ice_i[seq(from = j, to = j + days_in_a_row - 1, by = 1)] > halfway_sea_ice) )
      
    }
    # Second: starting from a date with many ice-free days in a row (day 250), go  
    # forward in time day by day until the number of ice-covered days out of 
    # days_in_a_row is reaches the maximum value it can take (i.e. days_in_a_row).
    j <- 1
    while(days_SI[j] < days_in_a_row) {
      j <- j + 1
    }
    advance_days <- c(advance_days, j + 249)
  }
  advance <- data.frame(year = year,
                        day_advance = advance_days)
  
  retreat_advance <- retreat %>%
    left_join(x = retreat,
              y = advance,
              by = "year") %>%
    mutate(ice_free_days = day_advance - day_retreat) 
  
  return(retreat_advance)
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

get_SI_retreat_advance_date_any_threshold <- function(daily_sea_ice, years, days_in_a_row, threshold) {
  # Retrieve the day number of retreat/advance for each year
  retreat_days <- c()
  year <- c()
  for (i in 2:length(daily_sea_ice)) {
    year <- c(year, years[i])
    days_SI <- c()
    for (j in 1:250) {
      days_SI <- c( days_SI,
                    sum(daily_sea_ice[[i]][seq(from = j, to = j + days_in_a_row - 1, by = 1)] < threshold) )
      
    }
    j <- 250
    while(days_SI[j] == days_in_a_row) {
      j <- j - 1
    }
    retreat_days <- c(retreat_days, j+1)
  }
  retreat <- data.frame(year = year,
                        day_retreat = retreat_days)
  
  advance_days <- c()
  year <- c()
  for (i in 2:( length(daily_sea_ice) - 1)) {
    daily_sea_ice_i <- c(daily_sea_ice[[i]][seq(from = 250, to = length(daily_sea_ice[[i]]), by = 1)],
                         daily_sea_ice[[i+1]][seq(from = 1, to = 135, by = 1)])
    year <- c(year, years[i])
    days_SI <- c()
    for (j in 1:length(daily_sea_ice_i)) {
      days_SI <- c( days_SI,
                    sum(daily_sea_ice_i[seq(from = j, to = j + days_in_a_row - 1, by = 1)] > threshold) )
      
    }
    j <- 1
    while(days_SI[j] < days_in_a_row) {
      j <- j + 1
    }
    advance_days <- c(advance_days, j + 249)
  }
  advance <- data.frame(year = year,
                        day_advance = advance_days)
  
  retreat_advance <- retreat %>%
    left_join(x = retreat,
              y = advance,
              by = "year") %>%
    mutate(ice_free_days = day_advance - day_retreat) 
  
  return(retreat_advance)
}


