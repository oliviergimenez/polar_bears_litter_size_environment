#==============================================================================#
#                                                                              #
#                           Denning data processing                            #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
Sys.setenv(LANG = "en")


# Load the data
denning_females_raw <- read_delim("04_raw_data/polar_bears/list_denning_females.csv", 
                                   ";", escape_double = FALSE, trim_ws = TRUE)
colnames(denning_females_raw)
head(denning_females_raw)

# Process rows
denning_females_processed <- denning_females_raw %>%
  rename(date = DATO,
         year = YEAR,
         ID_NR = `ID NR`,
         merking = MERKING,
         location = STED,
         lat = LAT,
         long = LON,
         age_group = `Age group`,
         cub_status = `Cub status`,
         cub_number = `Cub number`) %>%
  mutate(merking = ifelse(merking == "not", "no", merking)) %>%
  # count(cub_status, cub_number)
  filter(year >= 1992,
         cub_status != "c",    # Remove females with cubs
         cub_status != "y",    # Remove females with yearlings
         cub_status != "2y") %>% #  Remove females with 2y olds
  filter(merking %in% c("n", "r")) 

  
  
females_lost_litter <- denning_females_processed %>%
  filter(Den == 1) %>%
  mutate(ID_NR_year = paste0(ID_NR, "_", year))

  



# Creat dataset of females with cubs + females who lost their cubs

capture_data_females <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

capture_data_females_without_cubs <- capture_data_females %>%
  filter(cub_status == "n") %>%
  mutate(ID_NR_year = paste0(ID_NR, "_", year)) %>%
  filter(ID_NR_year %in% females_lost_litter$ID_NR_year) %>%
  select(-ID_NR_year)

capture_data_females_w_cubs <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")


capture_data_females_w_cubs_or_litter_loss <- rbind(capture_data_females_without_cubs,
                                                    capture_data_females_w_cubs)

write_csv(capture_data_females_w_cubs_or_litter_loss, "06_processed_data/CR_data/CR_f_w_cubs_or_litter_loss.csv")

