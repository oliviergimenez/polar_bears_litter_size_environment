#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#              Analysis of correlations between AO and temperature             #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(lubridate)
library(cowplot)
library(corrplot)
library(climwin)
Sys.setenv(LANG = "en")


# A. Linear correlations =======================================================

# ~ 1. Load and process the data -----------------------------------------------

# Cliamtic oscillation
AO <- read_csv("06_processed_data/AO_data/AO_data_winter_spring_1990-2019.csv")
NAO <- read_csv("06_processed_data/NAO_data/NAO_data_winter_spring_1990-2019.csv")

# Temperature
T_data <- read_csv("06_processed_data/temperature_data/temperature_processed_monthly") %>%
  mutate(season = ifelse(month %in% c("01", "02", "03"), "winter",
                         ifelse(month %in% c("04", "05", "06"), "spring", NA))) %>%
  filter(!is.na(season)) %>%
  group_by(year, season) %>%
  summarize(T_season = mean(T_monthly)) %>%
  pivot_wider(names_from = season, values_from = T_season) %>%
  rename(T_spring = spring,
         T_winter = winter)

# Sea ice
sea_ice <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
sea_ice <- data.frame(sea_ice,
                      ice_free_days_previous = c(NA, sea_ice$ice_free_days[-nrow(sea_ice)]),
                      ice_free_days_2y_ago = c(NA, NA, sea_ice$ice_free_days[-c(nrow(sea_ice), 
                                                                                nrow(sea_ice) - 1)]),
                      day_retreat_previous = c(NA, sea_ice$day_retreat[-nrow(sea_ice)]),
                      day_retreat_2y_ago = c(NA, NA, sea_ice$day_retreat[-c(nrow(sea_ice), 
                                                                            nrow(sea_ice) - 1)]),
                      day_advance_previous = c(NA, sea_ice$day_advance[-nrow(sea_ice)]),
                      day_advance_2y_ago = c(NA, NA, sea_ice$day_advance[-c(nrow(sea_ice), 
                                                                            nrow(sea_ice) - 1)])) %>%
  dplyr::select(year, day_retreat, day_retreat_previous, day_retreat_2y_ago,
                day_advance, day_advance_previous, day_advance_2y_ago,
                ice_free_days, ice_free_days_previous, ice_free_days_2y_ago)



AO_NAO_T_SI <- AO %>%
  left_join(x = AO,
            y = NAO,
            by = "year") %>%
  dplyr::select(year, winter_AO, prior_winter_AO, 
                winter_AO_2y_ago = two_year_prior_winter_AO,
                spring_AO, prior_spring_AO, 
                spring_AO_2y_ago = two_year_prior_spring_AO,
                winter_NAO, prior_winter_NAO, 
                winter_NAO_2y_ago = two_year_prior_winter_NAO,
                spring_NAO, prior_spring_NAO, 
                spring_NAO_2y_ago = two_year_prior_spring_NAO) %>%
  left_join(x = .,
            y = T_data,
            by = "year") %>%
  left_join(x = .,
            sea_ice,
            by = "year") %>%
  filter(year >= 1992)

AO_NAO_T_SI_not_2y <- AO_NAO_T_SI %>%
  dplyr::select(-ice_free_days_2y_ago, -day_retreat_2y_ago, -day_advance_2y_ago,
                -winter_AO_2y_ago, -spring_AO_2y_ago, 
                -winter_NAO_2y_ago, -spring_NAO_2y_ago)

colnames(AO_NAO_T_SI_not_2y)


# ~ 2. Visualize correlations --------------------------------------------------

# Correlations between climate oscillation indexes
res <- cor(AO_NAO_T_SI_not_2y[,c(2:9)], method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")


# Correlations among climate oscillation indexes + temperature
res <- cor(AO_NAO_T_SI_not_2y[,c(2:11)], method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")

library("Hmisc")
res2 <- rcorr(as.matrix(AO_NAO_T_SI_not_2y[,c(2:11)]))[[3]]
res2 <- res <- round(res2, 3)


# Correlations among sea ice metrics
res <- cor(AO_NAO_T_SI_not_2y[,c(12:17)], method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")


# Correlation between sea ice metrics + climate oscillations + temperature

res <- cor(AO_NAO_T_SI_not_2y[,-1], method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")


ggplot(AO_NAO_T_SI, aes(x = winter_AO, y = day_retreat)) +
  geom_point() +
  theme_bw() +
  labs(x = "Winter AO year t",
       y = "day of sea ice retreat year t")

x <- lm(day_retreat ~ winter_AO, data = AO_NAO_T_SI)
summary(x)
y <- lm(day_retreat ~ winter_NAO, data = AO_NAO_T_SI)
summary(y)
ggsave("10_meetings/2021-04-26 Meeting with Sarah/winter AO VS day retreat.png",
       height = 4, width = 6)

ggplot(AO_NAO_T_SI, aes(x = winter_NAO, y = day_retreat)) +
  geom_point() +
  theme_bw() +
  labs(x = "Winter NAO year t",
       y = "day of sea ice retreat year t")
ggsave("10_meetings/2021-04-26 Meeting with Sarah/winter NAO VS day retreat.png",
       height = 4, width = 6)


ggplot(AO_NAO_T_SI, aes(x = spring_AO, y = day_advance)) +
  geom_point()

ggplot(AO_NAO_T_SI, aes(x = spring_NAO, y = day_advance)) +
  geom_point()



# B. Principal Compoenent Analysis (PCA) =======================================

library(FactoMineR)
library(ade4)

AO_NAO_T_SI_not_2y <- AO_NAO_T_SI %>%
  dplyr::select(-ice_free_days_2y_ago, -day_retreat_2y_ago, -day_advance_2y_ago,
                -winter_AO_2y_ago, -spring_AO_2y_ago, 
                -winter_NAO_2y_ago, -spring_NAO_2y_ago)

# Run PCA
res.pca <- PCA(AO_NAO_T_SI_not_2y, graph = FALSE)

data.no.na <- AO_NAO_T_SI_not_2y[-is.na(AO_NAO_T_SI_not_2y), ]

res.pca <- dudi.pca(df = na.omit(AO_NAO_T_SI_not_2y), 
                    scale = T, scannf = FALSE, nf = 2)

# Plot
PCA(AO_NAO_T_SI_not_2y)

# Save plot
png("07_results/01_interim_results/correlations/PCA_no_2y_prior.png")
plot(res.pca, choix = "var", autoLab = "yes")
# Close device
dev.off()


# 
# Inertia captured by each axis
x <- inertia.dudi(res.pca)
inertia <- data.frame(axis = seq(1, 12, 1),
                      inertia = x[["tot.inertia"]][["inertia"]],
                      inertia_percent = 100*x[["tot.inertia"]][["inertia"]]/sum(x[["tot.inertia"]][["inertia"]]),
                      cum =  x[["tot.inertia"]][["cum"]],
                      cum_percent =  x[["tot.inertia"]][["cum(%)"]]) %>%
  arrange(-inertia_percent)

ggplot(inertia, aes(x = axis, y = inertia_percent)) + 
  geom_col(fill = "grey20", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "",
       y = "Inertia (%)")

ggsave("07_results/01_interim_results/correlations/inertia.png",
       width = 6, height = 4)


# Percentage of varians captured by each axis
eig.val <- as.data.frame(res.pca$eig) %>%
  mutate(variable = colnames(AO_NAO_T_SI_not_2y)) %>%
  rename(percentage_variance = 'percentage of variance',
         cumulative_percentage_variance = 'cumulative percentage of variance') %>% 
  arrange(-percentage_variance)

ggplot(eig.val, aes(x = variable, y = percentage_variance)) + 
  geom_col(fill = "grey20", color = "black") +
  scale_x_discrete(limits = eig.val$variable) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Variable",
       y = "Percentage of the variance")

ggsave("07_results/01_interim_results/correlations/percentage_variance.png",
       width = 6, height = 4)


# C. Climwin ===================================================================

library(climwin)

# Temperature
temperature.raw <- read_delim("04_raw_data/temperature/svalbard_temperature_daily_montly_1898-2018.csv", 
                              ";", escape_double = FALSE, col_types = cols(X6 = col_skip(), 
                                                                           X7 = col_skip()), trim_ws = TRUE)
T.processed <- temperature.raw %>%
  mutate(Month = ifelse(Month < 10, paste0(0, Month), 
                        Month),
         Day = ifelse(Day < 10, paste0(0, Day), 
                      Day),
         date = paste0(Day, "/", Month, "/", Year),
         T_daily = as.numeric(str_replace(T_daily, ",", "."))) %>%
  dplyr::select(month = Month,
                day = Day,
                year = Year,
                date,
                T_daily) %>%
  filter(year >= 1990)

# Arctic Oscillation (AO)
AO_raw <- read_table2("04_raw_data/arctic_oscillation/monthly_actic_oscillation_index.txt", 
                      col_types = cols(`1950` = col_skip()))

AO_processed <- AO_raw %>%
  filter(year > 1985,
         year < 2021) %>%
  mutate(month = ifelse(month < 10, paste0(0, month),
                        month),
         day = 01,
         date = paste0(day, "/", month, "/", year),
         positive = ifelse(index > 0, 1, 0)) %>%
  select(year, index, positive, date, month, day)


output <- slidingwin(xvar = list(index = AO_processed$index),
                     cdate = AO_processed$date, 
                     bdate = T.processed$date, 
                     baseline = glm(T_daily ~ 1, data = T.processed),
                     range = c(24, 0), 
                     type = "relative", stat = "mean", 
                     func = c("lin"), 
                     cmissing = FALSE, 
                     cinterval = "month") 


output[[1]]$BestModel
output[[1]]$BestModelData
output[[1]]$Dataset
output[["combos"]]

?glm

#533747, #5f506b, #6a6b83, #076949f, #86bbbd