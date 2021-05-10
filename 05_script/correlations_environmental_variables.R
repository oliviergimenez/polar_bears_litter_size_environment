#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#              Analysis of correlations between AO and temperature             #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(lubridate)
library(cowplot)
library(corrplot)
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

AO_NAO <- AO %>%
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
                spring_NAO_2y_ago = two_year_prior_spring_NAO) 






# Correlations between climate oscillation indexes
res <- cor(AO_NAO, method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")

# Correlations among climate oscillation indexes + temperature
AO_NAO_T <- AO_NAO %>%
  left_join(x = AO_NAO,
            y = T_data,
            by = "year") 

AO_NAO_T <- AO_NAO_T %>%
  dplyr::select(-year)

res <- cor(AO_NAO_T, method = "pearson", use = "complete.obs")
res <- round(res, 2)
corrplot(res, type = "lower")

library("Hmisc")
res2 <- rcorr(as.matrix(AO_NAO_T))[[3]]
res2 <- res <- round(res2, 3)




# Correlations among sea ice metrics
res <- cor(sea_ice[,-1], method = "pearson", use = "complete.obs")
res <- round(res, 2)
library(corrplot)
corrplot(res, type = "lower")


# Correlation between sea ice metrics and climate oscillations + temperature
AO_NAO_T_SI <- AO_NAO_T %>%
  left_join(x = AO_NAO,
            y = sea_ice,
            by = "year") %>%
  dplyr::select(-year)


res <- cor(AO_NAO_T_SI, method = "pearson", use = "complete.obs")
res <- round(res, 2)
library(corrplot)
corrplot(res, type = "lower")


ggplot(AO_NAO_T_SI, aes(x = winter_AO, y = day_retreat)) +
  geom_point() +
  theme_bw() +
  labs(x = "Winter AO year t",
       y = "day of sea ice retreat year t")
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

# Plot
PCA(AO_NAO_T_SI_not_2y)

# Save plot
png("07_results/01_interim_results/correlations/PCA_no_2y_prior.png")
plot(res.pca, choix = "var", autoLab = "yes")
# Close device
dev.off()




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


