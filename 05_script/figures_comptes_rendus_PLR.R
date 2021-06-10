#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#                     Figures pour le compte-rendu PLR2                        #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(lubridate)
library(raster)
library(grDevices)
library(mlogit)
library(nimble)
library(viridis)
Sys.setenv(LANG = "en")


# + 1. Arctic oscillation ------------------------------------------------------

AO_data <- read.csv("06_processed_data/AO_data/AO_data_1990-2019_processed.csv") %>%
  mutate(date = as.Date(date))



ggplot(data = AO_data, aes(x = date, y = index, fill = as.factor(positive))) +
  geom_col() +
  theme_bw() +
  theme(legend.position = "",
        axis.title.x = element_text(size = 7, vjust = 5),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, vjust = -57),
        plot.title.position = "plot") +
  scale_fill_manual(limits = c("1", "0"),
                    values = c("#132b43", "#56b1f7")) +
  labs(x = "Year",
       y = "Arctic oscillation index",
       title = "Figure 2. Monthly Arctic oscillation index") 

ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/AO_1990-2020.svg",
       width = 8, height = 4, unit = "cm")



# + 2. Temperature -------------------------------------------------------------

T_data <- read.csv("06_processed_data/temperature_data/temperature_processed_monthly") %>%
  mutate(date = as.Date(date))


ggplot(data = T_data, aes(x = date, y = T_monthly)) +
  geom_line(size = 0.3) +
  # geom_col() +
  theme_bw() +
  theme(legend.position = "",
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  labs(x = "Month/Year",
       y = "Monthly mean temperature") 

ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/T_1990-2020.svg",
       width = 8, height = 4, unit = "cm")

glimpse(T_data)



# Monthly
# ggplot(data = T_data, aes(x = month, y = T_monthly, alpha = as.factor(year))) +
#   geom_line(size = 0.3) +
#   theme_bw() +
#   theme(legend.position = "",
#         axis.title.x = element_text(size = 7),
#         axis.title.y = element_text(size = 7),
#         axis.text = element_text(size = 7)) +
#   scale_alpha_discrete(range = c(0.05, 1)) +
#   scale_x_continuous(breaks = c(1, 3, 6, 9, 12),
#                      labels = c("Jan", "Mar", "Jun", "Sep", "Dec")) +
#   labs(x = "Month",
#        y = "Monthly mean temperature")
# 
# ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/T_1990-2020_superposés.svg",
#        width = 8, height = 4, unit = "cm")


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
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, vjust = -57),
        plot.title.position = "plot") +
  scale_color_manual(values = my_colors,
                     breaks = c("1990", "2000", "2010", "2018")) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("Jan", "Apr", "Aug", "Oct")) +
  labs(x = "",
       y = "Monthly mean temperature",
       color = "Year",
       title = "Figure 3. Monthly mean temperature by year")

ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/T_1990-2020_superposés2.svg",
       width = 8, height = 4, unit = "cm")


# + 3. Sea ice concentration per year ------------------------------------------

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
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, vjust = -57),
        plot.title.position = "plot") +
  scale_color_manual(values = my_colors,
                     breaks = c("1990", "2000", "2010", "2018")) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("Janvier", "Avril", "Août", "Octobre")) +
  labs(x = "",
       y = "Concentration de banquise",
       color = "Année")
# title = "Figure 3. Monthly mean temperature by year")

ggsave("D:/ETUDES/4A 2020-2021/Mon PLR en 150 secondes/sea_ice_monthly.png",
       width = 8, height = 4, unit = "cm")


length(seq(from = 0.1, to = 0.96, by = 0.03))


scale_alpha_discrete()

# Daily
T_data_daily <- read_delim("04_raw_data/temperature/svalbard_temperature_daily_montly_1898-2018.csv", 
                              ";", escape_double = FALSE, col_types = cols(X6 = col_skip(), 
                                                                           X7 = col_skip()), trim_ws = TRUE) %>%
  mutate(Month = ifelse(Month < 10, paste0(0, Month), 
                        Month),
         Day = ifelse(Day < 10, paste0(0, Day), 
                      Day),
         date = as.Date(paste0(Year, "-", Month, "-", Day)),
         T_daily = as.numeric(str_replace(T_daily, ",", "."))) %>%
  dplyr::select(month = Month,
                day = Day,
                year = Year,
                date,
                T_daily,
                code = Code) %>%
  filter(date >= "1990-01-01") %>%
  mutate(day_nbr = (yday(date)))

ggplot(data = T_data_daily, aes(x = day_nbr, y = T_daily, alpha = as.factor(year))) +
  geom_line(size = 0.3) +
  theme_bw() +
  theme(legend.position = "",
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  labs(x = "Month",
       y = "Monthly mean temperature")

ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/T_1990-2020_superposés_daily.png",
       width = 8, height = 4, unit = "cm")




# + 4. Sea ice raster ----------------------------------------------------------


day.nbr <- lubridate::yday(as.Date("2020-09-15"))
day.nbr <- lubridate::yday(as.Date("2020-03-15"))



load("04_raw_data/sea_ice/hamburg/1992/raster_stack/sea_ice_raster_stack_1992.RData")
Sea_ice_stack_1992 <- sea_ice_raster_stack
rm(sea_ice_raster_stack)

# sea_ice_layer_mar <- Sea_ice_stack_1992[[75]] # March

sea_ice_layer_sep <- Sea_ice_stack_1992[[165]] # September


df_sea_ice_sep <- as.data.frame(sea_ice_layer_sep, 
                                xy = TRUE)

rects <- data.frame(xstart = 600, 
                            xend = 1600, 
                            ystart = -1000,
                            yend = 0)

ggplot() +
  geom_raster(data = df_sea_ice_sep, aes(x = x, y = y, fill = X1992.01.01.165)) +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend,
                              ymin = ystart, ymax = yend),
            color = "red", alpha = 0) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5, vjust = -90),
        plot.title.position = "plot") +
  
  labs(x = "",
       y = "",
       fill = "Sea ice \ncover",
       title = "Figure 1. Sea ice cover in the study area (15-09-1992)") +
  scale_fill_continuous(low = "#132b43", high = "white", 
                        guide = "colorbar", na.value = "#cccccc",
                        breaks = c(0, 50, 100),
                        labels = c("0%", "50%", "100%"))


ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/sea_ice_sep.svg",
       width = 9, height = 6, unit = "cm")


# + 5. Litter size years -------------------------------------------------------

capture_data_females_w_cubs <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

capture_data_females_w_cubs %>%
    mutate(CUB_NUMBER = as.numeric(CUB_NUMBER)) %>%
    ggplot(aes(x = as.numeric(YEAR), y = CUB_NUMBER)) +
    geom_jitter(width = 0.20,
                height = 0.10,
                size = 0.5,
                color = "grey20") +
    scale_y_continuous(breaks = c(1, 2, 3),
                       labels = c(1, 2, 3)) +
    scale_x_discrete(limits = c(2000, 2010),
                     labels = c("2000", "2010")) + 
    # geom_smooth(method = lm) +
    theme_bw() +
    theme(axis.text = element_text(size = 7),
          axis.title.x = element_text(size = 7, vjust = 5),
          axis.title.y = element_text(size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, vjust = -58),
          plot.title.position = "plot") +
    labs(x = "Year",
         y = "Litter size",       
         title = "Figure 4. Captured females' litter size")




ggsave("D:/ETUDES/4A 2020-2021/Compte rendu PLR/PLR2/litter_size_year.svg",
       width = 8, height = 4, unit = "cm")



# + 6. Model output ice-free-days ----------------------------------------------

load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")

N <- dim(fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]])[1]
res <- rbind(fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]][seq(1, N, by = 4), 
                                                                       c(1:5)],
             fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain2"]][seq(1, N, by = 4), 
                                                                       c(1:5)])
b1cub <- res[, c(1:4)]
b2_3cub <- res[, c(5, 2:4)] 

range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(prior_spring_AO_s) +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 4] * mean(day_s))
    
  }
}
end <- Sys.time() 
end - start

p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.mean.AO <- data.frame(var = grid,
                                  mean_p_0_cub = apply(p0cub, 2, mean),
                                  mean_p_1_cub = apply(p1cub, 2, mean),
                                  mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                  ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                  ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                  ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                  ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                  ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                  ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         AO = "mean")


ggplot(data = df.for.plot.mean.AO,
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("pas d'oursons", "1 ourson", "2-3 oursons")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Moyenne", "IC")) +
  theme_bw() +
  labs(x = "Jours sans banquise l'année passée",
       y = "Probabilité",
       color = "",
       linetype = "")

ggsave(filename = "D:/ETUDES/4A 2020-2021/Mon PLR en 150 secondes/graph_ice_free_days_0cubs_included.png",
       width = 5, height = 3)


# + 7. Model output only ice-free-days starting 2000 --------------------------------

CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
sea_ice_data <- data.frame(sea_ice_data,
                           ice_free_days_previous = c(NA, sea_ice_data$ice_free_days[-nrow(sea_ice_data)]),
                           ice_free_days_2y_prior = c(NA, NA, sea_ice_data$ice_free_days[-c(nrow(sea_ice_data),
                                                                                            nrow(sea_ice_data) - 1)]))


data_model <- CR_data %>%
  left_join(x = CR_data,
            y = sea_ice_data,
            by = "year")


# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Renumérotation des années
{year <- data_model$year
  year <- factor(year) # laisse tomber les modalites inutiles 
  year2 <- NULL
  for (i in 1:length(year)){
    year2 <- c(year2, which(year[i] == levels(year)))
  }
  year <- factor(year2)
  nbyear <- length(levels(year))
  year
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/1_sea_ice/1.1_Nimble_ice_free_days.R")
source("05_script/models/functions_for_models_Nimble.R")

{model_code <- "1.1.2_D"
  effect <- "common"
  
  # Predictor
  var <- data_model$ice_free_days_previous
  var_scaled <- (var - mean(var))/sd(var) 
  var_short_name <- "ice_free_days_previous_s"
  var_full_name <- "Ice-free days in previous year"
  
  # Are females without cubs taken into account ?
  mode <- ""       # Yes
  # mode <- "_bis"  # No
  
  my.constants <- list(N = length(y), # nb of females captured
                       J = length(levels(y)),
                       year = as.numeric(year),
                       nbyear = nbyear,
                       as.numeric(var_scaled)) 
  names(my.constants)[5] <- var_short_name
}



load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))


res <- rbind(get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain1,
             get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain2)

# Create grid of x values
range <- range(var_scaled)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(var) + mean(var)


q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)

b1cub <- res[, c(1, 2)]
b2_3cub <- res[, c(3, 2)] 

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid_scaled[i])
  }
}

# Backtransform
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot <- data.frame(var = grid,
                          mean_p_0_cub = apply(p0cub, 2, mean),
                          mean_p_1_cub = apply(p1cub, 2, mean),
                          mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                          ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                          ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", "credible_interval"))

color_labels <- c("no cubs", "1 cub", "2-3 cubs")

ggplot(data = df.for.plot, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = color_labels) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = var_full_name,
       y = "Probability", 
       color = "",
       linetype = "") +
  scale_x_continuous(limits = c(50, 325)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1))

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/special_models/model_", 
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)

# + 8. Model output only ice-free-days starting 2000 ---------------------------

CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
sea_ice_data <- data.frame(sea_ice_data,
                           ice_free_days_previous = c(NA, sea_ice_data$ice_free_days[-nrow(sea_ice_data)]),
                           ice_free_days_2y_prior = c(NA, NA, sea_ice_data$ice_free_days[-c(nrow(sea_ice_data),
                                                                                            nrow(sea_ice_data) - 1)]))


data_model <- CR_data %>%
  left_join(x = CR_data,
            y = sea_ice_data,
            by = "year") %>%
  filter(year > 1999)


# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Renumérotation des années
{year <- data_model$year
  year <- factor(year) # laisse tomber les modalites inutiles 
  year2 <- NULL
  for (i in 1:length(year)){
    year2 <- c(year2, which(year[i] == levels(year)))
  }
  year <- factor(year2)
  nbyear <- length(levels(year))
  year
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/1_sea_ice/1.1_Nimble_ice_free_days.R")
source("05_script/models/functions_for_models_Nimble.R")

{model_code <- "1.1.2_D"
  effect <- "common"
  
  # Predictor
  var <- data_model$ice_free_days_previous
  var_scaled <- (var - mean(var))/sd(var) 
  var_short_name <- "ice_free_days_previous_s"
  var_full_name <- "Ice-free days in previous year"
  
  # Are females without cubs taken into account ?
  mode <- ""       # Yes
  # mode <- "_bis"  # No
  
  my.constants <- list(N = length(y), # nb of females captured
                       J = length(levels(y)),
                       year = as.numeric(year),
                       nbyear = nbyear,
                       as.numeric(var_scaled)) 
  names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect, mode)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         a1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect)),     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits,          
                          monitors = params,   # parameters to monitor
                          thin = 10,
                          niter = 20000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE))
end <- Sys.time()
end - start

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/special_models/model_", 
                   model_code, "_effect_", effect, "_starting_2000.RData"))


load(file = paste0("07_results/01_interim_results/model_outputs/special_models/model_", 
                   model_code, "_effect_", effect, "_starting_2000.RData"))


res <- rbind(get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain1,
             get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain2)

# Create grid of x values
range <- range(var_scaled)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(var) + mean(var)


q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)

b1cub <- res[, c(1, 2)]
b2_3cub <- res[, c(3, 2)] 

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid_scaled[i])
  }
}

# Backtransform
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot <- data.frame(var = grid,
                          mean_p_0_cub = apply(p0cub, 2, mean),
                          mean_p_1_cub = apply(p1cub, 2, mean),
                          mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                          ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                          ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", "credible_interval"))

color_labels <- c("no cubs", "1 cub", "2-3 cubs")

ggplot(data = df.for.plot, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = color_labels) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = var_full_name,
       y = "Probability", 
       color = "",
       linetype = "") +
  scale_x_continuous(limits = c(50, 325)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1))

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/special_models/model_", 
                         model_code,  "_effect_", effect, "_starting_from_2000.png"),
       width = 6, height = 3)



