#==============================================================================#
#                                                                              #
#                      Models with environmental covariates                    #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(nimble)
library(viridis)
library(cowplot)
Sys.setenv(LANG = "en")


# A. READY VARIABLES ===========================================================

# ~~~ a. Load the data ---------------------------------------------------------
CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
sea_ice_data <- data.frame(sea_ice_data,
                           ice_free_days_previous = c(NA, sea_ice_data$ice_free_days[-nrow(sea_ice_data)]),
                           ice_free_days_2y_prior = c(NA, NA, sea_ice_data$ice_free_days[-c(nrow(sea_ice_data),
                                                                                            nrow(sea_ice_data) - 1)]))
AO_data <- read_csv("06_processed_data/AO_data/AO_data_winter_spring_1990-2019.csv")
NAO_data <- read_csv("06_processed_data/NAO_data/NAO_data_winter_spring_1990-2019.csv")


CR_data %>%
  left_join(x = CR_data,
            y = sea_ice_data,
            by = "year") %>% 
  left_join(.,
            y = AO_data,
            by = "year") %>% 
  left_join(.,
            y = NAO_data,
            by = "year") %>%
  filter(year >= 2000) -> data_model

rm(CR_data, sea_ice_data, AO_data, NAO_data)



# ~~~ b. Format data for Nimble ------------------------------------------------

y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- scale(day)

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- scale(ice_free_days_previous)

# Winter AO
winter_AO <- data_model$winter_AO
winter_AO_s <- scale(winter_AO)

# Spring AO t-1 
prior_spring_AO <- data_model$prior_spring_AO
prior_spring_AO_s <- scale(prior_spring_AO)


# Random effects ++++++++++++++++++++++++++++++++++++

# years
{year <- data_model$year
  year <- factor(year) # laisse tomber les modalites inutiles 
  year2 <- NULL
  for (i in 1:length(year)) {
    year2 <- c(year2, which(year[i] == levels(year)))
  }
  year <- factor(year2)
  rm(year2)
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/several_predictors/1_Nimble_environmental_covariates.R")
source("05_script/models/functions_for_models_Nimble.R")


# B. RUN THE MODELS ============================================================

# 1. Ice-free days t-1 + spring AO t-1: common slopes --------------------------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s  +  prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.2.2.1_common <- nimbleMCMC(code = model_1.1.2.D_2.2.2.1_common,     # model code  
                                               data = dat,                                   
                                               constants = my.constants,        
                                               inits = inits,          
                                               monitors = params,   # parameters to monitor
                                               thin = 10,
                                               niter = 80000,                  # nb iterations
                                               nburnin = 10000,              # length of the burn-in
                                               nchains = 2,
                                               summary = TRUE,
                                               WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.2.2.1_common$WAIC
# 775.5775

save(fit_1.1.2.D_2.2.2.1_common, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")



# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.2.2.1_common)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2.D_2.2.2.1_common.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")

N <- dim(fit_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]])[1]
res <- rbind(fit_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]][seq(1, N, by = 4), 
                                                                       c(1:5)],
             fit_1.1.2.D_2.2.2.1_common[["samples"]][["chain2"]][seq(1, N, by = 4), 
                                                                       c(1:5)])
b1cub <- res[, c(1:3)]
b2_3cub <- res[, c(4, 2:3)] 

# ~~~~~~ Sea ice ---------------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# Mean spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(prior_spring_AO_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(prior_spring_AO_s))
    
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


# ggplot(data = df.for.plot.mean.AO, 
#        aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE,                       
#                       labels = c("pas d'oursons", "1 ourson", "2-3 oursons")) +
#   scale_linetype_manual(limits = c("mean", "credible_interval"),
#                         values = c("solid", "dotted"),
#                         labels = c("Moyenne", "IC")) +
#   theme_bw() +
#   labs(x = "Jours sans banquise l'année passée",
#        y = "Probabilité", 
#        color = "",
#        linetype = "") 
# 
# ggsave(filename = "D:/ETUDES/4A 2020-2021/Mon PLR en 150 secondes/graph_ice_free_days_0cubs_included.png",
#        width = 5, height = 3)

# Low spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fst_decile <- quantile(x = prior_spring_AO_s, # = -0.4841
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * fst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * fst_decile)
    
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

df.for.plot.low_AO <- data.frame(var = grid,
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
         AO = "low")



# High spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lst_decile <- quantile(x = prior_spring_AO_s, # = 0.5719
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * lst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * lst_decile)
    
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


df.for.plot.high_AO <- data.frame(var = grid,
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
         AO = "high")


df.for.plot.spring_AO <- rbind(df.for.plot.low_AO, df.for.plot.mean.AO, df.for.plot.high_AO) %>%
  mutate(AO = factor(AO, levels = c("low", "mean", "high")))

write_csv(df.for.plot.spring_AO, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_ice_free_days_facet_spring_AO.csv")
df.for.plot.spring_AO <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_ice_free_days_facet_spring_AO.csv") %>%
  mutate(AO = factor(AO, levels = c("low", "mean", "high")))


labels <- c("low AO index", "mean AO index", "high AO index")
names(labels) <- c("low", "mean", "high")

ggplot(data = df.for.plot.spring_AO, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap( ~ AO,
            labeller = labeller(AO = labels))

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_ice_free_days_facet_spring_AO.png",
       width = 10, height = 4)







# ~~~~~~ Spring AO -------------------------------------------------------------
range <- range(prior_spring_AO_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(prior_spring_AO) + mean(prior_spring_AO)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# Mean number of ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * mean(ice_free_days_previous_s) +
                         b1cub[j, 3] * grid_scaled[i])
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * mean(ice_free_days_previous_s) +
                           b2_3cub[j, 3] * grid_scaled[i])
    
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


df.for.plot.mean_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "mean number of ice-free days")


# Few ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fst_decile <- quantile(x = ice_free_days_previous_s, # = 107 days
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * fst_decile +
                         b1cub[j, 3] * grid_scaled[i])
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * fst_decile +
                           b2_3cub[j, 3] * grid_scaled[i])
    
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


df.for.plot.few_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "few ice-free days")


# Many ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lst_decile <- quantile(x = ice_free_days_previous_s, # = 266 days
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * lst_decile +
                         b1cub[j, 3] * grid_scaled[i] +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * lst_decile +
                           b2_3cub[j, 3] * grid_scaled[i] +
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


df.for.plot.many_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "many ice-free days")


df.for.plot.ice_free_days <- rbind(df.for.plot.few_ice_free_days, 
                                   df.for.plot.mean_ice_free_days, 
                                   df.for.plot.many_ice_free_days) %>%
  mutate(ice_free_days = factor(ice_free_days, levels = c("few ice-free days", 
                                                          "mean number of ice-free days", 
                                                          "many ice-free days")))

write_csv(df.for.plot.ice_free_days, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_spring_AO_facet_ice_free_days.csv")
df.for.plot.ice_free_days <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_spring_AO_facet_ice_free_days.csv") %>%
  mutate(ice_free_days = factor(ice_free_days, levels = c("few ice-free days", 
                                                          "mean number of ice-free days", 
                                                          "many ice-free days")))


ggplot(data = df.for.plot.ice_free_days, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Previous spring AO index",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap(~ice_free_days)

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_spring_AO_facet_ice_free_days.png",
       width = 10, height = 4)





# 2. Ice-free days t-1 (com) + spring AO t-1 (com) + date of capture (dist) -------------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s  +  prior_spring_AO_s + day_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 2, 8)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.2.2.1_common_4_distinct <- nimbleMCMC(code = model_1.1.2.D_2.2.2.1_common_4_distinct,     # model code  
                                               data = dat,                                   
                                               constants = my.constants,        
                                               inits = inits,          
                                               monitors = params,   # parameters to monitor
                                               thin = 10,
                                               niter = 80000,                  # nb iterations
                                               nburnin = 10000,              # length of the burn-in
                                               nchains = 2,
                                               summary = TRUE,
                                               WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.2.2.1_common_4_distinct$WAIC
# 762.1947

save(fit_1.1.2.D_2.2.2.1_common_4_distinct, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.2.2.1_common_4_distinct)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.2.2.1_common_4_distinct.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct.RData")

N <- dim(fit_1.1.2.D_2.2.2.1_common_4_distinct[["samples"]][["chain1"]])[1]
res <- rbind(fit_1.1.2.D_2.2.2.1_common_4_distinct[["samples"]][["chain1"]][seq(1, N, by = 4), 
                                                                       c(1:6)],
             fit_1.1.2.D_2.2.2.1_common_4_distinct[["samples"]][["chain2"]][seq(1, N, by = 4), 
                                                                       c(1:6)])
b1cub <- res[, c(1:4)]
b2_3cub <- res[, c(5, 2:3, 6)] 

# ~~~~~~ Sea ice (facets AO)---------------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# Mean spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


# ggplot(data = df.for.plot.mean.AO, 
#        aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE,                       
#                       labels = c("pas d'oursons", "1 ourson", "2-3 oursons")) +
#   scale_linetype_manual(limits = c("mean", "credible_interval"),
#                         values = c("solid", "dotted"),
#                         labels = c("Moyenne", "IC")) +
#   theme_bw() +
#   labs(x = "Jours sans banquise l'année passée",
#        y = "Probabilité", 
#        color = "",
#        linetype = "") 
# 
# ggsave(filename = "D:/ETUDES/4A 2020-2021/Mon PLR en 150 secondes/graph_ice_free_days_0cubs_included.png",
#        width = 5, height = 3)

# Low spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fst_decile <- quantile(x = prior_spring_AO_s, # = -0.4841
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * fst_decile +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * fst_decile +
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

df.for.plot.low_AO <- data.frame(var = grid,
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
         AO = "low")



# High spring AO index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lst_decile <- quantile(x = prior_spring_AO_s, # = 0.5719
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * lst_decile +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * lst_decile +
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


df.for.plot.high_AO <- data.frame(var = grid,
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
         AO = "high")


df.for.plot.spring_AO <- rbind(df.for.plot.low_AO, df.for.plot.mean.AO, df.for.plot.high_AO) %>%
  mutate(AO = factor(AO, levels = c("low", "mean", "high")))

write_csv(df.for.plot.spring_AO, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_spring_AO.csv")
df.for.plot.spring_AO <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_spring_AO.csv") %>%
  mutate(AO = factor(AO, levels = c("low", "mean", "high")))

labels <- c("low AO index", "mean AO index", "high AO index")
names(labels) <- c("low", "mean", "high")

ggplot(data = df.for.plot.spring_AO, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap( ~ AO,
              labeller = labeller(AO = labels))

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_spring_AO.png",
       width = 10, height = 4)



# ~~~~~~ Sea ice (facets day capture)-------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# Mean day of capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

df.for.plot.mean.capture <- data.frame(var = grid,
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
         day = "mean")


# ggplot(data = df.for.plot.mean.AO, 
#        aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE,                       
#                       labels = c("pas d'oursons", "1 ourson", "2-3 oursons")) +
#   scale_linetype_manual(limits = c("mean", "credible_interval"),
#                         values = c("solid", "dotted"),
#                         labels = c("Moyenne", "IC")) +
#   theme_bw() +
#   labs(x = "Jours sans banquise l'année passée",
#        y = "Probabilité", 
#        color = "",
#        linetype = "") 
# 
# ggsave(filename = "D:/ETUDES/4A 2020-2021/Mon PLR en 150 secondes/graph_ice_free_days_0cubs_included.png",
#        width = 5, height = 3)

# Early date of capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fst_decile <- quantile(x = day_s, # = -0.4841
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(prior_spring_AO_s) + 
                         b1cub[j, 4] * fst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(prior_spring_AO_s) + 
                           b2_3cub[j, 4] * fst_decile)
    
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

df.for.plot.early_capture <- data.frame(var = grid,
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
         day = "early")



# Late date of capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lst_decile <- quantile(x = day_s, # = 0.5719
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(prior_spring_AO_s) +
                         b1cub[j, 4] * lst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 4] * lst_decile)
    
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


df.for.plot.late_capture <- data.frame(var = grid,
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
         day = "late")


df.for.plot.day_capture <- rbind(df.for.plot.early_capture, df.for.plot.mean.capture, df.for.plot.late_capture) %>%
  mutate(day = factor(day, levels = c("early", "mean", "late")))

write_csv(df.for.plot.day_capture, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_day_capture.csv")
df.for.plot.day_capture <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_day_capture.csv") %>%
  mutate(day = factor(day, levels = c("early", "mean", "late")))

labels <- c("early day of capture", "mean day of capture", "late day of capture")
names(labels) <- c("early", "mean", "late")

ggplot(data = df.for.plot.day_capture, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap( ~ day,
              labeller = labeller(day = labels))

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_ice_free_days_facet_day_capture.png",
       width = 10, height = 4)





# ~~~~~~ Spring AO -------------------------------------------------------------
range <- range(prior_spring_AO_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(prior_spring_AO) + mean(prior_spring_AO)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# Mean number of ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * mean(ice_free_days_previous_s) +
                         b1cub[j, 3] * grid_scaled[i] +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * mean(ice_free_days_previous_s) +
                           b2_3cub[j, 3] * grid_scaled[i] +
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


df.for.plot.mean_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "mean number of ice-free days")


# Few ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fst_decile <- quantile(x = ice_free_days_previous_s, # = 107 days
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * fst_decile +
                         b1cub[j, 3] * grid_scaled[i] +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * fst_decile +
                           b2_3cub[j, 3] * grid_scaled[i] +
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


df.for.plot.few_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "few ice-free days")


# Many ice-free days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lst_decile <- quantile(x = ice_free_days_previous_s, # = 266 days
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * lst_decile +
                         b1cub[j, 3] * grid_scaled[i] +
                         b1cub[j, 4] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * lst_decile +
                           b2_3cub[j, 3] * grid_scaled[i] +
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


df.for.plot.many_ice_free_days <- data.frame(var = grid,
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
         ice_free_days = "many ice-free days")


df.for.plot.ice_free_days <- rbind(df.for.plot.few_ice_free_days, 
                                   df.for.plot.mean_ice_free_days, 
                                   df.for.plot.many_ice_free_days) %>%
  mutate(ice_free_days = factor(ice_free_days, levels = c("few ice-free days", 
                                                          "mean number of ice-free days", 
                                                          "many ice-free days")))

write_csv(df.for.plot.ice_free_days, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_spring_AO_facet_ice_free_days.csv")
df.for.plot.ice_free_days <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_spring_AO_facet_ice_free_days.csv") %>%
  mutate(ice_free_days = factor(ice_free_days, levels = c("few ice-free days", 
                                                          "mean number of ice-free days", 
                                                          "many ice-free days")))


ggplot(data = df.for.plot.ice_free_days, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Previous spring AO index",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap(~ice_free_days)

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_4_distinct_spring_AO_facet_ice_free_days.png",
       width = 10, height = 4)



# 3. Ice-free days t-1 (com) * date of capture (dist) + spring AO t-1 (com) ----

# Here the coefficient for the interaction is common

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s  * day_s + prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 2, 8)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.2.2.1_common_4_distinct_2 <- nimbleMCMC(code = model_1.1.2.D_2.2.2.1_common_4_distinct_2,     # model code  
                                                    data = dat,                                   
                                                    constants = my.constants,        
                                                    inits = inits,          
                                                    monitors = params,   # parameters to monitor
                                                    thin = 10,
                                                    niter = 100000,                  # nb iterations
                                                    nburnin = 10000,              # length of the burn-in
                                                    nchains = 2,
                                                    summary = TRUE,
                                                    WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.2.2.1_common_4_distinct_2$WAIC
# 756.0809

save(fit_1.1.2.D_2.2.2.1_common_4_distinct_2, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct_2.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct_2.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.2.2.1_common_4_distinct_2)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.2.2.1_common_4_distinct_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)





# 4. Ice-free days t-1 (com) * date of capture (dist) + spring AO t-1 (com) ----

# Here the coefficients for the interaction are distinct


my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s * day_s +  prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 2, 8, 10)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.2.2.1_common_4_distinct_3 <- nimbleMCMC(code = model_1.1.2.D_2.2.2.1_common_4_distinct_3,     # model code  
                                                      data = dat,                                   
                                                      constants = my.constants,        
                                                      inits = inits,          
                                                      monitors = params,   # parameters to monitor
                                                      thin = 10,
                                                      niter = 100000,                  # nb iterations
                                                      nburnin = 10000,              # length of the burn-in
                                                      nchains = 2,
                                                      summary = TRUE,
                                                      WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.2.2.1_common_4_distinct_3$WAIC
# 756.6417

save(fit_1.1.2.D_2.2.2.1_common_4_distinct_3, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct_3.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common_4_distinct_3.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.2.2.1_common_4_distinct_3)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.2.2.1_common_4_distinct_3.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)






# 5. Ice-free days t-1 (com) + spring AO t-1 (com) + winter AO t (com) + date of capture (dist) -------------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     winter_AO_s = as.numeric(winter_AO_s),
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s  +  winter_AO_s + prior_spring_AO_s + 
                       day_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 2, 10)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct <- nimbleMCMC(code = model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct,     # model code  
                                                    data = dat,                                   
                                                    constants = my.constants,        
                                                    inits = inits,          
                                                    monitors = params,   # parameters to monitor
                                                    thin = 10,
                                                    niter = 100000,                  # nb iterations
                                                    nburnin = 10000,              # length of the burn-in
                                                    nchains = 2,
                                                    summary = TRUE,
                                                    WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct$WAIC
# 760.721 

save(fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)




# 6. Ice-free days t-1 (com) * date of capture (dist) + spring AO t-1 (com) + winter AO t (com) -------------

# /!\ The coefficient for the interaction is common

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     winter_AO_s = as.numeric(winter_AO_s),
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s * day_s + winter_AO_s + prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 11, 2, 10)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2 <- nimbleMCMC(code = model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2,     # model code  
                                                            data = dat,                                   
                                                            constants = my.constants,        
                                                            inits = inits,          
                                                            monitors = params,   # parameters to monitor
                                                            thin = 10,
                                                            niter = 100000,                  # nb iterations
                                                            nburnin = 10000,              # length of the burn-in
                                                            nchains = 2,
                                                            summary = TRUE,
                                                            WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2$WAIC
# 755.407 

save(fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2.RData")

N <- dim(fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2[["samples"]][["chain1"]])[1]
res <- rbind(fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2[["samples"]][["chain1"]][seq(1, N, by = 4), 
                                                                 c(1:8)],
             fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2[["samples"]][["chain2"]][seq(1, N, by = 4), 
                                                                 c(1:8)])
b1cub <- res[, c(1:6)]
b2_3cub <- res[, c(7, 2:4, 8, 6)] 

# ~~~~~~ Caterpillar plot ------------------------------------------------------

caterpilar <- as.data.frame(res[, c(1:8)]) %>%
  pivot_longer(cols = c("b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]"))

ggplot(data = caterpilar, 
       aes(x = name , y = value)) + 
  stat_summary(fun.data = get_mean_and_CI, 
               fun.args = list(lower = 0.05, upper = 0.95),
               geom = "pointrange", size = 0.5) +
  stat_summary(fun.data = get_mean_and_CI, 
               fun.args = list(lower = 0.20, upper = 0.80),
               geom = "pointrange", size = 1) +
  # stat_summary(fun.data = get_mean_and_CI, 
  #              fun.args = list(lower = 0.20, upper = 0.80),
  #              geom = "pointrange") +
  theme_bw() +
  labs(x = "",
       y = "") +
  coord_flip()

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_caterpillar.svg",
       width = 4, height = 3)



# ~~~~~~ Sea ice ---------------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)


# Mean everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * mean(day_s) +
                         b1cub[j, 6] * grid_scaled[i] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * mean(day_s) +
                           b2_3cub[j, 6] * grid_scaled[i] * mean(day_s))
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

df.for.plot.sea.ice <- data.frame(var = grid,
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
                       "credible_interval"))


ggplot(data = df.for.plot.sea.ice,
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("mean", "CI")) +
  theme_bw() +
  labs(x = "Ice-free days during the previous year",
       y = "Probability",
       color = "",
       linetype = "")

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_ice_free.png",
       width = 6, height = 3)

x <- apply(p0cub, 2, mean)
x <- c(apply(p0cub, 2, mean),
      apply(p1cub, 2, mean),
      apply(p2_3cub, 2, mean))

df.for.plot.sea.ice <- data.frame(var = rep(grid, times = 3),
                                  mean = c(apply(p0cub, 2, mean),
                                           apply(p1cub, 2, mean),
                                           apply(p2_3cub, 2, mean)),
                                  q2_5 = c(apply(p0cub, 2, quantile, probs = 0.025),
                                           apply(p1cub, 2, quantile, probs = 0.025),
                                           apply(p2_3cub, 2, quantile, probs = 0.025)),
                                  q97_5 = c(apply(p0cub, 2, quantile, probs = 0.975),
                                            apply(p1cub, 2, quantile, probs = 0.975),
                                            apply(p2_3cub, 2, quantile, probs = 0.975)),
                                  cub_number = c(rep(0, times = dim(p0cub)[2]),
                                                 rep(1, times = dim(p1cub)[2]),
                                                 rep(2, times = dim(p2_3cub)[2])))




ggplot(df.for.plot.sea.ice, 
       aes(x = var, y = mean, 
           color = as.factor(cub_number), group = as.factor(cub_number))) +
  geom_line(size = 0.75,
            lineend = "round") +
  geom_ribbon(aes(ymin = q2_5,
                  ymax = q97_5,
                  fill = as.factor(cub_number)), alpha = 0.3, color = NA) +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  scale_fill_viridis(discrete = TRUE,
                     labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  theme_bw() +
  labs(x = "Ice-free days during the previous year",
       y = "Probability",
       color = "",
       fill = "")

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/test.png",
       width = 6, height = 3)

################ Pour le PLR en 150 sec ###
df.for.plot.sea.ice <- data.frame(var = rep(grid, times = 2),
                                  mean = c(apply(p0cub, 2, mean),
                                           apply(p1cub, 2, mean) + apply(p2_3cub, 2, mean)),
                                  
                                  q2_5 = c(apply(p0cub, 2, quantile, probs = 0.025),
                                           apply(p1cub, 2, quantile, probs = 0.025) + apply(p2_3cub, 2, quantile, probs = 0.025)),
                                           
                                  q97_5 = c(apply(p0cub, 2, quantile, probs = 0.975),
                                            apply(p1cub, 2, quantile, probs = 0.975) + apply(p2_3cub, 2, quantile, probs = 0.975)),
                                            
                                  cub_number = c(rep("0", times = dim(p0cub)[2]),
                                                 rep("1 ou plus", times = dim(p1cub)[2])))


ggplot(df.for.plot.sea.ice, 
       aes(x = var, y = mean, 
           color = as.factor(cub_number), group = as.factor(cub_number))) +
  geom_line(size = 0.75,
            lineend = "round") +
  geom_ribbon(aes(ymin = q2_5,
                  ymax = q97_5,
                  fill = as.factor(cub_number)), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("tomato", "cornflowerblue")) +
  scale_fill_manual(values = c("tomato", "cornflowerblue")) +
  # scale_color_manual(values = c("tomato", "cornflowerblue"),
  #                    labels = "pas de portée", "portée") +
  # scale_fill_manual(values = c("tomato", "cornflowerblue"),
  #                    labels = "pas de portée", "portée")
  # scale_color_viridis(discrete = TRUE,
  #                     labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  # scale_fill_viridis(discrete = TRUE,
  #                    labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  theme_bw() +
  labs(x = "Jours sans banquise l'année n-1",
       y = "Probabilité",
       color = "oursons",
       fill = "oursons")

ggsave(filename = "D:/polar_bears_litter_size_environment/09_presentations/Mon PLR en 150 secondes/plot_mutli_2_couleurs.png",
       width = 4, height = 2.75)



# ~~~~~~ Day capture -----------------------------------------------------------
range <- range(day_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(day) + mean(day)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)


# Mean everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * mean(ice_free_days_previous_s) +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * grid_scaled[i] +
                         b1cub[j, 6] * mean(ice_free_days_previous_s) * grid_scaled[i])
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * grid_scaled[i] +
                           b2_3cub[j, 6] * mean(ice_free_days_previous_s) * grid_scaled[i])
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

df.for.plotday.capture <- data.frame(var = grid,
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
                       "credible_interval"))


ggplot(data = df.for.plotday.capture,
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Day of capture",
       y = "Probability",
       color = "",
       linetype = "")

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_date_capture.png",
       width = 6, height = 3)

# df.for.plot.mean.all <- data.frame(var = grid,
#                                    mean = c(apply(p0cub, 2, mean), 
#                                             apply(p1cub, 2, mean),  
#                                             apply(p2_3cub, 2, mean)),
#                                    q_2.5 = c(apply(p0cub, 2, quantile, probs = 0.025),
#                                              apply(p1cub, 2, quantile, probs = 0.025),
#                                              apply(p2_3cub, 2, quantile, probs = 0.025)),
#                                    q_97.5 = c(apply(p0cub, 2, quantile, probs = 0.975),
#                                               apply(p1cub, 2, quantile, probs = 0.975),
#                                               apply(p2_3cub, 2, quantile, probs = 0.975)),
#                                    cub_number = c(rep(0, times = dim(p0cub)[1]),
#                                                   rep(1, times = dim(p1cub)[1]),
#                                                   rep(2, times = dim(p2_3cub)[1]))) 
# 
# 
# 
# ggplot(data = df.for.plot.mean.all) +
#   geom_line(aes(x = var, y = mean, group = as.factor(cub_number), color = as.factor(cub_number))) +
#   # geom_ribbon(aes(x = var,
#   #                 ymin = q_2.5,
#   #                 ymax = q_97.5,
#   #                 group = as.factor(cub_number),
#   #                 fill = as.factor(cub_number))) +
#   scale_color_viridis(discrete = TRUE,
#                       labels = c("0 cubs", "1 cub", "2-3 cubs")) +
#   theme_bw() +
#   labs(x = "Jours sans banquise l'année passée",
#        y = "Probabilité",
#        color = "",
#        linetype = "")



# ~~~~~~ Winter AO -----------------------------------------------------------
range <- range(winter_AO_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(winter_AO) + mean(winter_AO)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)


# Mean everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * mean(ice_free_days_previous_s) +
                         b1cub[j, 3] * grid_scaled[i] +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * mean(day_s) +
                         b1cub[j, 6] * mean(ice_free_days_previous_s) * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * mean(ice_free_days_previous_s) +
                           b2_3cub[j, 3] * grid_scaled[i] +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * mean(day_s) +
                           b2_3cub[j, 6] * mean(ice_free_days_previous_s) * mean(day_s))
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

df.for.plot.winter.AO <- data.frame(var = grid,
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
                       "credible_interval"))


ggplot(data = df.for.plot.winter.AO,
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("mean", "CI")) +
  theme_bw() +
  labs(x = "Winter AO",
       y = "Probability",
       color = "",
       linetype = "")

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_winter_AO.png",
       width = 6, height = 3)


# ~~~~~~ Spring AO t-1-----------------------------------------------------------
range <- range(prior_spring_AO_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(prior_spring_AO) + mean(prior_spring_AO)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)


# Mean everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * mean(ice_free_days_previous_s) +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * grid_scaled[i] +
                         b1cub[j, 5] * mean(day_s) +
                         b1cub[j, 6] * mean(ice_free_days_previous_s) * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * mean(ice_free_days_previous_s) +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * grid_scaled[i] +
                           b2_3cub[j, 5] * mean(day_s) +
                           b2_3cub[j, 6] * mean(ice_free_days_previous_s) * mean(day_s))
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

df.for.plot.spring.AO <- data.frame(var = grid,
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
                       "credible_interval"))


ggplot(data = df.for.plot.spring.AO,
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0 cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Previous spring AO",
       y = "Probability",
       color = "",
       linetype = "")

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_spring_AO_t-1.png",
       width = 6, height = 3)







# ~~~~~~ Sea ice * day ---------------------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)


# Mean day capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df.for.plot.mean.day.capture <- df.for.plot.sea.ice %>%
  rename(day = AO)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * mean(day_s) +
                         b1cub[j, 6] * grid_scaled[i] * mean(day_s))
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * mean(day_s) +
                           b2_3cub[j, 6] * grid_scaled[i] * mean(day_s))
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

df.for.plot.mean.day.capture <- data.frame(var = grid,
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
         day = "mean")

# Early day capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fst_decile <- quantile(x = day_s, # = -0.4841
                       probs = 0.10, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * fst_decile +
                         b1cub[j, 6] * grid_scaled[i] * fst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * fst_decile +
                           b2_3cub[j, 6] * grid_scaled[i] * fst_decile)
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

df.for.plot.early.capture <- data.frame(var = grid,
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
         day = "early")



# Late capture day ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lst_decile <- quantile(x = day_s, # = 0.5719
                       probs = 0.90, 
                       na.rm = TRUE)

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(winter_AO_s) +
                         b1cub[j, 4] * mean(prior_spring_AO_s) +
                         b1cub[j, 5] * lst_decile +
                         b1cub[j, 6] * grid_scaled[i] * lst_decile)
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(winter_AO_s) +
                           b2_3cub[j, 4] * mean(prior_spring_AO_s) +
                           b2_3cub[j, 5] * lst_decile +
                           b2_3cub[j, 6] * grid_scaled[i] * lst_decile)
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

df.for.plot.late.capture <- data.frame(var = grid,
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
         day = "late")


df.for.plot.ice_free_facets_day <- rbind(df.for.plot.early.capture, 
                                         df.for.plot.mean.day.capture, 
                                         df.for.plot.late.capture) %>%
  mutate(day = factor(day, levels = c("early", "mean", "late")))

write_csv(df.for.plot.ice_free_facets_day, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_ice_free_days_facet_day_capture.csv")
df.for.plot.ice_free_facets_day <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_ice_free_days_facet_day_capture.csv") %>%
  mutate(day = factor(day, levels = c("early", "mean", "late")))


labels <- c("early capture", "mean", "late capture")
names(labels) <- c("early", "mean", "late")

ggplot(data = df.for.plot.ice_free_facets_day, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap( ~ day,
              labeller = labeller(day = labels))

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_ice_free_days_facet_day_capture.png",
       width = 10, height = 4)
ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2_ice_free_days_facet_day_capture.svg",
       width = 10, height = 4)





# 7. Ice-free days t-1 (com) * date of capture (dist) + spring AO t-1 (com) + winter AO t (com) -------------

# /!\ The coefficients for the interaction are distinct

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     winter_AO_s = as.numeric(winter_AO_s),
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s * day_s + winter_AO_s + prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 11, 2, 10, 12)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3 <- nimbleMCMC(code = model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3,     # model code  
                                                              data = dat,                                   
                                                              constants = my.constants,        
                                                              inits = inits,          
                                                              monitors = params,   # parameters to monitor
                                                              thin = 10,
                                                              niter = 100000,                  # nb iterations
                                                              nburnin = 10000,              # length of the burn-in
                                                              nchains = 2,
                                                              summary = TRUE,
                                                              WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3$WAIC
# 756.3347 

save(fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "b.9.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)





# 8. Ice-free days t-1 (com) * date of capture (dist) --------------------------

# Here the coefficient for the interaction is common

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s  * day_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 2, 6)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_common_4_distinct <- nimbleMCMC(code = model_1.1.2.D_common_4_distinct,     # model code  
                                                      data = dat,                                   
                                                      constants = my.constants,        
                                                      inits = inits,          
                                                      monitors = params,   # parameters to monitor
                                                      thin = 10,
                                                      niter = 100000,                  # nb iterations
                                                      nburnin = 10000,              # length of the burn-in
                                                      nchains = 2,
                                                      summary = TRUE,
                                                      WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_common_4_distinct$WAIC
# 756.0809

save(fit_1.1.2.D_common_4_distinct, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_common_4_distinct.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_common_4_distinct.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_common_4_distinct)


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_common_4_distinct.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)




# 9. Ice-free days t-1 (com) * date of capture (dist) + winter AO t (com) -------------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     winter_AO_s = as.numeric(winter_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s * day_s +  winter_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 2, 10)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2.D_2.1.1.1_common_4_distinct <- nimbleMCMC(code = model_1.1.2.D_2.1.1.1_common_4_distinct,     # model code  
                                                            data = dat,                                   
                                                            constants = my.constants,        
                                                            inits = inits,          
                                                            monitors = params,   # parameters to monitor
                                                            thin = 10,
                                                            niter = 100000,                  # nb iterations
                                                            nburnin = 10000,              # length of the burn-in
                                                            nchains = 2,
                                                            summary = TRUE,
                                                            WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_1.1.2.D_2.1.1.1_common_4_distinct$WAIC
# 758.0269 

save(fit_1.1.2.D_2.1.1.1_common_4_distinct, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_common_4_distinct.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.1.1.1_common_4_distinct.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2.D_2.1.1.1_common_4_distinct)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_1.1.2.D_2.1.1.1_common_4_distinct.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)




# 10. date of capture (dist) + spring AO t-1 (com) + winter AO t (com) -------------

# /!\ The coefficient for the interaction is common

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     winter_AO_s = as.numeric(winter_AO_s),
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| winter_AO_s + prior_spring_AO_s + day_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 2, 8)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_2.1.1.1_2.2.2.1_common_4_distinct <- nimbleMCMC(code = model_2.1.1.1_2.2.2.1_common_4_distinct,     # model code  
                                                              data = dat,                                   
                                                              constants = my.constants,        
                                                              inits = inits,          
                                                              monitors = params,   # parameters to monitor
                                                              thin = 10,
                                                              niter = 100000,                  # nb iterations
                                                              nburnin = 10000,              # length of the burn-in
                                                              nchains = 2,
                                                              summary = TRUE,
                                                              WAIC = TRUE)
end <- Sys.time()
end - start # 

fit_2.1.1.1_2.2.2.1_common_4_distinct$WAIC
# 759.7469

save(fit_2.1.1.1_2.2.2.1_common_4_distinct, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_2.1.1.1_2.2.2.1_common_4_distinct.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_2.1.1.1_2.2.2.1_common_4_distinct.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_2.1.1.1_2.2.2.1_common_4_distinct)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/model_2.1.1.1_2.2.2.1_common_4_distinct.RData.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)
