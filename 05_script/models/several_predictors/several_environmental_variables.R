#==============================================================================#
#                                                                              #
#                      Models with ice-free days + AO/NAO                      #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(nimble)
library(viridis)
library(ggmcmc)
library(gridExtra)
library(MCMCvis)
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
  filter(year != 1992)  %>% # Remove the captures from 1992 since I can't calculate the ice free days in 1991
  filter(year != 1993)  -> data_model

rm(CR_data, sea_ice_data, AO_data, NAO_data)



# Plot correlations

# Sea ice - sea ice
ggplot(data_model, aes(x = ice_free_days, y = day_retreat)) +
  geom_point()
summary(lm(ice_free_days ~ day_retreat, data = data_model))

ggplot(data_model, aes(x = ice_free_days, y = day_advance)) +
  geom_point()
summary(lm(ice_free_days ~ day_advance, data = data_model))

ggplot(data_model, aes(x = day_retreat, y = day_advance)) +
  geom_point()
summary(lm(day_retreat ~ day_advance, data = data_model))

# Sea ice - AO
ggplot(data_model, aes(x = ice_free_days, y = spring_AO)) +
  geom_point()
summary(lm(ice_free_days ~ spring_AO, data = data_model))

ggplot(data_model, aes(x = ice_free_days, y = winter_AO)) +
  geom_point()
summary(lm(ice_free_days ~ winter_AO, data = data_model))

# Sea ice - NAO
ggplot(data_model, aes(x = ice_free_days, y = spring_NAO)) +
  geom_point()
summary(lm(ice_free_days ~ spring_NAO, data = data_model))

ggplot(data_model, aes(x = ice_free_days, y = winter_NAO)) +
  geom_point()
summary(lm(ice_free_days ~ winter_NAO, data = data_model))

# AO - AO
ggplot(data_model, aes(x = spring_AO, y = winter_AO)) +
  geom_point()
summary(lm(spring_AO ~ winter_AO, data = data_model))

# NAO - NAO
ggplot(data_model, aes(x = spring_NAO, y = winter_NAO)) +
  geom_point()
summary(lm(spring_NAO ~ winter_NAO, data = data_model))



# ~~~ b. Format data for Nimble ------------------------------------------------

y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- (day - mean(day))/sd(day)

# Number of ice-free days t-1

ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- scale(ice_free_days_previous)

# Number of ice-free days t-1

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

# females
{id_fem <- data_model$ID_NR
  id_fem <- factor(id_fem) 
  id_fem2 <- NULL
  for (i in 1:length(id_fem)){
    id_fem2 <- c(id_fem2,which(id_fem[i]==levels(id_fem)))
  }
  id_fem <- factor(id_fem2)
  # nbind <- length(levels(id_fem))
  rm(id_fem2)
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/several_predictors/Nimble_several_envrionmental_predictors.R")
source("05_script/models/functions_for_models.R")

# B. RUN THE MODELS ============================================================

# 1. Previous year ice-free days + individual variables: common slopes ---------

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


inits <- function() list(b = coefs[c(1, 3, 5, 7, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_model_1.1.2.D_2.2.2.1_common <- nimbleMCMC(code = model_1.1.2.D_2.2.2.1_common,     # model code  
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
end - start # 24min

fit_model_1.1.2.D_2.2.2.1_common$WAIC
# 1073.887

save(fit_model_1.1.2.D_2.2.2.1_common, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "sigma1")

# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_model_1.1.2.D_2.2.2.1_common)
diagnostic_plot


# Save the plot
nrows = length(params.plot)
save_plot(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_model_1.1.2.D_2.2.2.1_common.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2.D_2.2.2.1_common.RData")

N <- dim(fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]])[1]
res <- rbind(fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain1"]][seq(1, N, by = 4), 
                                                  c(1:5)],
             fit_model_1.1.2.D_2.2.2.1_common[["samples"]][["chain2"]][seq(1, N, by = 4), 
                                                  c(1:5)])
b1cub <- res[, c(1:4)]
b2_3cub <- res[, c(5, 2:4)] 

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
  mutate(AO = factor(AO, levels = c("low AO index", "mean AO index", "high AO index")))

write_csv(df.for.plot.spring_AO, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_ice_free_days_facet_spring_AO.csv")
df.for.plot.spring_AO <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2.D_2.2.2.1_common_ice_free_days_facet_spring_AO.csv") %>%
  mutate(AO = factor(AO, levels = c("low AO index", "mean AO index", "high AO index")))

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
  facet_wrap( ~ AO)

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
