#==============================================================================#
#                                                                              #
#                Models with environmental covariates (binomial)               #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(nimble)
library(viridis)
library(cowplot)
Sys.setenv(LANG = "en")


# READY VARIABLES ===========================================================
# CR data
CR_data <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

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
  mutate(success = ifelse(cub_number == 1, 0, 1)) -> data_model





# ~~~ b. Format data for Nimble ------------------------------------------------

y <- data_model %>%
  pull(success)

# Individual traits +++++++++++++++++++++++++++++++++

# Size
size <- data_model$s_length
size_s <- (size - mean(size))/sd(size) 

# Size squarred
size_s2 <- size_s^2

# Age (in 3 categories)
{age <- as.numeric(data_model$age_for_analyses)
  age_factor <- as.factor(ifelse (age < 9, "b. young", 
                                  ifelse(age > 15, "c. old", "a. prime_age")))
  relevel(age_factor, ref = "a. prime_age")
  
  age_old <- ifelse (age > 15, 1, 0)
  age_young <- ifelse (age < 9, 1, 0)
}



# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- scale(day)


# Environmental variables +++++++++++++++++++++++++++++++++

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- scale(ice_free_days_previous)

# Winter AO
winter_AO <- data_model$winter_AO
winter_AO_s <- scale(winter_AO)

# Spring AO t-1 
prior_spring_AO <- data_model$prior_spring_AO
prior_spring_AO_s <- scale(prior_spring_AO)

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

# Load the Nimble models + the ancillary functions
source("05_script/models/several_predictors/2_Nimble_environmental_&_individual_covariates.R")
source("05_script/models/functions_for_models_Nimble.R")


# B. RUN THE MODELS ============================================================

# ~ 1. Ice-free days + day capture + individual ---------------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     day_s = as.numeric(day_s),
                     size_s = as.numeric(size_s),
                     size_s2 = as.numeric(size_s2),
                     age_young = age_young,
                     age_old = age_old)

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       size_s = as.numeric(size_s),
                       size_s2 = as.numeric(size_s2),
                       age_young = age_young,
                       age_old = age_old)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + day_s * (age_young + age_old) + size_s + size_s2, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5, 6, 7, 8, 9)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_individual_binomial <- nimbleMCMC(code = model_1.1.2_D_4_individual_binomial,     # model code  
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
end - start 

fit_1.1.2_D_4_individual_binomial$WAIC
# 310.7449

save(fit_1.1.2_D_4_individual_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "b.9.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_individual_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_individual_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)


# ~~~ c. Plot ------------------------------------------------------------------

load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial.RData")

N <- dim(fit_1.1.2_D_4_individual_binomial$samples$chain1)[1]
res <- rbind(fit_1.1.2_D_4_individual_binomial$samples$chain1[seq(1, N, by = 3), c(1:9)],
             fit_1.1.2_D_4_individual_binomial$samples$chain2[seq(1, N, by = 3), c(1:9)])

# ~~~~~~ Caterpillar plot

caterpilar <- as.data.frame(res[, c(1:9)]) %>%
  pivot_longer(cols = c("b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]", "b[9]"))

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

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_individual_binomial_caterpillar.svg",
       width = 5, height = 4)





# x axis: ice-free days 
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

p2_3cub <- matrix(data = NA, 
                  nrow = dim(res)[1], 
                  ncol = lengthgrid)

# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    p2_3cub[j, i] <- plogis(res[j, 1] + 
                              res[j, 2] * grid_scaled[i] + 
                              res[j, 3] * mean(day_s) +
                              # res[j, 4] * 0 +
                              # res[j, 5] * 0 +
                              res[j, 6] * mean(size_s) +
                              res[j, 7] * mean(size_s)^2)
                              # res[j, 8] * 0 * mean(day_s) +
                              # res[j, 9] * 0 * mean(day_s))
  }
}
df.for.plot.ice.free.days <- data.frame(var = grid,
                                        mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                        ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                        ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975))

ggplot(data = df.for.plot.ice.free.days) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Number of ice-free days",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "") 

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_individual_binomial_ice_free_days.svg",
       width = 5, height = 4)


ggplot(data = df.for.plot.ice.free.days) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Jours sans banquise l'année n-1",
       y = "Probabilité d'avoir >1 ourson",
       linetype = "") 


ggsave("D:/polar_bears_litter_size_environment/09_presentations/Mon PLR en 150 secondes/plot_bino.png",
       width = 2.75, height = 2.75)
# data.points <- data.frame(ice_free_days_previous = ice_free_days_previous,
#                           success = y)
# 
# ggplot() +
#   geom_line(data = df.for.plot.ice.free.days, aes(x = var, y = mean_p_2_3_cub)) +
#   geom_ribbon(data = df.for.plot.ice.free.days, aes(x = var, 
#                                                     ymin = ci_p_2_3_cub_2.5,
#                                                     ymax = ci_p_2_3_cub_97.5),
#               alpha = 0.15) +
#   theme_bw() +
#   labs(x = "Number of ice-free days",
#        y = "Probability of having 2 or 3 cubs", 
#        linetype = "")  +
#   geom_point(data = data.points, aes(x = ice_free_days_previous,
#                                      y = success))


############### Pour le PLR en 150 secondes

df.for.plot.ice.free.days <- data.frame(var = grid,
                                        mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                        ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                        ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975))












# x axis: day of capture
range <- range(day_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(day) + mean(day)

p2_3cub <- matrix(data = NA, 
                  nrow = dim(res)[1], 
                  ncol = lengthgrid)

# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
       p2_3cub[j, i] <- plogis(res[j, 1] + 
                                   res[j, 2] * mean(ice_free_days_previous_s) +
                                   res[j, 3] * grid_scaled[i] + 
                                   # res[j, 4] * 0 +
                                   # res[j, 5] * 0 +
                                   res[j, 6] * mean(size_s) +
                                   res[j, 7] * mean(size_s)^2)
    # res[j, 8] * 0 * mean(day_s) +
    # res[j, 9] * 0 * mean(day_s))
  }
}
df.for.plot.day.capture <- data.frame(var = grid,
                                      mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                      ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                      ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975))


ggplot(data = df.for.plot.day.capture) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Day of capture",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "") 

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_individual_binomial_day_capture.svg",
       width = 5, height = 4)



# x axis: size
range <- range(size_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(size) + mean(size)

p2_3cub <- matrix(data = NA, 
                  nrow = dim(res)[1], 
                  ncol = lengthgrid)

# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    p2_3cub[j, i] <- plogis(res[j, 1] + 
                              res[j, 2] * mean(ice_free_days_previous_s) +
                              res[j, 3] * mean(day_s) + 
                              # res[j, 4] * 0 +
                              # res[j, 5] * 0 +
                              res[j, 6] * grid_scaled[i] +
                              res[j, 7] * grid_scaled[i]^2)
    # res[j, 8] * 0 * mean(day_s) +
    # res[j, 9] * 0 * mean(day_s))
  }
}
df.for.plot.size <- data.frame(var = grid,
                                      mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                      ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                      ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975))




ggplot(data = df.for.plot.size) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Size",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "") 

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_individual_binomial_size.svg",
       width = 5, height = 4)





# x axis: day capture facetted
range <- range(day_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(day) + mean(day)

p2_3cub <- matrix(data = NA, 
                  nrow = dim(res)[1], 
                  ncol = lengthgrid)

# Mean Size
# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    p2_3cub[j, i] <- plogis(res[j, 1] + 
                              res[j, 2] * mean(ice_free_days_previous_s) +
                              res[j, 3] * grid_scaled[i] + 
                              # res[j, 4] * 0 +
                              # res[j, 5] * 0 +
                              res[j, 6] * mean(size_s) +
                              res[j, 7] * mean(size_s)^2)
    # res[j, 8] * 0 * mean(day_s) +
    # res[j, 9] * 0 * mean(day_s))
  }
}
df.for.plot.day.capture.mean <- data.frame(var = grid,
                                      mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                      ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                      ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975),
                                      age = "mean")


# Small Size
# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    p2_3cub[j, i] <- plogis(res[j, 1] + 
                              res[j, 2] * mean(ice_free_days_previous_s) +
                              res[j, 3] * grid_scaled[i] + 
                              res[j, 4] * 1 +
                              # res[j, 5] * 0 +
                              res[j, 6] * mean(size_s) +
                              res[j, 7] * mean(size_s)^2 +
                              res[j, 8] * 1 * grid_scaled[i])
    # res[j, 9] * 0 * grid_scaled[i])
  }
}
df.for.plot.day.capture.young <- data.frame(var = grid,
                                           mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                           ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                           ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975),
                                           age = "young")



# Large Size
# Back transform
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    p2_3cub[j, i] <- plogis(res[j, 1] + 
                              res[j, 2] * mean(ice_free_days_previous_s) +
                              res[j, 3] * grid_scaled[i] + 
                              # res[j, 4] * 0 +
                              res[j, 5] * 1 +
                              res[j, 6] * mean(size_s) +
                              res[j, 7] * mean(size_s)^2 +
                              # res[j, 8] * 0 * grid_scaled[i])
                              res[j, 9] * 1 * grid_scaled[i])
  }
}
df.for.plot.day.capture.old <- data.frame(var = grid,
                                            mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                            ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                            ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975),
                                            age = "old")

df.for.plot.day.capture.facet.age <- rbind(df.for.plot.day.capture.young, 
                                            df.for.plot.day.capture.mean,
                                            df.for.plot.day.capture.old) %>%
  mutate(age = factor(age, levels = c("young", "mean", "old")))
  

labels <- c("young", "prime-aged", "old")
names(labels) <- c("young", "mean", "old")

ggplot(data = df.for.plot.day.capture.facet.age) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Day of capture",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "") +
  facet_wrap(.~age,
             labeller = labeller(age = labels))

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_individual_binomial_day_capture_facet_age.svg",
       width = 10, height = 4)







# ~ 2. Ice-free days + day capture + individual without interaction ------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     day_s = as.numeric(day_s),
                     size_s = as.numeric(size_s),
                     size_s2 = as.numeric(size_s2),
                     age_young = age_young,
                     age_old = age_old)

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       size_s = as.numeric(size_s),
                       size_s2 = as.numeric(size_s2),
                       age_young = age_young,
                       age_old = age_old)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + day_s + age_young + age_old + size_s + size_s2, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5, 6, 7)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_individual_binomial_2 <- nimbleMCMC(code = model_1.1.2_D_4_individual_binomial_2,     # model code  
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
end - start 

fit_1.1.2_D_4_individual_binomial_2$WAIC
# 313.0591

save(fit_1.1.2_D_4_individual_binomial_2, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_2.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_2.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_individual_binomial_2)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_individual_binomial_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)




# ~ 3. Ice-free days + day capture + individual without quadratic ------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     day_s = as.numeric(day_s),
                     size_s = as.numeric(size_s),
                     age_young = age_young,
                     age_old = age_old)

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       size_s = as.numeric(size_s),
                       age_young = age_young,
                       age_old = age_old)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + day_s * (age_young + age_old) + size_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5, 6, 7, 8)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_individual_binomial_3 <- nimbleMCMC(code = model_1.1.2_D_4_individual_binomial_3,     # model code  
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
end - start 

fit_1.1.2_D_4_individual_binomial_3$WAIC
# 311.3152

save(fit_1.1.2_D_4_individual_binomial_3, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_3.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_3.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_individual_binomial_3)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_individual_binomial_3.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)




# ~ 4. Ice-free days + day capture + individual (without interac and quadratic) ---------------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     day_s = as.numeric(day_s),
                     size_s = as.numeric(size_s),
                     age_young = age_young,
                     age_old = age_old)

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       size_s = as.numeric(size_s),
                       age_young = age_young,
                       age_old = age_old)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + day_s + age_young + age_old + size_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5, 6)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_individual_binomial_4 <- nimbleMCMC(code = model_1.1.2_D_4_individual_binomial_4,     # model code  
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
end - start 

fit_1.1.2_D_4_individual_binomial_4$WAIC
# 310.7449

save(fit_1.1.2_D_4_individual_binomial_4, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_4.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_4.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_individual_binomial_4)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_individual_binomial_4.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)





# ~ 5. Day capture + individual (no ice-free days) -----------------------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     day_s = as.numeric(day_s),
                     size_s = as.numeric(size_s),
                     size_s2 = as.numeric(size_s2),
                     age_young = age_young,
                     age_old = age_old)

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       size_s = as.numeric(size_s),
                       size_s2 = as.numeric(size_s2),
                       age_young = age_young,
                       age_old = age_old)

mylogit <- glm(as.factor(y) ~ day_s * (age_young + age_old) + size_s + size_s2, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5, 6, 7, 8)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_individual_binomial_5 <- nimbleMCMC(code = model_1.1.2_D_4_individual_binomial_5,     # model code  
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
end - start 

fit_1.1.2_D_4_individual_binomial_5$WAIC
# 310.2714

save(fit_1.1.2_D_4_individual_binomial_5, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_5.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_individual_binomial_5.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.", "b.8.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_individual_binomial_5)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_individual_binomial_5.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)
