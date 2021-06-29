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
source("05_script/models/several_predictors/1_Nimble_environmental_covariates.R")
source("05_script/models/functions_for_models_Nimble.R")


# B. RUN THE MODELS ============================================================

# ~ 1. Ice-free days + day capture + winter AO + spring AO ---------------------

my.constants <- list(N = length(y), # nb of females captured
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

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + winter_AO_s + prior_spring_AO_s + day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4, 5)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial <- nimbleMCMC(code = model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial,     # model code  
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

fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial$WAIC
# 317.2758

save(fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)


# ~ 2. Ice-free days * day capture + winter AO + spring AO ---------------------

# Add the interaction between ice free days and day of capture

my.constants <- list(N = length(y), # nb of females captured
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

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s * day_s + winter_AO_s + prior_spring_AO_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 4, 5, 3, 6)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2 <- nimbleMCMC(code = model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2,     # model code  
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
end - start 

fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2$WAIC
# 319.3173

save(fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2)
diagnostic_plot

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)



# ~ 3. Ice-free days + day capture + spring AO ---------------------------------

# Remove winter AO

my.constants <- list(N = length(y), # nb of females captured
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

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + prior_spring_AO_s + day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_2.2.2.1_4_binomial <- nimbleMCMC(code = model_1.1.2_D_2.2.2.1_4_binomial,     # model code  
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

fit_1.1.2_D_2.2.2.1_4_binomial$WAIC
# 315.2216

save(fit_1.1.2_D_2.2.2.1_4_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.2.2.1_4_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.2.2.1_4_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_2.2.2.1_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_2.2.2.1_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)


# ~ 4. Ice-free days + day capture + winter AO ---------------------

# Remove spring AO t-1

my.constants <- list(N = length(y), # nb of females captured
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
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + winter_AO_s + day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))

# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_2.1.1.1_4_binomial <- nimbleMCMC(code = model_1.1.2_D_2.1.1.1_4_binomial,     # model code  
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

fit_1.1.2_D_2.1.1.1_4_binomial$WAIC
# 315.4223

save(fit_1.1.2_D_2.1.1.1_4_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_4_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_4_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_2.1.1.1_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_2.1.1.1_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)

# ~ 5. Ice-free days + winter AO + spring AO ---------------------

# Remove the day of capture

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     winter_AO_s = as.numeric(winter_AO_s),
                     prior_spring_AO_s = as.numeric(prior_spring_AO_s))

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + winter_AO_s + prior_spring_AO_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_2.1.1.1_2.2.2.1_binomial <- nimbleMCMC(code = model_1.1.2_D_2.1.1.1_2.2.2.1_binomial,     # model code  
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

fit_1.1.2_D_2.1.1.1_2.2.2.1_binomial$WAIC
# 318.0037

save(fit_1.1.2_D_2.1.1.1_2.2.2.1_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_2.1.1.1_2.2.2.1_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_2.1.1.1_2.2.2.1_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_2.1.1.1_2.2.2.1_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)

# ~ 6. Winter AO + spring AO + day of capture -----------------------------------

my.constants <- list(N = length(y), # nb of females captured
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

mylogit <- glm(as.factor(y) ~ winter_AO_s + prior_spring_AO_s + day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_2.1.1.1_2.2.2.1_4_binomial <- nimbleMCMC(code = model_2.1.1.1_2.2.2.1_4_binomial,     # model code  
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

fit_2.1.1.1_2.2.2.1_4_binomial$WAIC
# 316.8548

save(fit_2.1.1.1_2.2.2.1_4_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_2.1.1.1_2.2.2.1_4_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_2.1.1.1_2.2.2.1_4_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_2.1.1.1_2.2.2.1_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_2.1.1.1_2.2.2.1_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)



# ~ 7. Ice-free days + day capture ---------------------------------------------

# Remove winter AO and spring AO t-1

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       prior_spring_AO_s = prior_spring_AO_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s + day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))

# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_binomial <- nimbleMCMC(code = model_1.1.2_D_4_binomial,     # model code  
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

fit_1.1.2_D_4_binomial$WAIC
# 313.6873


save(fit_1.1.2_D_4_binomial, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_binomial.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)


# ~~~ c. Plot ------------------------------------------------------------------

load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_binomial.RData")

N <- dim(fit_1.1.2_D_4_binomial$samples$chain1)[1]
res <- rbind(fit_1.1.2_D_4_binomial$samples$chain1[seq(1, N, by = 3), c(1:3)],
             fit_1.1.2_D_4_binomial$samples$chain2[seq(1, N, by = 3), c(1:3)])

# ~~~~~~ Caterpillar plot

caterpilar <- as.data.frame(res[, c(1:3)]) %>%
  pivot_longer(cols = c("b[1]", "b[2]", "b[3]"))

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

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_binomial_caterpillar.svg",
       width = 3, height = 2)





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
    p2_3cub[j, i] <- plogis(res[j, 1] + res[j, 2] * grid_scaled[i] + res[j, 3] * mean(day_s))
    # p2_3cub[j, i] <- exp(res[j, 1] + res[j, 2] * grid_scaled[i])/(1 + exp(res[j, 1] + res[j, 2] * grid_scaled[i]))
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

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_binomial_ice_free_days.svg",
       width = 5, height = 4)


data.points <- data.frame(ice_free_days_previous = ice_free_days_previous,
                          success = y)

ggplot() +
  geom_line(data = df.for.plot.ice.free.days, aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(data = df.for.plot.ice.free.days, aes(x = var, 
                                                    ymin = ci_p_2_3_cub_2.5,
                                                    ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Number of ice-free days",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "")  +
  geom_point(data = data.points, aes(x = ice_free_days_previous,
                                     y = success))



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
    p2_3cub[j, i] <- plogis(res[j, 1] + res[j, 2] * mean(ice_free_days_previous_s) + res[j, 3] * grid_scaled[i])
    # p2_3cub[j, i] <- exp(res[j, 1] + res[j, 2] * grid_scaled[i])/(1 + exp(res[j, 1] + res[j, 2] * grid_scaled[i]))
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

ggsave("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1.1.2_D_4_binomial_day_capture.svg",
       width = 5, height = 4)

data.points <- data.frame(day = day,
                          success = y)

ggplot() +
  geom_line(data = df.for.plot.day.capture, aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(data = df.for.plot.day.capture, aes(x = var, 
                                                    ymin = ci_p_2_3_cub_2.5,
                                                    ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = "Number of ice-free days",
       y = "Probability of having 2 or 3 cubs", 
       linetype = "")  +
  geom_point(data = data.points, aes(x = day,
                                       y = success))





# ~ 8. Ice-free days * day capture ---------------------------------------------

# Remove winter AO and spring AO t-1

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s),
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s * day_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2, 3, 4)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))

# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_4_binomial_2 <- nimbleMCMC(code = model_1.1.2_D_4_binomial_2,     # model code  
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

fit_1.1.2_D_4_binomial_2$WAIC
# 315.6758


save(fit_1.1.2_D_4_binomial_2, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_binomial_2.RData")


# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.1.2_D_4_binomial_2.RData")
params.plot <- c("b.1.", "b.2.", "b.3.", "b.4.", "sigma1")


# Run the function
diagnostic_plot <- check_convergence_several_predictors(params.plot = params.plot, 
                                                        nimble_output = fit_1.1.2_D_4_binomial)

# Save the plot
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/several_predictors/graphs/diagnostic_plots/fit_1.1.2_D_4_binomial.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

rm(diagnostic_plot)

