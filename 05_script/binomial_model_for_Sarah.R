#==============================================================================#
#                                                                              #
#                          Binomial model for Sarah                           #
#                                                                              #
#==============================================================================#


library(tidyverse)
library(mlogit)
library(viridis)
library(ggmcmc)
library(gridExtra)
library(nimble)
Sys.setenv(LANG = "en")

# READY VARIABLES ===========================================================
# CR data
CR_data <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

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

# Response variable
y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Number of litters per year
n <- data_model %>%
  mutate(ones = 1) %>%
  group_by(year) %>%
  summarise(n = sum(ones)) %>%
  left_join(x = data_model,
            y = .,
            by = "year") %>%
  pull(n)


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


# RUN MODELS ===================================================================

# A. Null model ================================================================

# ~~~ a. Run the model ---------------------------------------------------------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     n = n) 

params <- c("b0", "sigma1", "eps1") 

# Find initial value for the intercept
binomial.frequ <- glm(y ~ 1, 
                      family = "binomial", 
                      data = data.frame(y = as.factor(y)))

coefs <- as.vector(summary(binomial.frequ)$coefficients[1])
inits_null <- function() list(b0 = coefs + round(runif(n = 1, -1, 1))/10)


start <- Sys.time()
fit_0.0.0.0_binomial <- nimbleMCMC(code = model_0.0.0.0_binomial,     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits_null,          
                          monitors = params,   # parameters to monitor
                          niter = 20000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE)
end <- Sys.time()
end - start

fit_0.0.0.0_binomial$WAIC

save(list = fit_0.0.0.0_binomial, 
     file = "07_results/model_0.0.0.0_binomial.RData")



# ~~~ b. Check convergence -----------------------------------------------------
nimble_output <- fit_0.0.0.0_binomial

# Process Nimble output into dataframe
chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
  dplyr::select(params[-length(params)]) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
  dplyr::select(params[-length(params)]) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params[-length(params)], names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarize(m = mean(value))

param.running.mean <- chains_l %>%
  arrange(parameter, iteration) %>%
  group_by(parameter, chain) %>%
  mutate(rm = cumsum(value)/iteration)

trace.plots <- ggplot(data = chains_l, 
                      aes(x = iteration, y = value, color = chain)) +
  geom_line() +
  labs(y = "trace") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
              ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free_y",
              ncol = 1)

running.mean.plot <- ggplot(param.running.mean, 
                            aes(x = iteration, y = rm, color = chain)) + 
  geom_line() + 
  geom_hline(aes(yintercept = m), param.mean,
             colour = "black", alpha = 0.5) + 
  ylab("running Mean") +
  facet_grid(parameter ~ chain, scales = "free")

# Plot all the plots together
diagnostic_plot <- plot_grid(trace.plots,
                             density.plots, 
                             running.mean.plot,
                             ncol = 3, nrow = 1)
diagnostic_plot




# B. Model with only ice-free days =============================================

# ~~~ a. Run the model ---------------------------------------------------------

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"

# Are females without cubs taken into account ?

my.constants <- list(N = length(y),            # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     n = n,                    # nb of females captured the same year
                     as.numeric(var_scaled)) 
names(my.constants)[6] <- var_short_name


# Define the parameters to estimate
params <- c("b0", "b1", "sigma1", "eps1")
temp.data <- data.frame(y = y, 
                        var_scaled = var_scaled)
mylogit <- glm(as.factor(y) ~ var_scaled, data = temp.data, family = "binomial")

inits <- function() list(b0 = mylogit$coefficients[1] + round(runif(n = 1, -1, 1))/10, 
                         b1 = mylogit$coefficients[2] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()

fit_1.1.2_D_binomial <- nimbleMCMC(code = model_1.1.2_D_binomial,     # model code  
                                   data = dat,                                   
                                   constants = my.constants,        
                                   inits = inits,          
                                   monitors = params,   # parameters to monitor
                                   thin = 10,
                                   niter = 15000,                  # nb iterations
                                   nburnin = 5000,              # length of the burn-in
                                   nchains = 2,
                                   summary = TRUE,
                                   WAIC = TRUE)
end <- Sys.time()
end - start


fit_1.1.2_D_binomial$WAIC

save(list = fit_1.1.2_D_binomial, 
     file = "07_results/model_1.1.2_D_binomial.RData")


# ~~~ b. Check convergence -----------------------------------------------------

nimble_output <- fit_1.1.2_D_binomial

# Process Nimble output into dataframe
chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
  dplyr::select(params[-length(params)]) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
  dplyr::select(params[-length(params)]) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params[-length(params)], names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarize(m = mean(value))

param.running.mean <- chains_l %>%
  arrange(parameter, iteration) %>%
  group_by(parameter, chain) %>%
  mutate(rm = cumsum(value)/iteration)

trace.plots <- ggplot(data = chains_l, 
                      aes(x = iteration, y = value, color = chain)) +
  geom_line() +
  labs(y = "trace") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
              ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free_y",
              ncol = 1)

running.mean.plot <- ggplot(param.running.mean, 
                            aes(x = iteration, y = rm, color = chain)) + 
  geom_line() + 
  geom_hline(aes(yintercept = m), param.mean,
             colour = "black", alpha = 0.5) + 
  ylab("running Mean") +
  facet_grid(parameter ~ chain, scales = "free")

# Plot all the plots together
diagnostic_plot <- plot_grid(trace.plots,
                             density.plots, 
                             running.mean.plot,
                             ncol = 3, nrow = 1)
diagnostic_plot
