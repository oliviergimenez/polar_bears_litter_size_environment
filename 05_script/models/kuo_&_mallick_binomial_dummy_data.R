#==============================================================================#
#                                                                              #
#                Test of covariate selection with dummy data                   #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(nimble)
library(viridis)
library(cowplot)
Sys.setenv(LANG = "en")


# READY VARIABLES ==============================================================

CR_data <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

sea_ice_data <- read_csv("06_processed_data/sea_ice_data/SI_metrics_D.csv")
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


# RUN THE MODELS ===============================================================

# ~ 1. Day capture + Ice-free days with indicator ------------------------------

# ~~~ a. Generate the data -----------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))


# ~~~ b. Generate the data
set.seed(666)
z <- 1 - 0.5*day_number_s - 2*ice_free_days_previous_s # + rnorm(1, 0, 1)       # linear combination 

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

source("05_script/models/functions_for_models_Nimble.R")


# ~~~ b. Build the model -------------------------------------------------------

model_dummy_data_test <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      beta.day*day_number_s[i] +
      w.ice*beta.ice*ice_free_days_previous_s[i] 
  } # 'i'
  
  
  beta0 ~ dnorm(0, V.mu)
  beta.day ~ dnorm(0, V.beta.day)
  beta.ice ~ dnorm(0, V.beta.ice)		
  
  V ~ dunif(0,100)
  
  ## Equations ##
  V.mu <- V/nbr_param_weight[M, 1]
  V.beta.day <- V/nbr_param_weight[M, 2]
  V.beta.ice <- V/nbr_param_weight[M, 3]
  
  
  M ~ dcat(p.M[])
  
  w.ice <- equals(M, 2) # w.alpha <- 1 if M = 2, otherwise, w.alpha = 0
  
  for (i in 1:dim.M){
    indicator_model[i] <- equals(M,i)
  }
})


# ~~~ c. Constants and data ----------------------------------------------------

# Number of Parameter weight for each Model
nbr_param_weight <- matrix(data = c(0.5, 0.3333333, 0.5, 0.3333333, 1, 0.3333333),
                           nrow = 2)
nbr_param_weight

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)


# Bundle data
my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            nbr_param_weight = nbr_param_weight)


# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         beta.day = rnorm(1, 0, 10),
                         beta.ice = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model", "w.ice",
            "M") 

# ~~~ d. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_dummy_data_test <-  nimbleMCMC(code = model_dummy_data_test,     # model code  
                                  data = dat,                                   
                                  constants = my.constants,        
                                  inits = inits,          
                                  monitors = params,   # parameters to monitor
                                  niter = 100000,                  # nb iterations
                                  nburnin = 10000,              # length of the burn-in
                                  nchains = 2,
                                  summary = TRUE,
                                  WAIC = TRUE)
end <- Sys.time()
end - start


# ~~~ e. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model.1.",
                 "indicator_model.2.", "w.ice", "M") 

subsample <- seq(from = 1, to = dim(fit_dummy_data_test[["samples"]][["chain1"]])[1], by = 2)


chain1 <- data.frame(fit_dummy_data_test[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_dummy_data_test[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_dummy_data_test[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_dummy_data_test[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarise(m = mean(value))

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

# density.plots <- ggplot(data = chains_l, 
#                         aes(x = value, color = chain, fill = chain)) +
#   geom_density(alpha = 0.25) +
#   labs(x = "density") +
#   theme(legend.position = "none") +
#   facet_wrap( ~ parameter,
#               scales = "free",
#               ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_histogram(alpha = 0.25,
                 position="identity") +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
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

nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_dummy_data_test.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)



# ~ 2. Day capture  with indicator + Ice-free days with indicator --------------

# ~~~ a. Generate the data -----------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# ~~~ b. Generate the data
set.seed(666)
z <- 1 - 0.5*day_number_s - 2*ice_free_days_previous_s# + rnorm(1, 0, 1)       # linear combination 

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

source("05_script/models/functions_for_models_Nimble.R")


# ~~~ b. Build the model -------------------------------------------------------

model_dummy_data_test_2 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * beta.day * day_number_s[i] +
      w[2] * beta.ice * ice_free_days_previous_s[i] 
  } # 'i'
  
  
  beta0 ~ dnorm(0, V.mu)
  beta.day ~ dnorm(0, V.beta.day)
  beta.ice ~ dnorm(0, V.beta.ice)		
  
  V ~ dunif(0,100)
  
  ## Equations ##
  V.mu <- V/nbr_param_weight[M, 1]
  V.beta.day <- V/nbr_param_weight[M, 2]
  V.beta.ice <- V/nbr_param_weight[M, 3]
  
  
  M ~ dcat(p.M[])
  
  for (i in 1:2) {
    w[i] <- values_w[M, i]
  }
  
  for (i in 1:dim.M){
    indicator_model[i] <- equals(M, i)
  }
})

# ~~~ c. Constants and data ----------------------------------------------------

nbr_param_weight <- as.matrix(rbind(c(1, 1, 1),
                                    c(0.5, 0.5, 0.5),
                                    c(0.5, 0.5, 0.5),
                                    c(0.3333333, 0.3333333, 0.3333333)))
nbr_param_weight

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)

# Values of w[i]
values_w <- as.matrix(rbind(c(0, 0),
                            c(1, 0),
                            c(0, 1),
                            c(1, 1)))


# Bundle data
my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            nbr_param_weight = nbr_param_weight,
            values_w = values_w)
#values_indicator = values_indicator)


# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         beta.day = rnorm(1, 0, 10),
                         beta.ice = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model", "w[1]", "w[2]",
            "M") 

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_dummy_data_test_2 <-  nimbleMCMC(code = model_dummy_data_test_2,     # model code  
                                  data = dat,                                   
                                  constants = my.constants,        
                                  inits = inits,          
                                  monitors = params,   # parameters to monitor
                                  niter = 100000,                  # nb iterations
                                  nburnin = 10000,              # length of the burn-in
                                  nchains = 2,
                                  summary = TRUE,
                                  WAIC = TRUE)
end <- Sys.time()
end - start

fit_dummy_data_test_2$summary$all.chains

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "w.1.", "w.2.", "M",
                 "indicator_model.1.", "indicator_model.2.", 
                 "indicator_model.3.", "indicator_model.4.") 

subsample <- seq(from = 1, to = dim(fit_dummy_data_test_2[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_dummy_data_test_2[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_dummy_data_test_2[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_dummy_data_test_2[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_dummy_data_test_2[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarise(m = mean(value))

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

# density.plots <- ggplot(data = chains_l, 
#                         aes(x = value, color = chain, fill = chain)) +
#   geom_density(alpha = 0.25) +
#   labs(x = "density") +
#   theme(legend.position = "none") +
#   facet_wrap( ~ parameter,
#               scales = "free",
#               ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_histogram(alpha = 0.25,
                 position="identity") +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
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

nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_dummy_data_test_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)




# ~ 3. Day capture with indicator + Ice-free days with indicator + size with indicator --------

# ~~~ a. Generate the data -----------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

size <- data_model$s_length
size_s <- as.vector(scale(size))

# ~~~ b. Generate the data
set.seed(666)
z <- 1 - 0.5*day_number_s - 2*ice_free_days_previous_s + 0.75*size# + rnorm(1, 0, 1)       # linear combination 

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

# ~~~ b. Build the model ------------------------------------------------------

model_dummy_data_test_3 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * beta.day * day_number_s[i] +
      w[2] * beta.ice * ice_free_days_previous_s[i] +
      w[3] * beta.size * size_s[i]
  } # 'i'
  
  
  beta0 ~ dnorm(0, V.mu)
  beta.day ~ dnorm(0, V.beta.day)
  beta.ice ~ dnorm(0, V.beta.ice)		
  beta.size ~ dnorm(0, V.beta.size)	
  
  V ~ dunif(0,100)
  
  ## Equations ##
  V.mu <- V/nbr_param_weight[M, 1]
  V.beta.day <- V/nbr_param_weight[M, 2]
  V.beta.ice <- V/nbr_param_weight[M, 3]
  V.beta.size <- V/nbr_param_weight[M, 4]
  
  
  M ~ dcat(p.M[])
  
  for (i in 1:3) {
    w[i] <- values_w[M, i]
  }
  
  for (i in 1:dim.M){
    indicator_model[i] <- equals(M, i)
  }
  
  # for (i in 1:dim.M){
  #   indicator_model[i] <- values_indicator[M, i]
  # }
})



# ~~~ c. Constants and data ----------------------------------------------------

# Number of Parameter weight for each Model
nbr_param_weight <- as.matrix(rbind(c(1, 1, 1, 1),
                                    c(0.5, 0.5, 0.5, 0.5),
                                    c(0.5, 0.5, 0.5, 0.5),
                                    c(0.5, 0.5, 0.5, 0.5),
                                    c(0.3333333, 0.3333333, 0.3333333, 0.3333333),
                                    c(0.3333333, 0.3333333, 0.3333333, 0.3333333),
                                    c(0.3333333, 0.3333333, 0.3333333, 0.3333333),
                                    c(0.25, 0.25, 0.25, 0.25)))
nbr_param_weight

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)

# Values of w[i]
values_w <- as.matrix(rbind(c(0, 0, 0),
                            c(1, 0, 0),
                            c(0, 1, 0),
                            c(0, 0, 1),
                            c(1, 1, 0),
                            c(1, 0, 1),
                            c(0, 1, 1),
                            c(1, 1, 1)))

# Bundle data
my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     size_s = size_s,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            nbr_param_weight = nbr_param_weight,
            values_w = values_w)


# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         beta.day = rnorm(1, 0, 10),
                         beta.ice = rnorm(1, 0, 10),
                         beta.size = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "beta.day", "beta.ice", "beta.size", "V", "indicator_model", "w", # "w[1]", "w[2]",
            "M") 

# ~~~ d. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_dummy_data_test_3 <-  nimbleMCMC(code = model_dummy_data_test_3,     # model code  
                                  data = dat,                                   
                                  constants = my.constants,        
                                  inits = inits,          
                                  monitors = params,   # parameters to monitor
                                  niter = 150000,                  # nb iterations
                                  nburnin = 10000,              # length of the burn-in
                                  nchains = 2,
                                  summary = TRUE,
                                  WAIC = TRUE)
end <- Sys.time()
end - start

fit_dummy_data_test_3$summary$all.chains


# ~~~ e. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "beta.size", "V", 
                 "indicator_model.1.", "indicator_model.2.", "indicator_model.3.",
                 "indicator_model.4.", "indicator_model.5.", "indicator_model.6.",
                 "indicator_model.7.", "indicator_model.8.",
                 "w.1.", "w.2.", "w.3.", "M") 

subsample <- seq(from = 1, to = dim(fit_dummy_data_test_3[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_dummy_data_test_3[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_dummy_data_test_3[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_dummy_data_test_3[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_dummy_data_test_3[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarise(m = mean(value))

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

# density.plots <- ggplot(data = chains_l, 
#                         aes(x = value, color = chain, fill = chain)) +
#   geom_density(alpha = 0.25) +
#   labs(x = "density") +
#   theme(legend.position = "none") +
#   facet_wrap( ~ parameter,
#               scales = "free",
#               ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_histogram(alpha = 0.25,
                 position="identity") +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
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
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_dummy_data_test_3.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)



# ~ 4. With 4 indicators without cnst variance --------------------------------

# ~~~ a. Generate the data -----------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# Age
{age <- as.numeric(data_model$age_for_analyses)
  age_factor <- as.factor(ifelse (age < 9, "b. young", 
                                  ifelse(age > 15, "c. old", "a. prime_age")))
  relevel(age_factor, ref = "a. prime_age")
  
  age_old <- ifelse (age > 15, 1, 0)
  age_young <- ifelse (age < 9, 1, 0)
}

# ~~~ b. Generate the data
set.seed(666)
z <- 1 - 0.5*day_number_s - 2*ice_free_days_previous_s - 0.75*age_young - 0.75*age_old

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

# ~~~ b. Build the model -------------------------------------------------------

model_dummy_data_test_4 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * beta.day * day_number_s[i] +
      w[2] * beta.ice * ice_free_days_previous_s[i] +
      w[3] * beta.age.y * age_young[i] +
      w[4] * beta.age.o * age_old[i]
  } # 'i'
  
  
  beta0 ~ dnorm(0, 1.5)
  beta.day ~ dnorm(0, 1.5)
  beta.ice ~ dnorm(0, 1.5)		
  beta.age.y ~ dnorm(0, 1.5)	
  beta.age.o ~ dnorm(0, 1.5)
  
  ## Equations ##
  M ~ dcat(p.M[])
  
  for (i in 1:4) {
    w[i] <- values_w[M, i]
  }
  
  for (i in 1:dim.M){
    indicator_model[i] <- equals(M, i)
  }
  
})


# ~~~ b. Constants and data ----------------------------------------------------

dim.M <- 16
p.M <- rep(1/dim.M, dim.M)


# Values of w[i]
all_vectors <- c()
for (k in 1:1000) {
  vector_i <- sample(0:1, size = 4, replace = TRUE)
  all_vectors <- rbind(all_vectors, vector_i)
}
values.df <- unique(as.data.frame(all_vectors)[ , 1:4]) %>%
  mutate(sum = V1 + V2 + V3 + V4) %>%
  arrange(sum, -V1, -V2, -V3, -V4)
values_w <- as.matrix(values.df[1:4])


# Bundle data
my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     age_young = age_young,
                     age_old = age_old,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            values_w = values_w)


# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         beta.day = rnorm(1, 0, 10),
                         beta.ice = rnorm(1, 0, 10),
                         beta.age.y = rnorm(1, 0, 10),
                         beta.age.o = rnorm(1, 0, 10),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o", "M", 
            "indicator_model", "w") # "w[1]", "w[2]",) 

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_dummy_data_test_4 <-  nimbleMCMC(code = model_dummy_data_test_4,     # model code  
                                   data = dat,                                   
                                   constants = my.constants,        
                                   inits = inits,          
                                   monitors = params,   # parameters to monitor
                                   niter = 150000,                  # nb iterations
                                   nburnin = 10000,              # length of the burn-in
                                   nchains = 2,
                                   summary = TRUE,
                                   WAIC = TRUE)
end <- Sys.time()
end - start

fit_dummy_data_test_4$summary$all.chains
save(fit_dummy_data_test_4, 
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/fit_dummy_data_test_4.RData")

# ~~~ d. Check convergence ----------------------------------------------------

load("07_results/01_interim_results/model_outputs/covariate_model_selection/fit_dummy_data_test_4.RData")


params.plot <- c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o", 
                 # "indicator_model.1.", "indicator_model.2.", "indicator_model.3.", 
                 # "indicator_model.4.", "indicator_model.5.", "indicator_model.6.", 
                 # "indicator_model.7.", "indicator_model.8.", 
                 "w.1.", "w.2.", "w.3.", "w.4.", "M") 

subsample <- seq(from = 1, to = dim(fit_dummy_data_test_4[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_dummy_data_test_4[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_dummy_data_test_4[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_dummy_data_test_4[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_dummy_data_test_4[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

param.mean <- chains_l %>%
  group_by(parameter, chain) %>%
  summarise(m = mean(value))

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

# density.plots <- ggplot(data = chains_l, 
#                         aes(x = value, color = chain, fill = chain)) +
#   geom_density(alpha = 0.25) +
#   labs(x = "density") +
#   theme(legend.position = "none") +
#   facet_wrap( ~ parameter,
#               scales = "free",
#               ncol = 1)

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_histogram(alpha = 0.25,
                 position="identity") +
  labs(x = "density") +
  theme(legend.position = "none") +
  facet_wrap( ~ parameter,
              scales = "free",
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
nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_dummy_data_test_4.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)


# ~~~ e. Plot ------------------------------------------------------------------

load("07_results/01_interim_results/model_outputs/covariate_model_selection/fit_dummy_data_test_4.RData")

params.plot.2 <- c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o")

subsample <- seq(from = 1, to = dim(fit_dummy_data_test_4[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_dummy_data_test_4[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_dummy_data_test_4[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_dummy_data_test_4[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_dummy_data_test_4[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

# ~~~~~~ Caterpillar plot

caterpilar <- as.data.frame(chains) %>%
  pivot_longer(cols = c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o"))
source("05_script/models/functions_for_models_Nimble.R")

ggplot(data = caterpilar, 
       aes(x = name , y = value)) + 
  scale_x_discrete(limits = params.plot.2) +
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

ggsave(filename = "10_meetings/2021-07-02 Meeting with Olivier & Sarah/caterpilar_plot_test_9.png",
       width = 3, height = 4)



