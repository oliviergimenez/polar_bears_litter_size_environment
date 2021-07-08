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


# READY VARIABLES ==============================================================
# ~~~ a. Load the data ---------------------------------------------------------
# CR data
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

# Load the Nimble models + the ancillary functions
source("05_script/models/functions_for_models.R")



# RUN THE MODELS ============================================================

# ~ 1. Ice-free days ---------------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:2) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
  M ~ dcat(p.M[])
})


# ~~~ b. Constants and data ----------------------------------------------------

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s),
                     p.M = c(0.5, 0.5))

dat <- list(y = as.numeric(y))

params <- c("b", "sigma1" ,"eps1", "M")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       ice_free_days_previous_s = ice_free_days_previous_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1),
                         M = sample(0:1, 1))


# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test <- nimbleMCMC(code = model_1.1.2_D_test,     # model code  
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

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("b.1.", "b.2.", "sigma1", "M")


chain1 <- data.frame(fit_1.1.2_D_test[["samples"]][["chain1"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test[["samples"]][["chain2"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test[["samples"]][["chain2"]])[1], by = 1))
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

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
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


# ~ 2. Ice-free days with indicator ---------------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test_2 <- nimbleCode({
 
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b[1] + 
      w2 * b[2] * ice_free_days_previous_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:dim.M) {
    b[i] ~ dnorm(0.00000E+00, sd = V[i])
  }
  
  v ~ dunif(0, 100)
  for (i in 1:dim.M) {
    V[i] <- v*nbr_param_weight[M, i]
  }
  
  M ~ dcat(p.M[])
  
  w2 <- equals(M, as.integer(2))
  
  # for (i in 1:dim.M){
  #   indicator_model[i] <- equals(M,i)
  # }
})


# ~~~ b. Constants and data ----------------------------------------------------

nbr_param_weight <- matrix(data = c(1, 0.5, 0, 0.5),
                     nrow = 2)

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s),
                     dim.M = dim.M)

dat <- list(y = as.numeric(y),
            nbr_param_weight = nbr_param_weight,
            p.M = p.M)

# ~~~ c. Params and inits

params <- c("b", "sigma1" ,"eps1", "M", "v", "w2") #, "indicator_model")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       ice_free_days_previous_s = ice_free_days_previous_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1),
                         v = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))




# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_2 <- nimbleMCMC(code = model_1.1.2_D_test_2,     # model code  
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


# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("b.1.", "b.2.", "sigma1" , "M", "v", "w2") 
                 # "indicator_model.1.", "indicator_model.2.")


chain1 <- data.frame(fit_1.1.2_D_test_2[["samples"]][["chain1"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_2[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_2[["samples"]][["chain2"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_2[["samples"]][["chain2"]])[1], by = 1))
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

density.plots <- ggplot(data = chains_l, 
                        aes(x = value, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
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


# ~ 3. Ice-free days with indicator ---------------------

# ~~~ a. Build the model ------------------------------------------------------

# model_1.1.2_D_test_3 <- nimbleCode({
#   for(i in 1:N) {
#     # y[i] ~ dbin(p[i], 1)
#     y[i] ~ dbern(p[i])
#     logit(p[i]) <- b[1] + 
#       b[2] * ice_free_days_previous_s[i] +
#       eps1[year[i]]
#   }
#   for (i in 1:nbyear) {
#     eps1[i] ~ dnorm(0, sd = sigma1)
#   }
#   sigma1 ~ dunif(0.00000E+00, 10)
#   
#   for (i in 1:2) {
#     b[i] ~ dnorm(0.00000E+00, sd = 1.5)
#   }
#   M ~ dcat(p.M[])
#   w2 <- M - 1
#   v ~ dunif(1, 100)
# })

model_1.1.2_D_test_3 <- nimbleCode({
  for(i in 1:N) {
    # y[i] ~ dbin(p[i], 1)
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b[1] + 
      w2 * b[2] * ice_free_days_previous_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b[1] ~ dnorm(0.00000E+00, v1)
  b[2] ~ dnorm(0.00000E+00, v2)

  M ~ dcat(p.M[])
  w2 <- M - 1
  
  v ~ dunif(1, 100)
  v1 <- v/M
  v2 <- v/M
  
})

# ~~~ b. Constants and data ----------------------------------------------------

nbr_param_weight <- matrix(data = c(1, 2, 1, 2),
                     nrow = 2)

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)

my.constants <- list(N = length(y), # nb of females captured
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s))
                     # dim.M = dim.M)

dat <- list(y = as.numeric(y),
            # nbr_param_weight = nbr_param_weight,
            p.M = p.M)

params <- c("b", "sigma1" ,"eps1", "M", "w2", "v")


# Generate initial values
temp.dat <- data.frame(y = y, 
                       ice_free_days_previous_s = ice_free_days_previous_s)

mylogit <- glm(as.factor(y) ~ ice_free_days_previous_s, 
               data = temp.dat, 
               family = "binomial")

coefs <- as.vector(mylogit$coefficients)

inits <- function() list(b = coefs[c(1, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1),
                         M = sample(0:1, 1),
                         v = runif(1, 1, 100))

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test <- nimbleMCMC(code = model_1.1.2_D_test_3,     # model code  
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

fit_1.1.2_D_test$summary$all.chains

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("b.1.", "b.2.", "sigma1", "M", "w2", "v")


chain1 <- data.frame(fit_1.1.2_D_test[["samples"]][["chain1"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test[["samples"]][["chain2"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test[["samples"]][["chain2"]])[1], by = 1))
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




# ~ 4. Ice-free days with indicator, starting from Thierry's code --------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test_4 <- nimbleCode({

  for (i in 1:nind){
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


# ~~~ b. Constants and data ----------------------------------------------------
data <- data_model

# Number of rows
nind <- dim(data)[1]


#-------------------------------------------------------------------------------#
### COVARIATE: Number of ice-free days t-1 

ice_free_days_previous <- data$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

#-------------------------------------------------------------------------------#
### COVARIATE: date of capture

day_number <- data$day_number
day_number_s <- as.vector(scale(day_number))


### ********************************************************** ###
y <- data$success
head(y) ; length(y)


### ********************************************************** ###
# Number of Parameter weight for each Model
nbr_param_weight <- as.matrix(read.csv("03_methods_protocols/Model selection Bayesian Kuo & Mallick/Model-Recr-v2/Model-Matrix.csv", sep = ";", header = F))
nbr_param_weight

dim.M <- dim(nbr_param_weight)[1] 
p.M <- rep(1/dim.M, dim.M)


### ********************************************************** ###
# Bundle data
my.constants <- list(nind = nind,
                     day_number_s = day_number_s,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            nbr_param_weight = nbr_param_weight)

### **********************************************************

inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         beta.day = rnorm(1, 0, 10),
                         beta.ice = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))

params <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model", "w.ice",
            "M") 



# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_4 <-  nimbleMCMC(code = model_1.1.2_D_test_4,     # model code  
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


# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model.1.",
                 "indicator_model.2.", "w.ice", "M") 

chain1 <- data.frame(fit_1.1.2_D_test_4[["samples"]][["chain1"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_4[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_4[["samples"]][["chain2"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_4[["samples"]][["chain2"]])[1], by = 1))
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_4.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

# ~ 5. Ice-free days with indicator, adapted ot my style -----------------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test_5 <- nimbleCode({
  
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


# ~~~ b. Constants and data ----------------------------------------------------

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# y
y <- data_model$success
head(y) ; length(y)

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

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_5 <-  nimbleMCMC(code = model_1.1.2_D_test_5,     # model code  
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


# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model.1.",
                 "indicator_model.2.", "w.ice", "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_5[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_5[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_5[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_5[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_5[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_5.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)



# ~ 6. Ice-free days with 2 indicators -----------------------------------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test_6 <- nimbleCode({
  
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
  
  # for (i in 1:dim.M){
  #   indicator_model[i] <- values_indicator[M, i]
  # }
})

# equals(M,i)
##### i = 1
# M = 1 => indicator_model[1] <- 1
# M = 2 => indicator_model[1] <- 0
##### i = 2
# M = 1 => indicator_model[2] <- 0
# M = 2 => indicator_model[2] <- 1

# ~~~ b. Constants and data ----------------------------------------------------

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# y
y <- data_model$success
head(y) ; length(y)

# Number of Parameter weight for each Model
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

# Values of indicator_model[i]
values_indicator <- as.matrix(rbind(c(1, 0, 0, 0),
                                    c(0, 1, 0, 0),
                                    c(0, 0, 1, 0),
                                    c(0, 0, 0, 1)))

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
fit_1.1.2_D_test_6 <-  nimbleMCMC(code = model_1.1.2_D_test_6,     # model code  
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

fit_1.1.2_D_test_6$summary$all.chains

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model.1.",
                 "indicator_model.2.", "indicator_model.3.", "indicator_model.4.", 
                 "w.1.", "w.2.", "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_6[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_6[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_6[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_6[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_6[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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





# ~ 7. Ice-free days with 2 indicators + random effect -------------------------

# ~~~ a. Build the model -------------------------------------------------------

model_1.1.2_D_test_7 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * beta.day * day_number_s[i] +
      w[2] * beta.ice * ice_free_days_previous_s[i] +
      eps1[year[i]]
      
  } # 'i'
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
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

# ~~~ b. Constants and data ----------------------------------------------------

# year
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

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# y
y <- data_model$success
head(y) ; length(y)

# Number of Parameter weight for each Model
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

# Values of indicator_model[i]
values_indicator <- as.matrix(rbind(c(1, 0, 0, 0),
                                    c(0, 1, 0, 0),
                                    c(0, 0, 1, 0),
                                    c(0, 0, 0, 1)))

# Bundle data
my.constants <- list(N = length(y),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
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
                         M = sample(dim.M, size = 1),
                         sigma1 = runif(1))


# Parameters to monitor
params <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model", "w[1]", "w[2]",
            "M", "sigma1", "eps1") 

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_7 <-  nimbleMCMC(code = model_1.1.2_D_test_7,     # model code  
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

fit_1.1.2_D_test_7$summary$all.chains

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "V", "indicator_model.1.",
                 "indicator_model.2.", "indicator_model.3.", "indicator_model.4.", 
                 "w.1.", "w.2.", "M", "sigma1") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_7[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_7[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_7[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_7[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_7[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
# diagnostic_plot

nrows = length(params.plot)
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_7.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)





# ~ 8. Ice-free days with 3 indicators -----------------------------------------

# ~~~ a. Build the model ------------------------------------------------------

model_1.1.2_D_test_8 <- nimbleCode({
  
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

# equals(M,i)
##### i = 1
# M = 1 => indicator_model[1] <- 1
# M = 2 => indicator_model[1] <- 0
##### i = 2
# M = 1 => indicator_model[2] <- 0
# M = 2 => indicator_model[2] <- 1

# ~~~ b. Constants and data ----------------------------------------------------

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# Size
size <- data_model$s_length
size_s <- as.vector(scale(size)) 

# y
y <- data_model$success
head(y) ; length(y)

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
#values_indicator = values_indicator)


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

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_8 <-  nimbleMCMC(code = model_1.1.2_D_test_8,     # model code  
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

fit_1.1.2_D_test_8$summary$all.chains

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "beta.size", "V", "indicator_model.1.",
                 "indicator_model.2.", "indicator_model.3.", "indicator_model.4.", 
                 "indicator_model.5.", "indicator_model.6.", "indicator_model.7.", 
                 "indicator_model.8.", "w.1.", "w.2.", "w.3.", "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_8[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_8[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_8[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_8[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_8[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_8.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)




# ~ 10. With 4 indicators without cnst variance --------------------------------

# ~~~ a. Build the model -------------------------------------------------------

model_1.1.2_D_test_10 <- nimbleCode({
  
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

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))


# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))


# Age (in 3 categories)
{age <- as.numeric(data_model$age_for_analyses)
  age_factor <- as.factor(ifelse (age < 9, "b. young", 
                                  ifelse(age > 15, "c. old", "a. prime_age")))
  relevel(age_factor, ref = "a. prime_age")
  
  age_old <- ifelse (age > 15, 1, 0)
  age_young <- ifelse (age < 9, 1, 0)
}


# y
y <- data_model$success
head(y) ; length(y)

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
fit_1.1.2_D_test_10 <-  nimbleMCMC(code = model_1.1.2_D_test_10,     # model code  
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

fit_1.1.2_D_test_10$summary$all.chains
save(fit_1.1.2_D_test_10, 
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/fit_1.1.2_D_test_10.RData")

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o", 
                 # "indicator_model.1.", "indicator_model.2.", "indicator_model.3.", 
                 # "indicator_model.4.", "indicator_model.5.", "indicator_model.6.", 
                 # "indicator_model.7.", "indicator_model.8.", 
                 "w.1.", "w.2.", "w.3.", "w.4.", "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_10[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_10[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_10[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_10[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_10[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_10.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)


# ~~~ e. Plot ------------------------------------------------------------------

params.plot.2 <- c("beta0", "beta.day", "beta.ice", "beta.age.y", "beta.age.o")

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_10[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_10[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_10[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_10[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_10[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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


# ~ 11. With indicators for all indiv variables without cnst variance ----------

# ~~~ a. Build the model -------------------------------------------------------

model_1.1.2_D_test_11 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * b[1] * day_number_s[i] +
      w[2] * b[2] * age_young[i] +
      w[3] * b[3] * age_old[i] +
      w[4] * b[4] * size_s[i] +
      w[5] * b[5] * size_s2[i] +
      w[6] * b[6] * day_number_s[i] * age_young[i] +
      w[7] * b[7] * day_number_s[i] * age_old[i]
  } # 'i'
  
  
  beta0 ~ dnorm(0, 1.5)
  
  M ~ dcat(p.M[])
  for (i in 1:7) {
    b[i] ~ dnorm(0, 1.5)
    w[i] <- values_w[M, i]
  }
  
  ## Equations ##
  
})


# ~~~ b. Constants and data ----------------------------------------------------

nbr_indicators <- 7

# Date of capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))


# Age (in 3 categories)
{age <- as.numeric(data_model$age_for_analyses)
  age_factor <- as.factor(ifelse (age < 9, "b. young", 
                                  ifelse(age > 15, "c. old", "a. prime_age")))
  relevel(age_factor, ref = "a. prime_age")
  
  age_old <- ifelse (age > 15, 1, 0)
  age_young <- ifelse (age < 9, 1, 0)
}

# Size
size <- data_model$s_length
size_s <- as.vector(scale(size))

# Size squarred
size_s2 <- size_s^2

# Number of ice-free days t-1 
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# y
y <- data_model$success
head(y) ; length(y)

dim.M <- 2^nbr_indicators
p.M <- rep(1/dim.M, dim.M)


# Values of w[i]
all_vectors <- c()
for (k in 1:5000) {
  vector_i <- sample(0:1, size = nbr_indicators, replace = TRUE)
  all_vectors <- rbind(all_vectors, vector_i)
}
values.df <- unique(as.data.frame(all_vectors)[ , 1:nbr_indicators]) %>%
  mutate(sum = V1 + V2 + V3 + V4 + V5 + V6 + V7) %>%
  arrange(sum, -V1, -V2, -V3, -V4, -V5, -V6, -V7)
values_w <- as.matrix(values.df[1:nbr_indicators])


# Bundle data
my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     age_young = age_young,
                     age_old = age_old,
                     size_s = size_s,
                     size_s2 = size_s2,
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            values_w = values_w)


# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         b = rnorm(7, 0, 10),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "b", "w", "M") 

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_11 <-  nimbleMCMC(code = model_1.1.2_D_test_11,     # model code  
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

fit_1.1.2_D_test_11$summary$all.chains
save(fit_1.1.2_D_test_11, 
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/fit_1.1.2_D_test_11.RData")

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0",
                 "b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.",
                 "w.1.", "w.2.", "w.3.", "w.4.", "w.5.", "w.6.", "w.7.", 
                 "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_11[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_11[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_11[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_11[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_11[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_11.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)



# ~~~ e. Plot ------------------------------------------------------------------

params.plot.2 <- c("beta0",
                  "b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7.") 
subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_11[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_11[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_11[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_11[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_11[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

# ~~~~~~ Caterpillar plot

caterpilar <- as.data.frame(chains) %>%
  pivot_longer(cols = c("beta0",
                        "b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", "b.7."))
source("05_script/models/functions_for_models.R")

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

ggsave(filename = "07_results/01_interim_results/model_outputs/covariate_model_selection/graphs/caterpilar_plot_test_11.png",
       width = 3, height = 4)



# ~ 12. With 10 indicators without cnst variance --------------------------------

# ~~~ a. Build the model -------------------------------------------------------

model_1.1.2_D_test_12 <- nimbleCode({
  
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    
    logit(p[i]) <- beta0 +
      w[1] * b[1] * day_number_s[i] +
      w[2] * b[2] * age_young[i] +
      w[3] * b[3] * age_old[i] +
      w[4] * b[4] * size_s[i] +
      w[5] * b[5] * size_s2[i] +
      w[6] * b[6] * day_number_s[i] * age_young[i] +
      w[7] * b[7] * day_number_s[i] * age_old[i] +
      w[8] * b[8] * ice_free_days_previous_s[i] +
      w[9] * b[9] * winter_AO_s[i] +
      w[10] * b[10] * prior_spring_AO_s[i]
      
  } # 'i'
  beta0 ~ dnorm(0, 1.5)
  
  M ~ dcat(p.M[])
  
  for (i in 1:10) {
    b[i] ~ dnorm(0, 1.5)
    w[i] <- values_w[M, i]
  } 
  
})


# ~~~ b. Constants and data ----------------------------------------------------

# Date of capture
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

# Size
size <- data_model$s_length
size_s <- as.vector(scale(size))

# Size squarred
size_s2 <- size_s^2

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Winter AO
winter_AO <- data_model$winter_AO
winter_AO_s <- as.vector(scale(winter_AO))

# Spring AO t-1 
prior_spring_AO <- data_model$prior_spring_AO
prior_spring_AO_s <- as.vector(scale(prior_spring_AO))

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


# y
y <- data_model %>%
  pull(success)

dim.M <- 1024
p.M <- rep(1/dim.M, dim.M)


# Values of w[i]
all_vectors <- c()
for (k in 1:100000) {
  vector_i <- sample(0:1, size = 10, replace = TRUE)
  all_vectors <- rbind(all_vectors, vector_i)
}
values.df <- unique(as.data.frame(all_vectors)[ , 1:10]) %>%
  mutate(sum = V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10) %>%
  arrange(sum, -V1, -V2, -V3, -V4, -V5, -V6, -V7, -V8, -V9,  -V10)
values_w <- as.matrix(values.df[1:10])


# Bundle data
# my.constants <- list(N = length(y),
#                      day_number_s = day_number_s,
#                      ice_free_days_previous_s = ice_free_days_previous_s,
#                      age_young = age_young,
#                      age_old = age_old,
#                      dim.M = dim.M,  
#                      p.M = p.M)
# 
# dat <- list(y = as.numeric(y),
#             values_w = values_w)

my.constants <- list(N = length(y),
                     day_number_s = day_number_s,
                     age_young = age_young,
                     age_old = age_old,
                     size_s = size_s,
                     size_s2 = size_s2,
                     ice_free_days_previous_s = ice_free_days_previous_s,
                     winter_AO_s = winter_AO_s,
                     prior_spring_AO_s = prior_spring_AO_s,
                     dim.M = dim.M,  
                     p.M = p.M)

dat  <- list(y = y, 
             values_w = values_w)

# Initial values
inits <- function() list(beta0 = rnorm(1, 0, 10), 
                         b = rnorm(10, 0, 10),
                         M = sample(dim.M, size = 1))


# Parameters to monitor
params <- c("beta0", "b", "M", "w") # "w[1]", "w[2]",) 

# ~~~ c. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_1.1.2_D_test_12 <-  nimbleMCMC(code = model_1.1.2_D_test_12,     # model code  
                                   data = dat,                                   
                                   constants = my.constants,        
                                   inits = inits,          
                                   monitors = params,   # parameters to monitor
                                   niter = 200000,                  # nb iterations
                                   nburnin = 10000,              # length of the burn-in
                                   nchains = 2,
                                   summary = TRUE,
                                   WAIC = TRUE)
end <- Sys.time()
end - start

fit_1.1.2_D_test_12$summary$all.chains
save(fit_1.1.2_D_test_12, 
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/fit_1.1.2_D_test_12.RData")

# ~~~ d. Check convergence ----------------------------------------------------

params.plot <- c("beta0",
                 "b.1.", "b.2.", "b.3.", "b.4.", "b.5.", "b.6.", 
                 "b.7.", "b.8.", "b.9.", "b.10.",
                 "w.1.", "w.2.", "w.3.", "w.4.", "w.5.", "w.6.", 
                 "w.7.", "w.8.", "w.9.", "w.10.", 
                 "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_12[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_12[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_12[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_12[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_12[["samples"]][["chain2"]][subsample, ])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, 
                         cols = all_of(params.plot), 
                         names_to = "parameter") 

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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/fit_1.1.2_D_test_12.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)


# ~~~ e. Plot ------------------------------------------------------------------

params.plot.2 <- c("beta0",
                   "b.1.", "b.2.", "b.3.", "b.4.", "b.5.", 
                   "b.6.", "b.7.", "b.8.", "b.9.", "b.10.",
                    "M") 

subsample <- seq(from = 1, to = dim(fit_1.1.2_D_test_12[["samples"]][["chain1"]])[1], by = 5)


chain1 <- data.frame(fit_1.1.2_D_test_12[["samples"]][["chain1"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_1.1.2_D_test_12[["samples"]][["chain1"]][subsample, ])[1], by = 1))
chain2 <- data.frame(fit_1.1.2_D_test_12[["samples"]][["chain2"]][subsample, ]) %>%
  dplyr::select(all_of(params.plot.2)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_1.1.2_D_test_12[["samples"]][["chain2"]][subsample, ])[1], by = 1))
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
