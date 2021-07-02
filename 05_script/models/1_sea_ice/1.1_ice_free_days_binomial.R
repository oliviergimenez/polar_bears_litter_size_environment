#==============================================================================#
#                                                                              #
#                     Models with only sea ice-free days                       #
#                              (Binomial model)                                #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(viridis)
library(gridExtra)
library(nimble)
Sys.setenv(LANG = "en")

# READY VARIABLES ===========================================================
# CR data
CR_data <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/SI_metrics_D.csv")
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

# /!\ test another n
# n <- unique(n)

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

# on renumerote les individus
{id_fem <- data_model$ID_NR
  id_fem <- factor(id_fem) 
  id_fem2 <- NULL
  for (i in 1:length(id_fem)){
    id_fem2 <- c(id_fem2,which(id_fem[i]==levels(id_fem)))
  }
  id_fem <- factor(id_fem2)
  nbind <- length(levels(id_fem))
  id_fem
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/1_sea_ice/1.1_Nimble_ice_free_days.R")
source("05_script/models/functions_for_models_Nimble.R")



# RUN THE MODELS ===============================================================

# A. Null model ================================================================

# ~~~ a. Run the model ---------------------------------------------------------

model_code <- "0.0.0.0"
mode <- "_binomial"

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
assign(x = paste0("fit_", model_code, mode),
       value = nimbleMCMC(code = get(paste0("model_", model_code, mode)),     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits_null,          
                          monitors = params,   # parameters to monitor
                          niter = 20000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE))
end <- Sys.time()
end - start

get(paste0("fit_", model_code, mode))$WAIC
# 1083.937


save(list = paste0("fit_", model_code, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/", 
                   model_code, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
nimble_output <- get(paste0("fit_", model_code, mode))
# rm(list(get(paste0("fit_", model_code, "_effect_", effect))))
print(paste0("check convergence of model_", model_code, "_effect_", effect))

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






# B. Ice-free days t-1 =========================================================

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"
effect <- "_binomial"

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
}

# Define the parameters to estimate
params <- c("b0", "b1", "sigma1", "eps1")
temp.data <- data.frame(y = y, 
                        var_scaled = var_scaled)
mylogit <- glm(y ~ var_scaled, data = temp.data, family = "binomial")

inits <- function() list(b0 = mylogit$coefficients[1] + round(runif(n = 1, -1, 1))/10, 
                         b1 = mylogit$coefficients[2] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_", effect),
       value = nimbleMCMC(code = get(paste0("model_", model_code, effect)),     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits,          
                          monitors = params,   # parameters to monitor
                          thin = 10,
                          niter = 15000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE))
end <- Sys.time()
end - start


get(paste0("fit_", model_code, "_", effect))$WAIC
# 1085.991

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------

nimble_output <- get(paste0("fit_", model_code, effect))
# rm(list(get(paste0("fit_", model_code, "_effect_", effect))))
print(paste0("check convergence of model_", model_code, "_effect_", effect))

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


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
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
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode,
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"))))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)



# Plot 
range <- range(var_scaled)
lengthgrid <- 238
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(var) + mean(var)


logit_prob  <- matrix(data = NA, 
                      nrow = lengthgrid)
for (i in 1:lengthgrid) {
 logit_prob[i] <- mylogit$coefficients[1] + 
                  mylogit$coefficients[2] * grid_scaled[i]
  
}

prob <- matrix(data = NA, 
               nrow = lengthgrid)
for (i in 1:lengthgrid){
  # prob[i] <- plogis(logit_prob[i])
  prob[i] <- exp(logit_prob[i])/(1 + exp(logit_prob[i]))
}

df.for.plot <- data.frame(var = grid,
                          prob1 = prob) %>%
  mutate(prob2_3 = 1 - prob1) %>%
  pivot_longer(data = .,
               cols = c(prob1, prob2_3),
               names_to = "cub_number",
               values_to = "probability")


ggplot(df.for.plot, aes(x = var, y = probability, group = cub_number, color = cub_number)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,
                      labels = c("1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Moyenne", "IC")) +
  theme_bw() +
  labs(x = "number of ice-free days during the previous year",
       y = "proability",
       color = "")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#














# C. Ice-free days t-2 =========================================================

# ~ 1. Effect only on 1cub VS 0cubs (1.1.3_E_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_E"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"

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

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect, mode)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         a1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect, mode),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect, mode)),     # model code  
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




get(paste0("fit_", model_code, "_effect_", effect, mode))$WAIC
# 1085.51

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))

temp <- get_probabilities(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
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
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode,
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"))))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (1.1.3_E_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "1.1.3_E"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"

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

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect, mode)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         b1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect, mode),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect, mode)),     # model code  
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




get(paste0("fit_", model_code, "_effect_", effect, mode))$WAIC
# 1084.646

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))

temp <- get_probabilities(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
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
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode,
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"))))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (1.1.3_E_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_E"
effect <- "common"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"

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

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect, mode)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         a1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect, mode),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect, mode)),     # model code  
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




get(paste0("fit_", model_code, "_effect_", effect, mode))$WAIC
# 1083.797

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))

temp <- get_probabilities(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
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
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode,
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"))))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (1.1.3_E_distinct) ----------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_E"
effect <- "distinct"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"

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
                         b1 = coefs[4] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))



# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect, mode),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect, mode)),     # model code  
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




get(paste0("fit_", model_code, "_effect_", effect, mode))$WAIC
# 1085.799

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))

temp <- get_probabilities(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
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
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode,
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"))))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#















