#==============================================================================#
#                                                                              #
#                         Models with only ice-free days                       #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(mlogit)
library(viridis)
library(gridExtra)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")

# READY VARIABLES ===========================================================
# CR data
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
  filter(year >= 2000)


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

{model_code <- "0.0.0.0"

dat <- list(y = as.numeric(y))
params <- c("a0", "b0", "sigma1", "eps1") 

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear) 
}

temp.dat <- data.frame(y = y)
temp.dat$yfac <- as.factor(temp.dat$y)   # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)

inits_null <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                              b0 = coefs[2] + round(runif(n = 1, -1, 1))/10)


start <- Sys.time()
assign(x = paste0("fit_", model_code),
       value = nimbleMCMC(code = get(paste0("model_", model_code)),     # model code  
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

# Get and save wAIC
wAIC_table <- save_wAIC_null_model(model_code)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code), 
     file = paste0("07_results/01_interim_results/model_outputs/", 
                   model_code, ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/", 
                   model_code, ".RData"))
nimble_output <- get(paste0("fit_", model_code))

print(paste0("check convergence of model_", model_code))

# Process Nimble output into dataframe
chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
  select(params[-length(params)]) %>%
  mutate(chain = "chain 1",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
  select(params[-length(params)]) %>%
  mutate(chain = "chain 2",
         iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
chains <- rbind(chain1, chain2)


trace_a0 <- ggplot(data = chains, aes(x = iteration, y = a0, color = chain)) +
  geom_line() +
  labs(y = "a0")
density_a0 <- ggplot(data = chains, 
                     aes(x = a0, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
  labs(x = "a0") +
  theme(legend.position = "none")

trace_b0 <- ggplot(data = chains, aes(x = iteration, y = b0, color = chain)) +
  geom_line() +
  labs(y = "b0")
density_b0 <- ggplot(data = chains, 
                     aes(x = b0, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
  labs(x = "b0") +
  theme(legend.position = "none")

trace_sigma1 <- ggplot(data = chains, aes(x = iteration, y = sigma1, color = chain)) +
  geom_line() +
  labs(y = "sigma1")
density_sigma1 <- ggplot(data = chains, 
                         aes(x = sigma1, color = chain, fill = chain)) +
  geom_density(alpha = 0.25) +
  labs(x = "sigma1") +
  theme(legend.position = "none")

diagnostic_plot <- plot_grid(trace_a0, density_a0,
                             trace_b0, density_b0,
                             trace_sigma1, density_sigma1,
                             ncol = 2, nrow = 3)

save_plot(filename = paste0("07_results/01_interim_results/model_outputs/graph/diagnostic_plots/fit_",
                            model_code, ".png"), 
          plot = diagnostic_plot,
          ncol = 2,
          nrow = 3)




# B. Ice-free days t-1 =========================================================

# ~ 1. Effect only on 1cub VS 0cubs (1.1.2_D_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"


my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs
  
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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))


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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (1.1.2_D_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "1.1.2_D"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"




my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}



# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         b1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)

# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (1.1.2_D_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"
effect <- "common"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (1.1.2_D_distinct) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"
effect <- "distinct"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         a1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         b1 = coefs[4] + round(runif(n = 1, -1, 1))/10,
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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






# C. Ice-free days t-2 =========================================================

# ~ 1. Effect only on 1cub VS 0cubs (1.1.3_D_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_D"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (1.1.3_D_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "1.1.3_D"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         b1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (1.1.3_D_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_D"
effect <- "common"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (1.1.3_D_distinct) ----------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.3_D"
effect <- "distinct"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days two years before"





my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[5] <- var_short_name
}


# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

# Generate starting values
coefs <- get_coefs_and_params(y, var_scaled, effect)$coefs

inits <- function() list(a0 = coefs[1] + round(runif(n = 1, -1, 1))/10, 
                         b0 = coefs[2] + round(runif(n = 1, -1, 1))/10, 
                         a1 = coefs[3] + round(runif(n = 1, -1, 1))/10,
                         b1 = coefs[4] + round(runif(n = 1, -1, 1))/10,
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


# Get and save wAIC
wAIC_table <- save_wAIC(model_code, var_short_name, effect)
write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table.csv")

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#















