#==============================================================================#
#                                                                              #
#                       Models with only date of capture                       #
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

data_model <- read_csv("06_processed_data/CR_data/CR_f_clean.csv") %>%
  filter(year >= 2000)


# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- (day - mean(day))/sd(day)

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
source("05_script/models/4_day_capture/4_Nimble_day_capture.R")
source("05_script/models/functions_for_models_Nimble.R")

# A. Null model ================================================================

# see in 1.1_ice_free_days_w_Nimble.R


# B. Day of capture ============================================================

# ~ 1. Effect only on 1cub VS 0cubs (4_1c_VS_0c) -------------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "4"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$day_number
var_scaled <- scale(var)
var_short_name <- "day_s"
var_full_name <- "Day of capture"


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


temp <- get_probabilities(model_code, effect, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, "_for_plot")), 
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
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
            paste0("fit_", model_code, "_effect_", effect, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (4_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "4"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$day_number
var_scaled <- scale(var)
var_short_name <- "day_s"
var_full_name <- "Day of capture"



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


temp <- get_probabilities(model_code, effect, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, "_for_plot")), 
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
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
            paste0("fit_", model_code, "_effect_", effect, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (4_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "4"
effect <- "common"

# Predictor
var <- data_model$day_number
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "day_s"
var_full_name <- "Day of capture"


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
                          niter = 25000,                  # nb iterations
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


temp <- get_probabilities(model_code, effect, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, "_for_plot")), 
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
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
            paste0("fit_", model_code, "_effect_", effect, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (4_distinct) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "4"
effect <- "distinct"

# Predictor
var <- data_model$day_number
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "day_s"
var_full_name <- "Day of capture"





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
                          niter = 25000,                  # nb iterations
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


temp <- get_probabilities(model_code, effect, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, "_for_plot")), 
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
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect),
            paste0("fit_", model_code, "_effect_", effect, "_for_plot")))
rm(model_code, effect, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#