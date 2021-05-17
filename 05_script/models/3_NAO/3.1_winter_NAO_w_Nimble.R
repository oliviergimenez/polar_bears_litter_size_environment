#==============================================================================#
#                                                                              #
#                          Models with only winter NAO                         #
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
CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
NAO_data <- read_csv("06_processed_data/NAO_data/NAO_data_winter_spring_1990-2019.csv")

data_model <- CR_data %>%
  left_join(x = CR_data,
            y = NAO_data,
            by = "year")


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

N <- length(y) # nb of reproductive events
J <- length(levels(y)) # number of categories

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear) 

# Load the JAGS models + the ancillary functions
source("05_script/models/3_NAO/3.1_Nimble_winter_NAO.R")
source("05_script/models/functions_for_models_Nimble.R")



# RUN THE MODELS ===============================================================

# A. Null model ================================================================

model_code <- "null_model"
mode <- ""

load(file = paste0("07_results/01_interim_results/model_outputs/", 
                   model_code, toupper(mode), ".RData"))

get(paste0("fit_", model_code, mode))$WAIC
# 1083.937

# B. Winter NAO t ===============================================================

# ~ 1. Effect only on 1cub VS 0cubs (3.1.1.1_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.1"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_s"
var_full_name <- "Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1115.816

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (3.1.1.1_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "3.1.1.1"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_s"
var_full_name <- "Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1116.146

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (3.1.1.1_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.1"
effect <- "common"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_s"
var_full_name <- "Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (3.1.1.1_distinct) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.1"
effect <- "distinct"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_s"
var_full_name <- "Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


# ~ 5. Effect only on 1cub VS 0cubs, FACTOR (3.1.1.2_1c_VS_0c) -----------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.2"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_pos"
var_full_name <- "Phase of the Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1117.115

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 6. Effect only on 2-3cub VS 0cubs, FACTOR (3.1.1.2_2-3c_VS_0c) -------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.2"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_pos"
var_full_name <- "Phase of the Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1116.878

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 7. Common effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.1.2_common) ------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.2"
effect <- "common"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_pos"
var_full_name <- "Phase of the Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



# ~ 8. Distinct effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.1.2_distinct) ----

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.1.2"
effect <- "distinct"

# Predictor
var <- data_model$winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_pos"
var_full_name <- "Phase of the Winter NAO"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1082.624

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


# C. Winter NAO t-1 =============================================================

# ~ 1. Effect only on 1cub VS 0cubs (3.1.2.1_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.1"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_s"
var_full_name <- "Winter NAO in previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1115.406

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (3.1.2.1_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "3.1.2.1"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_s"
var_full_name <- "Winter NAO in previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1114.781

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




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (3.1.2.1_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.1"
effect <- "common"

# Predictor
var <- data_model$prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_s"
var_full_name <- "Winter NAO in previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1115.16

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







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (3.1.2.1_distinct) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.1"
effect <- "distinct"

# Predictor
var <- data_model$prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_s"
var_full_name <- "Winter NAO in previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1115.843

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
                   paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


# ~ 5. Effect only on 1cub VS 0cubs, FACTOR (3.1.2.2_1c_VS_0c) -----------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.2"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_pos"
var_full_name <- "Phase of the Winter NAO during previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1116.25

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 6. Effect only on 2-3cub VS 0cubs, FACTOR (3.1.2.2_2-3c_VS_0c) -------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.2"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_pos"
var_full_name <- "Phase of the Winter NAO during previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1115.799

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 7. Common effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.2.2_common) ------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.2"
effect <- "common"

# Predictor
var <- data_model$prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_pos"
var_full_name <- "Phase of the Winter NAO during previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 1115.9

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 8. Distinct effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.2.2_distinct) ----

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.2.2"
effect <- "distinct"

# Predictor
var <- data_model$prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_prior_pos"
var_full_name <- "Phase of the Winter NAO during previous year"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1116.915

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)





# D. Winter NAO t-2 =========================================================

# ~ 1. Effect only on 1cub VS 0cubs (3.1.3.1_1c_VS_0c) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.1"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$two_year_prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_s"
var_full_name <- "Winter NAO two years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1116.826

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


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 2. Effect only on 2-3cub VS 0cubs (3.1.3.1_2-3c_VS_0c) ---------------------

# ~~~ a. Run the model ---------------------------------------------------------


{model_code <- "3.1.3.1"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$two_year_prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_s"
var_full_name <- "Winter NAO two years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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




# ~ 3. Common effect of 1c VS 0 and of 2-3c VS 0 (3.1.3.1_common) --------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.1"
effect <- "common"

# Predictor
var <- data_model$two_year_prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_s"
var_full_name <- "Winter NAO two years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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







# ~ 4. Distinct effect of 1c VS 0 and of 2-3c VS 0 (3.1.3.1_distinct) ----------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.1"
effect <- "distinct"

# Predictor
var <- data_model$two_year_prior_winter_NAO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_s"
var_full_name <- "Winter NAO two years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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




# ~ 5. Effect only on 1cub VS 0cubs, FACTOR (3.1.1.2_1c_VS_0c) -----------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.2"
effect <- "1c_VS_0c"

# Predictor
var <- data_model$two_year_prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_pos"
var_full_name <- "Phase of the Winter NAO 2 years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1116.179


save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 6. Effect only on 2-3cub VS 0cubs, FACTOR (3.1.1.2_2-3c_VS_0c) -------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.2"
effect <- "2_3c_VS_0c"

# Predictor
var <- data_model$two_year_prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_pos"
var_full_name <- "Phase of the Winter NAO 2 years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#




# ~ 7. Common effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.1.2_common) ------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.2"
effect <- "common"

# Predictor
var <- data_model$two_year_prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_pos"
var_full_name <- "Phase of the Winter NAO 2 years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

# Define the parameters to estimate
params <- get_coefs_and_params(y, var_scaled, effect, mode)$params
# params <- c("a0", "b0", "b1", "sigma1", "eps1") 

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
# 

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#







# ~ 8. Distinct effect of 1c VS 0 and of 2-3c VS 0, FACTOR (3.1.1.2_distinct) ----

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "3.1.3.2"
effect <- "distinct"

# Predictor
var <- data_model$two_year_prior_winter_NAO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_NAO_2y_prior_pos"
var_full_name <- "Phase of the Winter NAO 2 years before"

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

dat <- list(as.numeric(y), as.numeric(var_scaled))
names(dat) <- c("y", var_short_name)

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
# 1116.835

save(list = paste0("fit_", model_code, "_effect_", effect, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, toupper(mode), ".RData"))


temp <- get_probabilities_factor(model_code, effect, mode, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")), 
       aes(x = as.factor(var), y = probability, fill = nbr_cub)) +
  geom_violin() +
  theme_bw() +
  scale_x_discrete(labels = c("negative", "positive")) +
  scale_fill_viridis(discrete = TRUE, 
                     limits = color_labels,
                     labels = color_labels) +
  theme(legend.position = "none") +
  labs(x = var_full_name,
       y = "Probability") +
  facet_wrap( ~ nbr_cub)

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, mode, ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_effect_", effect, mode),
            paste0("fit_", model_code, "_effect_", effect, mode, "_for_plot")))
rm(model_code, effect, dat, params, coefs, inits, 
   var, var_scaled, var_short_name, var_full_name,
   color_labels)











