#==============================================================================#
#                                                                              #
#                         Models with only ice-free days                       #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(R2jags)
library(viridis)
library(ggmcmc)
library(gridExtra)
Sys.setenv(LANG = "en")

# A. READY VARIABLES ===========================================================
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
  filter(year != 1992)  %>% # Remove the captures from 1992 since I can't calculate the ice free days in 1991
  filter(year != 1993) # Same 


# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number))
summary(y)

# Renumérotation des années
year <- data_model$year
year <- factor(year) # laisse tomber les modalites inutiles 
year2 <- NULL
for (i in 1:length(year)){
  year2 <- c(year2, which(year[i] == levels(year)))
}
year <- factor(year2)
nbyear <- length(levels(year))
year

# on renumerote les individus
id_fem <- data_model$ID_NR
id_fem <- factor(id_fem) 
id_fem2 <- NULL
for (i in 1:length(id_fem)){
  id_fem2 <- c(id_fem2,which(id_fem[i]==levels(id_fem)))
}
id_fem <- factor(id_fem2)
nbind <- length(levels(id_fem))

id_fem

N <- length(y) # nb of reproductive events
J <- length(levels(y)) # number of categories


# Load the JAGS models + the ancillary functions
source("05_script/models/1_sea_ice/1.1_JAGS_ice_free_days.R")
source("05_script/models/functions_for_models.R")





# B. RUN THE MODELS ============================================================

# 1. Previous year ice-free days: common slope (1.1.2_D_common) ------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"

# Number of slopes
slope <- "common"
# slope <- "distinct" 

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var_short_name)

# Define the parameters to estimate
params <- define_params(mode, slope)[[1]]
nb.beta <- define_params(mode, slope)[[2]]

# Generate starting values
inits <- get_inits(y, var_scaled, mode, slope)


# Run the model
assign(x = paste0("fit_", model_code, "_", slope, mode),
       value = jags(data = dat, 
                    inits = inits, 
                    parameters.to.save = params, 
                    n.chains = 2, 
                    n.iter = 30000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, "_slope", mode))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode),
#        value = get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


temp <- get_probabilities(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode, "_for_plot")), 
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

ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", 
                         model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




# 2. Previous year ice-free days: distinct slopes (1.1.2_D_distinct) -----------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "1.1.2_D"

# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"

# Number of slopes
slope <- "distinct"
# slope <- "distinct" 

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var_short_name)

# Define the parameters to estimate
params <- define_params(mode, slope)[[1]]
nb.beta <- define_params(mode, slope)[[2]]

# Generate starting values
inits <- get_inits(y, var_scaled, mode, slope)


# Run the model
assign(x = paste0("fit_", model_code, "_", slope, mode),
       value = jags(data = dat, 
                    inits = inits, 
                    parameters.to.save = params, 
                    n.chains = 2, 
                    n.iter = 30000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, "_slope", mode))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode),
#        value = get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)




# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


temp <- get_probabilities(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode, "_for_plot")), 
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

ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", 
                         model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 3. 2y prior ice-free days: common slope (1.1.3_D_common) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <- "1.1.3_D"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days 2y prior to capture"

# Number of slopes
slope <- "common"
# slope <- "distinct" 

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var_short_name)

# Define the parameters to estimate
params <- define_params(mode, slope)[[1]]
nb.beta <- define_params(mode, slope)[[2]]

# Generate starting values
inits <- get_inits(y, var_scaled, mode, slope)


# Run the model
assign(x = paste0("fit_", model_code, "_", slope, mode),
       value = jags(data = dat, 
                    inits = inits, 
                    parameters.to.save = params, 
                    n.chains = 2, 
                    n.iter = 30000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, "_slope", mode))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode),
#        value = get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


temp <- get_probabilities(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode, "_for_plot")), 
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

ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 4. 2y prior ice-free days: distinct slopes (1.1.3_D_distinct) ====================

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <-  "1.1.3_D"

# Predictor
var <- data_model$ice_free_days_2y_prior
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_2y_prior_s"
var_full_name <- "Ice-free days 2y prior to capture"

# Number of slopes
slope <- "distinct"
# slope <- "distinct" 

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var_short_name)

# Define the parameters to estimate
params <- define_params(mode, slope)[[1]]
nb.beta <- define_params(mode, slope)[[2]]

# Generate starting values
inits <- get_inits(y, var_scaled, mode, slope)


# Run the model
assign(x = paste0("fit_", model_code, "_", slope, mode),
       value = jags(data = dat, 
                    inits = inits, 
                    parameters.to.save = params, 
                    n.chains = 2, 
                    n.iter = 30000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, "_slope", mode))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode),
#        value = get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))



# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)




# ~~~ c. Plot the model --------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

temp <- get_probabilities(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode, "_for_plot")), 
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

ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode, "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 5. Nimble test ---------------------

library(nimble)
Nimble_ie_free_days_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
})

{model_code <- "1.1.2_D_Nimble"
  
# Predictor
var <- data_model$ice_free_days_previous
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "ice_free_days_previous_s"
var_full_name <- "Ice-free days in previous year"
  
# Number of slopes
slope <- "common"
# slope <- "distinct" 
  
# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var_short_name)

# Define the parameters to estimate
params <- define_params(mode, slope)[[1]]

# Generate starting values
inits <- get_inits(y, var_scaled, mode, slope)

temp.dat <- data.frame(y = y, var_scaled = var_scaled)
temp.dat$yfac <- as.factor(temp.dat$y)   # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| var_scaled, 
                     data = mnl.dat, 
                     reflevel = ifelse(mode == "_bis", "1", "0"))

coefs <- as.vector(summary(mlogit.mod)$coefficients)
inits1 <- list(a0 = coefs[1], b0 = coefs[2], c0 = coefs[3],
               a1 = coefs[4], eps1 = 0,
               sigma1 = runif(1))
inits2 <- list(a0 = coefs[1] + 0.1, b0 = coefs[2] - 0.1, c0 = coefs[3] + 0.1,
               a1 = coefs[4] -0.1,
               sigma1 = runif(1), eps1 = 0)
  
inits <- list(inits1, inits2)


my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y))) 

mcmc.output <- nimbleMCMC(code = Nimble_ie_free_days_common,     # model code  
                          data = dat,                  # data
                          constants = my.constants,        # constants
                          inits = get_inits(y, var_scaled, mode, slope),          # initial values
                          monitors = params,   # parameters to monitor
                          niter = 30000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 1)
initial.values <- function() {
  list(phi = runif(1,0,1))
}
