#==============================================================================#
#                                                                              #
#                         Models with only winter AO                           #
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

# AO data
AO_data <- read_csv("06_processed_data/AO_data/AO_data_winter_spring_1990-2019.csv")

data_model <- CR_data %>%
  left_join(x = CR_data,
            y = AO_data,
            by = "year")


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
source("05_script/models/2_AO/2.1_JAGS_winter_AO.R")
source("05_script/models/functions_for_models.R")





# B. RUN THE MODELS ============================================================

# 1. winter AO: common slope (2.1.1.1_common) ----------------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "2.1.1.1"

# Predictor
var <- data_model$winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_s"
var_full_name <- "winter AO index"

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
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_", slope, mode),
       value = jags(data = dat, 
                    inits = inits, 
                    parameters.to.save = params, 
                    n.chains = 2, 
                    n.iter = 30000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))
end <- Sys.time()
end - start
# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




# 2. winter AO: distinct slopes (2.1.1.1_distinct) ----------------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "2.1.1.1"

# Predictor
var <- data_model$winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_s"
var_full_name <- "winter AO index"

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
                    n.iter = 50000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 3. Previous winter AO: common slope (2.1.2.1_common) ------------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "2.1.2.1"

# Predictor
var <- data_model$prior_winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_prior_s"
var_full_name <- "Previous winter AO index"

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
                    n.iter = 50000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




# 4. Previous winter AO: distinct slopes (2.1.2.1_distinct) -------------------

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "2.1.2.1"

# Predictor
var <- data_model$prior_winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_prior_s"
var_full_name <- "Previous winter AO index"

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
                    n.iter = 50000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 5. 2y prior winter AO: common slope (2.1.3.1_common) -------------------------

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <- "2.1.3.1"

# Predictor
var <- data_model$two_year_prior_winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_2y_prior_s"
var_full_name <- "2 year prior winter AO index"

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
                    n.iter = 50000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)



# 6. 2y prior winter AO: distinct slopes (2.1.3.1_distinct) ====================

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <-  "2.1.3.1"

# Predictor
var <- data_model$two_year_prior_winter_AO
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_2y_prior_s"
var_full_name <- "2 year prior winter AO index"

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
                    n.iter = 50000, n.burnin = 5000, 
                    n.thin = 10,
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

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
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)

# Plot
ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode,
   dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)





# 7. winter AO : distinct slopes, FACTOR (2.1.1.2_distinct) --------------------

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <-  "2.1.1.2"

# Predictor
var <- data_model$winter_AO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_pos"
var_full_name <- "winter AO phase"

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
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model ---------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

temp <- get_probabilities_factor(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

df.plot <- get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"))

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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

# Save
ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", 
                         model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode, dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




# 8. Previous winter AO : distinct slopes, FACTOR (2.1.2.2_distinct) ----------

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <-  "2.1.2.2"

# Predictor
var <- data_model$prior_winter_AO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_prior_pos"
var_full_name <- "Previous winter AO phase"

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
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model ---------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

temp <- get_probabilities_factor(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

df.plot <- get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"))

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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

# Save
ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", 
                         model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode, dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




# 9. Previous winter AO : distinct slopes, FACTOR (2.1.2.2_distinct) ----------

# ~~~ a. Run the model ---------------------------------------------------------
{model_code <-  "2.1.3.2"

# Predictor
var <- data_model$two_year_prior_winter_AO_pos
var_scaled <- (var - mean(var))/sd(var) 
var_short_name <- "w_AO_2y_prior_pos"
var_full_name <- "2 year prior winter AO phase"

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
                    model.file = get(paste0("model_", model_code, "_", slope, mode[1], "_slope"))))


# assign(x = paste0("DICw_", model_code, "_", slope, mode[1]),
#        value = get(paste0("fit_", model_code, "_", slope, mode[1]))$BUGSoutput$DIC)

save(list = paste0("fit_", model_code, "_", slope, mode), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))


# ~~~ b. Check convergence -----------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model ---------------------------------------------------------
load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_", slope, "_slope", toupper(mode), ".RData"))

temp <- get_probabilities_factor(model_code, mode, slope, var_scaled, var)
# Get df for ggplot
assign(x = paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
       value = temp[[1]])

df.plot <- get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"))

# Get the legend labels
color_labels <- temp[[2]]
rm(temp)


ggplot(data = get(paste0("fit_", model_code, "_", slope, mode[1], "_for_plot")), 
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

# Save
ggsave(filename = paste0("07_results/01_interim_results/model_outputs/graphs/model_", 
                         model_code, "_", slope, "_slope", toupper(mode), ".png"),
       width = 6, height = 3)


rm(list = c(paste0("fit_", model_code, "_", slope, mode),
            paste0("fit_", model_code, "_", slope, mode[1], "_for_plot"),
            paste0("DICw_", model_code, "_", slope, mode)))
rm(model_code, slope, mode, dat, params, nb.beta, inits, var, var_scaled, 
   var_short_name, var_full_name,  nb.beta)




