#==============================================================================#
#                                                                              #
#                  Models with ice-free days and spring AO                     #
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

setwd("D:/Litter_size_environment_project/")

# A. READY VARIABLES ===========================================================
# CR data
CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
sea_ice_data <- data.frame(sea_ice_data,
                           ice_free_days_previous = c(NA, sea_ice_data$ice_free_days[-nrow(sea_ice_data)]),
                           ice_free_days_2y_prior = c(NA, NA, sea_ice_data$ice_free_days[-c(nrow(sea_ice_data),
                                                                                            nrow(sea_ice_data) - 1)]))
AO_data <- read_csv("06_processed_data/AO_data/AO_data_winter_spring_1990-2019.csv")

CR_data %>%
  left_join(x = CR_data,
            y = sea_ice_data,
            by = "year") %>% 
  left_join(.,
            y = AO_data,
            by = "year") %>% 
  filter(year != 1992)  %>% # Remove the captures from 1992 since I can't calculate the ice free days in 1991
  filter(year != 1993)  -> data_model


# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number))
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


# Load the JAGS models + the ancillary functions
source("05_script/models/several_predictors/JAGS_multi_predictors.R")
source("05_script/models/functions_for_models.R")





# B. RUN THE MODELS ============================================================

# 1. Previous year ice-free days + prior spring AO: common slope ---------------

# ~~~ a. Run the model ---------------------------------------------------------

{pred1 <- "1.1.2_D"
pred2 <- "2.2.2.1"
model_code <- paste0(pred1, "__", pred2)

# Predictors
var1 <- data_model$ice_free_days_previous
var1_scaled <- (var1 - mean(var1))/sd(var1) 
var1_short_name <- "ice_free_days_previous_s"
var1_full_name <- "Ice-free days in previous year"

var2 <- data_model$prior_spring_AO
var2_scaled <- (var2 - mean(var2))/sd(var2) 
var2_short_name <- "s_AO_prior_s"
var2_full_name <- "Spring AO index in previous year"

# Number of slopes
slope1 <- "common"
slope2 <- "common"
slope <- paste0(slope1, "_", slope2)
# slope <- "distinct" 

# Are females without cubs taken into account ?
mode <- ""       # Yes
# mode <- "_bis"  # No
}

# Put all the data in a list
dat <- list(y, N, J, year, nbyear, var1_scaled, var2_scaled)
names(dat) <- c("y", "N", "J","year", "nbyear", var1_short_name, var2_short_name)

# Define the parameters to estimate
params <- c("a0", "a1", "a2", "b0", "c0", "sigma1", "eps1")
nb.beta <- length(params) - 2

# Generate starting values
temp.dat <- data.frame(y = y, var1_scaled = var1_scaled, var2_scaled = var2_scaled)
temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 

mlogit.mod <- mlogit(yfac ~ 1| var1_scaled + var2_scaled, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)
inits1 <- list(a0 = coefs[1], b0 = coefs[2], c0 = coefs[3],
               a1 = coefs[4], a2 = coefs[5],
               sigma1 = runif(1))
inits2 <- list(a0 = coefs[1] + 0.1, b0 = coefs[2] - 0.1, c0 = coefs[3] + 0.1,
               a1 = coefs[4] - 0.1, a2 = coefs[5] + 0.1,
               sigma1 = runif(1))
inits <- list(inits1, inits2)




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

res <- get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$sims.matrix

b1cub <- res[, c(1:3)]
b2cub <- res[, c(4, 2, 3)] 
b3cub <- res[, c(5, 2, 3)] 

range_var1 <- range(var1_scaled)

lengthgrid <- 100
grid_var1_scaled <- seq(from = range_var1[1] - 0.1*(range_var1[2] - range_var1[1]), 
                        to = range_var1[2] + 0.1*(range_var1[2] - range_var1[1]), 
                        length = lengthgrid)

grid_var1 <- grid_var1_scaled * sd(var1) + mean(var1)

q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(b1cub)[1], 
                                                               ncol = lengthgrid)
for (i in 1:lengthgrid) {
  for (j in 1:dim(b2cub)[1]) {
    q0cub_var1[j, i] <- exp(0)
    q1cub_var1[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_var1_scaled[i] + b1cub[j, 3] * mean(var2_scaled))	
    q2cub_var1[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_var1_scaled[i] + b2cub[j, 3] * mean(var2_scaled))	
    q3cub_var1[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_var1_scaled[i] + b3cub[j, 3] * mean(var2_scaled))		
  }
}

# Backtransform
p0cub_var1 <- p1cub_var1 <- p2cub_var1 <- p3cub_var1 <- matrix(NA, dim(b2cub)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(b2cub)[1]){
    norm <- (q0cub_var1[j, i] + q1cub_var1[j, i] + q2cub_var1[j, i] + q3cub_var1[j,i])
    p0cub_var1[j, i] <- q0cub_var1[j, i]/norm
    p1cub_var1[j, i] <- q1cub_var1[j, i]/norm
    p2cub_var1[j, i] <- q2cub_var1[j, i]/norm
    p3cub_var1[j, i] <- q3cub_var1[j, i]/norm
  }
}
df.for.plot.var1 <- data.frame(var = grid_var1,
                          mean_p_0_cub = apply(p0cub_var1, 2, mean),
                          mean_p_1_cub = apply(p1cub_var1, 2, mean),
                          mean_p_2_cub = apply(p2cub_var1, 2, mean),
                          mean_p_3_cub = apply(p3cub_var1, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub_var1, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub_var1, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub_var1, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub_var1, 2, quantile, probs = 0.975),
                          ci_p_2_cub_2.5 = apply(p2cub_var1, 2, quantile, probs = 0.025),
                          ci_p_2_cub_97.5 = apply(p2cub_var1, 2, quantile, probs = 0.975),
                          ci_p_3_cub_2.5 = apply(p3cub_var1, 2, quantile, probs = 0.025),
                          ci_p_3_cub_97.5 = apply(p3cub_var1, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_cub", "mean_p_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_cub_2.5", "ci_p_2_cub_97.5", 
                        "ci_p_3_cub_2.5", "ci_p_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_cub", "ci_p_2_cub_2.5", "ci_p_2_cub_97.5"), 2, 3))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_cub", "mean_p_3_cub"), "mean", "credible_interval"))

color_labels <- c("no cubs", "1 cub", "2 cubs", "3 cubs")



ggplot(data = df.for.plot.var1, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = color_labels) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = var1_full_name,
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






# 2. Previous year ice-free days + individual variables (Dorinda's model5A) ---------------

# ~~~ a. Run the model ---------------------------------------------------------

model <- "model_1.1.2_D_common_5A"
# Get variables
{var1 <- data_model$ice_free_days_previous
var1_scaled <- (var1 - mean(var1))/sd(var1) 
var1_short_name <- "ice_free_days_previous_s"
var1_full_name <- "Ice-free days in previous year"

var2 <- data_model$age_for_analyses
var2_scaled <- (var2 - mean(var2))/sd(var2) 
var2_short_name <- "age_s"
var2_full_name <- "Age"

var3 <- data_model$age_for_analyses^2
var3_scaled <- (var3 - mean(var3))/sd(var3) 
var3_short_name <- "age2_s"
var3_full_name <- "Age squared"

var4 <- data_model$s_length
var4_scaled <- (var4 - mean(var4))/sd(var4) 
var4_short_name <- "size_s"
var4_full_name <- "Size"

var5 <- data_model$s_length^2
var5_scaled <- (var5 - mean(var5))/sd(var5) 
var5_short_name <- "size2_s"
var5_full_name <- "Size squared"

var6 <- data_model$day_number
var6_scaled <- (var6 - mean(var6))/sd(var6) 
var6_short_name <- "day_s"
var6_full_name <- "Day of capture"
}
dat <- list(y, N, J, year, nbyear, var1_scaled, var2_scaled, var3_scaled, var4_scaled, 
            var5_scaled, var6_scaled)
names(dat) <- c("y", "N", "J","year","nbyear", var1_short_name, var2_short_name, 
                var3_short_name, var4_short_name, var5_short_name, var6_short_name)
params <- c("b","sigma1","eps1") 
nb.beta = 9


# Generate starting values
temp.dat <- data.frame(y = y, var1_scaled = var1_scaled, var2_scaled = var2_scaled,
                       var3_scaled = var3_scaled, var4_scaled = var4_scaled,
                       var5_scaled = var5_scaled, var6_scaled = var6_scaled)
temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 

mlogit.mod <- mlogit(yfac ~ 1| var1_scaled + var2_scaled + var3_scaled + var4_scaled + var5_scaled + var6_scaled, 
                     data = mnl.dat, 
                     reflevel = "0")
coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits1 <- list(b = coefs[c(1, 4, 7, 10, 13, 16, 19, 2, 3)], sigma1 = runif(1)) 
inits2 <- list(b = (coefs[c(1, 4, 7, 10, 13, 16, 19, 2, 3)] + 0.1), sigma1 = runif(1)) 
inits <- list(inits1, inits2)

start <- Sys.time()
fit_1.1.2_D_common_5A <- jags(data = dat, 
                              inits = inits, 
                              params, 
                              n.chains = 2, 
                              n.iter = 200000, 
                              n.burnin = 40000, 
                              n.thin = 10,
                              model.file = model_1.1.2_D_common_5A)
end <- Sys.time()
end - start

save(fit_1.1.2_D_common_5A, file = "07_results/01_interim_results/model_outputs/fit_1.1.2_D_common_5A.RData")

temp.fit_1.1.2_D_common_5A <- fit_1.1.2_D_common_5A



# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/fit_1.1.2_D_common_5A.RData")

processed_output <- ggs(as.mcmc(fit_1.1.2_D_common_5A)) %>%
  filter(Parameter %in% c("b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]",
                          "b[7]", "b[8]", "b[9]", "b[10]", "deviance"))

f1 <- ggs_traceplot(processed_output) + 
  theme_bw() +
  theme(legend.position = "none")
f2 <- ggs_density(processed_output) + 
  theme_bw()
f3 <- ggs_running(processed_output) + 
  theme_bw() +
  theme(legend.position = "none")
x <- grid.arrange(f1, f2, f3, ncol = 3, nrow = 1)

nbr_rows <- length(unique(processed_output$Parameter))

ggsave(x,
       filename = paste0("07_results/01_interim_results/model_outputs/graphs/diagnostic_plots/model_1.1.2_D_common_5A.png"),
       width = 17.5, height = 2.5*nbr_rows)

rm(f1, f2, f3, x, processed_output)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/fit_1.1.2_D_common_5A.RData")

res <- fit_1.1.2_D_common_5A$BUGSoutput$sims.matrix


b1cub <- res[, c(1:7)]
b2cub <- res[, c(8, 2:7)] 
b3cub <- res[, c(9, 2, 5)] 

range_var1 <- range(var1_scaled)

lengthgrid <- 100
grid_var1_scaled <- seq(from = range_var1[1] - 0.1*(range_var1[2] - range_var1[1]), 
                        to = range_var1[2] + 0.1*(range_var1[2] - range_var1[1]), 
                        length = lengthgrid)

grid_var1 <- grid_var1_scaled * sd(var1) + mean(var1)

q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(b1cub)[1], 
                                                               ncol = lengthgrid)

# Mean for all variables but sea ice
{specificity <- ""
for (i in 1:lengthgrid) {
  for (j in 1:dim(b2cub)[1]) {
    q0cub_var1[j, i] <- exp(0)
    
    q1cub_var1[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_var1_scaled[i] + b1cub[j, 3] * mean(var2_scaled) +
                              b1cub[j, 4] * mean(var3_scaled) + b1cub[j, 5] * mean(var4_scaled) +
                              b1cub[j, 6] * mean(var5_scaled) + b1cub[j, 7] * mean(var6_scaled))
    
    q2cub_var1[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_var1_scaled[i] + b2cub[j, 3] * mean(var2_scaled) +
                              b2cub[j, 4] * mean(var3_scaled) + b2cub[j, 5] * mean(var4_scaled) +
                              b2cub[j, 6] * mean(var5_scaled) + b2cub[j, 7] * mean(var6_scaled))	
    q3cub_var1[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_var1_scaled[i] + b3cub[j, 3] * mean(var4_scaled))		
  }
}
}


# Mean for all variables but sea ice and but AGE (and age2): first quartile
specificity <- "age_first_quartile"
fst_quartile_age <- third.quart.adequacy <- quantile(x = var2_scaled, 
                                                      probs = 0.25, 
                                                      na.rm = TRUE)
fst_quartile_age2 <- third.quart.adequacy <- quantile(x = var3_scaled, 
                                                     probs = 0.25, 
                                                     na.rm = TRUE)
(6 - mean(var2))/sd(var2)
for (i in 1:lengthgrid) {
  for (j in 1:dim(b2cub)[1]) {
    q0cub_var1[j, i] <- exp(0)
    
    q1cub_var1[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_var1_scaled[i] + b1cub[j, 3] * fst_quartile_age +
                              b1cub[j, 4] * fst_quartile_age2 + b1cub[j, 5] * mean(var4_scaled) +
                              b1cub[j, 6] * mean(var5_scaled) + b1cub[j, 7] * mean(var6_scaled))
    
    q2cub_var1[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_var1_scaled[i] + b2cub[j, 3] * fst_quartile_age +
                              b2cub[j, 4] * fst_quartile_age2 + b2cub[j, 5] * mean(var4_scaled) +
                              b2cub[j, 6] * mean(var5_scaled) + b2cub[j, 7] * mean(var6_scaled))	
    q3cub_var1[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_var1_scaled[i] + b3cub[j, 3] * mean(var4_scaled))		
  }
}

# Mean for all variables but sea ice and but AGE (and age2): third quartile
{specificity <- "age_third_quartile"
third_quartile_age <- third.quart.adequacy <- quantile(x = var2_scaled, 
                                                     probs = 0.75, 
                                                     na.rm = TRUE)
third_quartile_age2 <- third.quart.adequacy <- quantile(x = var3_scaled, 
                                                       probs = 0.75, 
                                                       na.rm = TRUE)
start <- Sys.time()
for (i in 1:lengthgrid) {
  for (j in 1:dim(b2cub)[1]) {
    q0cub_var1[j, i] <- exp(0)
    
    q1cub_var1[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_var1_scaled[i] + b1cub[j, 3] * third_quartile_age +
                              b1cub[j, 4] * third_quartile_age2 + b1cub[j, 5] * mean(var4_scaled) +
                              b1cub[j, 6] * mean(var5_scaled) + b1cub[j, 7] * mean(var6_scaled))
    
    q2cub_var1[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_var1_scaled[i] + b2cub[j, 3] * third_quartile_age +
                              b2cub[j, 4] * third_quartile_age2 + b2cub[j, 5] * mean(var4_scaled) +
                              b2cub[j, 6] * mean(var5_scaled) + b2cub[j, 7] * mean(var6_scaled))	
    q3cub_var1[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_var1_scaled[i] + b3cub[j, 3] * mean(var4_scaled))		
  }
}
end <- Sys.time()
end - start 
}


# Mean for all variables but sea ice and but SIZE (and size2): third quartile
{specificity <- "size_third_quartile"
third_quartile_size <- third.quart.adequacy <- quantile(x = var4_scaled, 
                                                       probs = 0.75, 
                                                       na.rm = TRUE)
third_quartile_size2 <- third.quart.adequacy <- quantile(x = var5_scaled, 
                                                       probs = 0.75, 
                                                       na.rm = TRUE)
  start <- Sys.time()
  for (i in 1:lengthgrid) {
    for (j in 1:dim(b2cub)[1]) {
      q0cub_var1[j, i] <- exp(0)
      
      q1cub_var1[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_var1_scaled[i] + b1cub[j, 3] * mean(var2_scaled) +
                                b1cub[j, 4] * mean(var3_scaled) + b1cub[j, 5] * third_quartile_size +
                                b1cub[j, 6] * third_quartile_size2 + b1cub[j, 7] * mean(var6_scaled))
      
      q2cub_var1[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_var1_scaled[i] + b2cub[j, 3] * mean(var2_scaled) +
                                b2cub[j, 4] * mean(var3_scaled) + b2cub[j, 5] * third_quartile_size +
                                b2cub[j, 6] * third_quartile_size2 + b2cub[j, 7] * mean(var6_scaled))	
      q3cub_var1[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_var1_scaled[i] + b3cub[j, 3] * third_quartile_size)		
    }
  }
  end <- Sys.time()
  end - start 
}






# Back-transform
p0cub_var1 <- p1cub_var1 <- p2cub_var1 <- p3cub_var1 <- matrix(NA, dim(b2cub)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(b2cub)[1]){
    norm <- (q0cub_var1[j, i] + q1cub_var1[j, i] + q2cub_var1[j, i] + q3cub_var1[j,i])
    p0cub_var1[j, i] <- q0cub_var1[j, i]/norm
    p1cub_var1[j, i] <- q1cub_var1[j, i]/norm
    p2cub_var1[j, i] <- q2cub_var1[j, i]/norm
    p3cub_var1[j, i] <- q3cub_var1[j, i]/norm
  }
}
df.for.plot.var1 <- data.frame(var = grid_var1,
                               mean_p_0_cub = apply(p0cub_var1, 2, mean),
                               mean_p_1_cub = apply(p1cub_var1, 2, mean),
                               mean_p_2_cub = apply(p2cub_var1, 2, mean),
                               mean_p_3_cub = apply(p3cub_var1, 2, mean),
                               ci_p_0_cub_2.5 = apply(p0cub_var1, 2, quantile, probs = 0.025),
                               ci_p_0_cub_97.5 = apply(p0cub_var1, 2, quantile, probs = 0.975),
                               ci_p_1_cub_2.5 = apply(p1cub_var1, 2, quantile, probs = 0.025),
                               ci_p_1_cub_97.5 = apply(p1cub_var1, 2, quantile, probs = 0.975),
                               ci_p_2_cub_2.5 = apply(p2cub_var1, 2, quantile, probs = 0.025),
                               ci_p_2_cub_97.5 = apply(p2cub_var1, 2, quantile, probs = 0.975),
                               ci_p_3_cub_2.5 = apply(p3cub_var1, 2, quantile, probs = 0.025),
                               ci_p_3_cub_97.5 = apply(p3cub_var1, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_cub", "mean_p_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_cub_2.5", "ci_p_2_cub_97.5", 
                        "ci_p_3_cub_2.5", "ci_p_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_cub", "ci_p_2_cub_2.5", "ci_p_2_cub_97.5"), 2, 3))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_cub", "mean_p_3_cub"), "mean", "credible_interval"))

color_labels <- c("no cubs", "1 cub", "2 cubs", "3 cubs")



ggplot(data = df.for.plot.var1, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = color_labels) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = var1_full_name,
       y = "Probability", 
       color = "",
       linetype = "")

ggsave(paste0("07_results/01_interim_results/model_outputs/graphs/", model, "_", specificity, ".png"),
       width = 6, height = 3)
