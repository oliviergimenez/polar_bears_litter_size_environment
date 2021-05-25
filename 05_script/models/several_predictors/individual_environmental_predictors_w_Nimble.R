#==============================================================================#
#                                                                              #
#              Models with individual + envrionmental predictors               #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(nimble)
library(viridis)
library(ggmcmc)
library(gridExtra)
library(MCMCvis)
Sys.setenv(LANG = "en")


# A. READY VARIABLES ===========================================================

# ~~~ a. Load the data ---------------------------------------------------------
CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Sea ice data
sea_ice_data <- read_csv("06_processed_data/sea_ice_data/retreat_advance_ice_free_days_D.csv")
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
  filter(year != 1992)  %>% # Remove the captures from 1992 since I can't calculate the ice free days in 1991
  filter(year != 1993)  -> data_model

rm(CR_data, sea_ice_data, AO_data, NAO_data)

# ~~~ b. Format data for Nimble ------------------------------------------------

y <- factor(as.numeric(data_model$cub_number_2))
summary(y)

# Individual traits +++++++++++++++++++++++++++++++++

# Size
size <- data_model$s_length
size_s <- (size - mean(size))/sd(size) 

# Size squarred
size_s2 <- size_s^2


# Age (in 3 categories)
age <- as.numeric(data_model$age_for_analyses)
age_factor <- as.factor(ifelse (age < 9, "b. young", 
                                ifelse(age > 15, "c. old", "a. prime_age")))
relevel(age_factor, ref = "a. prime_age")

age_old <- ifelse (age > 15, 1, 0)
age_young <- ifelse (age < 9, 1, 0)


# Other variables +++++++++++++++++++++++++++++++

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- (day - mean(day))/sd(day)

# Number of ice-free days t-1

ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- (ice_free_days_previous - mean(ice_free_days_previous))/sd(ice_free_days_previous)


# Random effects ++++++++++++++++++++++++++++++++++++

# years
{year <- data_model$year
  year <- factor(year) # laisse tomber les modalites inutiles 
  year2 <- NULL
  for (i in 1:length(year)) {
    year2 <- c(year2, which(year[i] == levels(year)))
  }
  year <- factor(year2)
  rm(year2)
}

# females
{id_fem <- data_model$ID_NR
  id_fem <- factor(id_fem) 
  id_fem2 <- NULL
  for (i in 1:length(id_fem)){
    id_fem2 <- c(id_fem2,which(id_fem[i]==levels(id_fem)))
  }
  id_fem <- factor(id_fem2)
  # nbind <- length(levels(id_fem))
  rm(id_fem2)
}

dat <- list(y = as.numeric(y))

# Load the JAGS models + the ancillary functions
source("05_script/models/several_predictors/Nimble_individual_envrionmental_predictors.R")
source("05_script/models/functions_for_models.R")


# B. RUN THE MODELS ============================================================

# 1. Previous year ice-free days + individual variables: common slopes ---------

my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = length(levels(year)),
                     ice_free_days_previous_s = as.numeric(ice_free_days_previous_s), 
                     size_s = as.numeric(size_s), size_s2 = as.numeric(size_s2), 
                     age_young = as.factor(age_young), age_old = as.factor(age_old), 
                     day_s = as.numeric(day_s))

params <- c("b", "sigma1" ,"eps1")

# Generate initial values
temp.dat <- data.frame(y = y, 
                       age_factor = age_factor,
                       size_s = size_s, 
                       size_s2 = size_s2,
                       day_s = day_s,
                       ice_free_days_previous_s = ice_free_days_previous_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| ice_free_days_previous_s + size_s + size_s2 + age_factor * day_s, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)


inits <- function() list(b = coefs[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 2)] + round(runif(n = 1, -1, 1))/10, 
                         sigma1 = runif(1))


# ~~~ a. Run the model ---------------------------------------------------------

start <- Sys.time()
fit_model_1 <- nimbleMCMC(code = model_1,     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits,          
                          monitors = params,   # parameters to monitor
                          thin = 10,
                          niter = 220000,                  # nb iterations
                          nburnin = 44000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE)
end <- Sys.time()
end - start # 24min

fit_model_1$WAIC
# 1074.887
save(fit_model_1, 
     file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.RData")

# ~~~ b. Check convergence -----------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.RData")

MCMCsummary(object = fit_model_1, round = 2)
MCMCtrace(object = fit_model_1, params = b, pdf = FALSE)
?MCMCtrace
check_convergence(jags_output = get(paste0("fit_", model_code, "_", slope, mode)),
                  model_code = model_code, 
                  slope = slope, mode = mode)



# ~~~ c. Plot the model --------------------------------------------------------
load(file = "07_results/01_interim_results/model_outputs/several_predictors/model_1.RData")

res <- rbind(fit_model_1$samples$chain1,
             fit_model_1$samples$chain2)
b1cub <- res[, c(1:9)]
b2_3cub <- res[, c(10, 2:9)] 

# ~~~~~~ Sea ice ---------------------------------------------------
range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)

q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                    nrow = dim(res)[1], 
                                    ncol = lengthgrid)
q0cub_var1 <- q1cub_var1 <- q2cub_var1 <- q3cub_var1 <- matrix(data = NA, 
                                                               nrow = dim(res)[1], 
                                                               ncol = lengthgrid)

# ~~~~~~~~~ Ages -----------------------------------
# Prime-aged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 0 * mean(day_s) + b1cub[j, 9] * 0 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 0 * mean(day_s) + b2_3cub[j, 9] * 0 * mean(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.prime.aged <- data.frame(var = grid,
                          mean_p_0_cub = apply(p0cub, 2, mean),
                          mean_p_1_cub = apply(p1cub, 2, mean),
                          mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                          ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                          ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         age = "prime-aged")


# Young ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

start <- Sys.time() 
for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 1 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 1 * mean(day_s) + b1cub[j, 9] * 0 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 1 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 1 * mean(day_s) + b2_3cub[j, 9] * 0 * mean(day_s))	
  }
}
end <- Sys.time() 
end - start

start <- Sys.time() 
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}
end <- Sys.time() 
end - start

df.for.plot.young <- data.frame(var = grid,
                          mean_p_0_cub = apply(p0cub, 2, mean),
                          mean_p_1_cub = apply(p1cub, 2, mean),
                          mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                          ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                          ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         age = "young")


# Old ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 1 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 0 * mean(day_s) + b1cub[j, 9] * 1 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 1 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 0 * mean(day_s) + b2_3cub[j, 9] * 1 * mean(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.old <- data.frame(var = grid,
                                mean_p_0_cub = apply(p0cub, 2, mean),
                                mean_p_1_cub = apply(p1cub, 2, mean),
                                mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", "credible_interval"),
         age = "old")



df.for.plot.age <- rbind(df.for.plot.young, df.for.plot.prime.aged, df.for.plot.old)
write_csv(df.for.plot.age, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_age.csv")
df.for.plot.age <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_age.csv")

ggplot(data = df.for.plot.age, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap(~ age)

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet age.png",
       width = 14, height = 6)


# ~~~~~~~~~ Sizes -----------------------------------

# Medium sized ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 0 * mean(day_s) + b1cub[j, 9] * 0 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 0 * mean(day_s) + b2_3cub[j, 9] * 0 * mean(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.medium.sized <- data.frame(var = grid,
                                     mean_p_0_cub = apply(p0cub, 2, mean),
                                     mean_p_1_cub = apply(p1cub, 2, mean),
                                     mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                     ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                     ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                     ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                     ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                     ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                     ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         size = "medium")

# df.for.plot.medium.sized <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_age.csv") %>%
#   filter(age == "prime-aged") %>%
#   select(-age) %>%
#   mutate(size = "medium")


# Small ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
fst_decile <- quantile(x = size_s, # = 187 kg
                       probs = 0.10, 
                       na.rm = TRUE)

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * fst_decile + b1cub[j, 4] * fst_decile^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 0 * mean(day_s) + b1cub[j, 9] * 0 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * fst_decile + b2_3cub[j, 4] * fst_decile^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 0 * mean(day_s) + b2_3cub[j, 9] * 0 * mean(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}


df.for.plot.small <- data.frame(var = grid,
                                     mean_p_0_cub = apply(p0cub, 2, mean),
                                     mean_p_1_cub = apply(p1cub, 2, mean),
                                     mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                     ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                     ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                     ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                     ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                     ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                     ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         size = "small")


# Large ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
lst_decile <- quantile(x = size_s, # = 187 kg
                       probs = 0.90, 
                       na.rm = TRUE)


for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * lst_decile + b1cub[j, 4] * lst_decile^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * mean(day_s) +
                         b1cub[j, 8] * 0 * mean(day_s) + b1cub[j, 9] * 0 * mean(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * lst_decile + b2_3cub[j, 4] * lst_decile^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * mean(day_s) +
                           b2_3cub[j, 8] * 0 * mean(day_s) + b2_3cub[j, 9] * 0 * mean(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.large <- data.frame(var = grid,
                                mean_p_0_cub = apply(p0cub, 2, mean),
                                mean_p_1_cub = apply(p1cub, 2, mean),
                                mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         size = "large")



df.for.plot.size <- rbind(df.for.plot.small, df.for.plot.medium.sized, df.for.plot.large) %>%
  mutate(size = factor(size, levels = c("small", "medium", "large")))
         
write_csv(df.for.plot.size, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_size.csv")
df.for.plot.size <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_size.csv")

ggplot(data = df.for.plot.size, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap(~ size)

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet size.png",
       width = 14, height = 6)



# ~~~~~~~~~ Day of capture -----------------------------------

# Median day of capture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * median(day_s) +
                         b1cub[j, 8] * 0 * median(day_s) + b1cub[j, 9] * 0 * median(day_s))	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * median(day_s) +
                           b2_3cub[j, 8] * 0 * median(day_s) + b2_3cub[j, 9] * 0 * median(day_s))	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.mid.season <- data.frame(var = grid,
                                     mean_p_0_cub = apply(p0cub, 2, mean),
                                     mean_p_1_cub = apply(p1cub, 2, mean),
                                     mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                     ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                     ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                     ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                     ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                     ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                     ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         day = "mid season")


# Early season ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
fst_decile <- quantile(x = day_s, # = day 90 of the year
                       probs = 0.10, 
                       na.rm = TRUE)

for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * fst_decile +
                         b1cub[j, 8] * 0 * fst_decile + b1cub[j, 9] * 0 * fst_decile)	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * fst_decile +
                           b2_3cub[j, 8] * 0 * fst_decile + b2_3cub[j, 9] * 0 * fst_decile)	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}


df.for.plot.early <- data.frame(var = grid,
                                mean_p_0_cub = apply(p0cub, 2, mean),
                                mean_p_1_cub = apply(p1cub, 2, mean),
                                mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         day = "early")


# Large ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
lst_decile <- quantile(x = day_s, # = day 116 of the year
                       probs = 0.90, 
                       na.rm = TRUE)


for (i in 1:lengthgrid) {
  for (j in 1:dim(res)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + 
                         b1cub[j, 2] * grid_scaled[i] +
                         b1cub[j, 3] * mean(size_s) + b1cub[j, 4] * mean(size_s)^2 +
                         b1cub[j, 5] * 0 + b1cub[j, 5] * 0 +
                         b1cub[j, 7] * lst_decile +
                         b1cub[j, 8] * 0 * lst_decile + b1cub[j, 9] * 0 * lst_decile)	
    q2_3cub[j, i] <- exp(b2_3cub[j, 1] + 
                           b2_3cub[j, 2] * grid_scaled[i] +
                           b2_3cub[j, 3] * mean(size_s) + b2_3cub[j, 4] * mean(size_s)^2 +
                           b2_3cub[j, 5] * 0 + b2_3cub[j, 5] * 0 +
                           b2_3cub[j, 7] * lst_decile +
                           b2_3cub[j, 8] * 0 * lst_decile + b2_3cub[j, 9] * 0 * lst_decile)	
  }
}
p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(res)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2_3cub[j, i] <- q2_3cub[j, i]/norm
  }
}

df.for.plot.late <- data.frame(var = grid,
                                mean_p_0_cub = apply(p0cub, 2, mean),
                                mean_p_1_cub = apply(p1cub, 2, mean),
                                mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                                ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                                ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                                ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                                ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                                ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                                ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
  pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                        "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                        "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                        "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
  mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                             ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                    ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
         type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", 
                       "credible_interval"),
         day = "late")



df.for.plot.size <- rbind(df.for.plot.small, df.for.plot.medium.sized, df.for.plot.large) %>%
  mutate(size = factor(size, levels = c("small", "medium", "large")))

write_csv(df.for.plot.size, 
          "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_size.csv")
df.for.plot.size <- read_csv("07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet_size.csv")

ggplot(data = df.for.plot.size, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = c("no cubs", "1 cub", "2-3 cubs")) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Number of ice-free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") +
  facet_wrap(~ size)

ggsave(filename = "07_results/01_interim_results/model_outputs/several_predictors/graphs/model_1_ice_free_days_facet size.png",
       width = 14, height = 6)











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
fst_quartile_age <- quantile(x = var2_scaled, 
                                                     probs = 0.25, 
                                                     na.rm = TRUE)
fst_quartile_age2 <- quantile(x = var3_scaled, 
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
