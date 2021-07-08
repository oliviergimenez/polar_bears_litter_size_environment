#==============================================================================#
#                                                                              #
#            Variable selection using the reversible jump MCMC method          #
#                                 (multinomial)                                #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")

# Load the data ----------------------------------------------------------------

CR_data <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

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


# Format data for Nimble ------------------------------------------------

y <- data_model %>%
  pull(cub_number_2)

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- as.vector(scale(day))

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

# ~ 1. With 1 predictors -------------------------------------------------------

# ~~~ a. Write the model -------------------------------------------------------

glmIndicatorCode <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- a0 +                             
      zbeta[1]*day_s[i]                                           
    
    log(q[i, 3]) <- b0 + 
      zbeta[2]*day_s[i]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  a0 ~ dnorm(0, sd = 1.5)
  b0 ~ dnorm(0, sd = 1.5)
  
  psi ~ dunif(0, 1) 
  for(i in 1:2){
    z[i] ~ dbern(psi)               
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]      
    
  } 
})


# ~~~ b. Constants and data ----------------------------------------------------

## Set up the model.
constants <- list(N = length(y),
                  J = length(unique(y)))

dat  <- list(y = y,
             day_s = day_s)

# Inits
temp.dat <- data.frame(y = y, 
                       day_s = day_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur

mnl.dat <- mlogit.data(temp.dat, varying = NULL, 
                       choice = "yfac", shape = "wide")

mlogit.mod <- mlogit(yfac ~ 1| day_s, 
                     data = mnl.dat, 
                     reflevel = "0")
summary(mlogit.mod)
coefs <- as.vector(summary(mlogit.mod)$coefficients)

inits <- list(psi = 0.5,
              a0 = coefs[1],
              b0 = coefs[2],
              beta = coefs[c(3, 4)],
              z = sample(x = 0:1, size = 2, replace = TRUE))

glmIndicatorModel <- nimbleModel(code = glmIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)


# ~~~ c. Configuring RJMCMC  ---------------------------------------------------
glmIndicatorConf <- configureMCMC(glmIndicatorModel)
glmIndicatorConf$addMonitors('z')

configureRJ(conf = glmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))


# ~~~ d. Build and run the RJMCMC ----------------------------------------------

mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(reversible_jump <- runMCMC(CMCMCIndicatorRJ, 
                                       niter = 200000, 
                                       nburnin = 10000,
                                       nchains = 2))



# ~~~ e. Individual inclusion probabilities ------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("intercept", "day number", "age young", "age old",
                                                     "size", "size2", "day*age young",
                                                     "day*age old"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("intercept", "day number", "age young", "age old",
                              "size", "size2", "day*age young",
                              "day*age old"))


# ~~~ f. Parameter distributions -----------------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])
betaCols <- grep("beta\\[", colnames(allchains))

chains_l <- as.data.frame(allchains[, c(betaCols, max(betaCols) + 1)]) %>%
  mutate(iteration = rep(seq(1:dim(reversible_jump[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(reversible_jump[["chain1"]])[1]),
                   rep(2, times = dim(reversible_jump[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
                        "beta[6]", "beta[7]", "beta[8]", "psi"),
               names_to = "parameter")

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(reversible_jump[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(reversible_jump[["chain1"]])[1]),
                   rep(2, times = dim(reversible_jump[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = colnames(allchains),
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

nrows = length(colnames(allchains))
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)




#______________________----
# ~ 1. With 10 predictors -------------------------------------------------------

# ~~~ a. Write the model -------------------------------------------------------

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0,1)                  # prior on inclusion probability
  
  for(i in 1:(numVars)) {
    z[i] ~ dbern(psi)               # indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]      # indicator * beta
  }
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- zbeta[1] +
      zbeta[2] * day_s[i] +
      zbeta[3] * age_young[i] +
      zbeta[4] * age_old[i] +
      zbeta[5] * size_s[i] +
      zbeta[6] * size_s2[i] +
      zbeta[7] * day_s[i] * age_young[i] +
      zbeta[8] * day_s[i] * age_old[i] +
      zbeta[9] * ice_free_days_previous_s[i] +
      zbeta[10] * winter_AO_s[i] + 
      zbeta[11] * prior_spring_AO_s[i]
    
    log(q[i, 3]) <- zbeta[12] +
      zbeta[13] * day_s[i] +
      zbeta[14] * age_young[i] +
      zbeta[15] * age_old[i] +
      zbeta[16] * size_s[i] +
      zbeta[17] * size_s2[i] +
      zbeta[18] * day_s[i] * age_young[i] +
      zbeta[19] * day_s[i] * age_old[i] +
      zbeta[20] * ice_free_days_previous_s[i] +
      zbeta[21] * winter_AO_s[i] + 
      zbeta[22] * prior_spring_AO_s[i]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
})


# ~~~ b. Constants and data ----------------------------------------------------

## Set up the model.
constants <- list(N = length(y),
                  J = length(unique(y)),
                  numVars = 22)

dat  <- list(y = y, 
             day_s = day_s,
             age_young = age_young,
             age_old = age_old,
             size_s = size_s,
             size_s2 = size_s2,
             ice_free_days_previous_s = ice_free_days_previous_s,
             winter_AO_s = winter_AO_s,
             prior_spring_AO_s = prior_spring_AO_s)

# Inits
temp.dat <- data.frame(y = y, 
                       day_s = day_s,
                       age_young = age_young,
                       age_old = age_old,
                       size_s = size_s,
                       size_s2 = size_s2,
                       ice_free_days_previous_s = ice_free_days_previous_s,
                       winter_AO_s = winter_AO_s,
                       prior_spring_AO_s = prior_spring_AO_s)

temp.dat$yfac <- as.factor(temp.dat$y) # Ajout de Y en facteur

mnl.dat <- mlogit.data(temp.dat, varying = NULL, 
                       choice = "yfac", shape = "wide")

mlogit.mod <- mlogit(yfac ~ 1| day_s * (age_young + age_old) + size_s + size_s2 +
                       ice_free_days_previous_s  + winter_AO_s + prior_spring_AO_s, 
                     data = mnl.dat, 
                     reflevel = "0")
summary(mlogit.mod)
coefs <- as.vector(summary(mlogit.mod)$coefficients)

inits <- list(psi = 0.5,
              beta = coefs[c(seq(1, 21, 2), seq(2, 22, 2))],
              z = sample(0:1, constants$numVars, 0.5))

glmIndicatorModel <- nimbleModel(code = glmIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)


# ~~~ c. Configuring RJMCMC  ---------------------------------------------------
glmIndicatorConf <- configureMCMC(glmIndicatorModel)
glmIndicatorConf$addMonitors('z')

configureRJ(conf = glmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

# ~~~ d. Build and run the RJMCMC ----------------------------------------------

mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(reversible_jump <- runMCMC(CMCMCIndicatorRJ, 
                                       niter = 200000, 
                                       nburnin = 10000,
                                       nchains = 2))



# ~~~ e. Individual inclusion probabilities ------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("intercept", "day number", "age young", "age old",
                                                     "size", "size2", "day*age young",
                                                     "day*age old"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("intercept", "day number", "age young", "age old",
                              "size", "size2", "day*age young",
                              "day*age old"))


# ~~~ f. Parameter distributions -----------------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])
betaCols <- grep("beta\\[", colnames(allchains))

chains_l <- as.data.frame(allchains[, c(betaCols, max(betaCols) + 1)]) %>%
  mutate(iteration = rep(seq(1:dim(reversible_jump[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(reversible_jump[["chain1"]])[1]),
                   rep(2, times = dim(reversible_jump[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
                        "beta[6]", "beta[7]", "beta[8]", "psi"),
               names_to = "parameter")

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(reversible_jump[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(reversible_jump[["chain1"]])[1]),
                   rep(2, times = dim(reversible_jump[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = colnames(allchains),
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

nrows = length(colnames(allchains))
save_plot(filename = 
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)
#______________________----
# ~ 2. With 6 predictors -------------------------------------------------------






#______________________----
# ~ 2. With 6 predictors -------------------------------------------------------





