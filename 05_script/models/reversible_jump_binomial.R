#==============================================================================#
#                                                                              #
#            Variable selection using the reversible jump MCMC method          #
#                                   (binomial)                                 #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")

# Load the data ----------------------------------------------------------------

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


# Format data for Nimble ------------------------------------------------

y <- data_model %>%
  pull(success)

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


# ~ 1. With 8 predictors -------------------------------------------------------

# ~~~ a. Write the model -------------------------------------------------------

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- zbeta[1] +
      zbeta[2] * day_s[i] +
      zbeta[3] * age_young[i] +
      zbeta[4] * age_old[i] +
      zbeta[5] * size_s[i] +
      zbeta[6] * size_s2[i] +
      zbeta[7] * day_s[i] * age_young[i] +
      zbeta[8] * day_s[i] * age_old[i]
    
  }
})



# ~~~ b. Constants and data ----------------------------------------------------

## Set up the model.
constants <- list(N = length(y), 
                  numVars = 8)
dat  <- list(y = y, 
             day_s = day_s,
             age_young = age_young,
             age_old = age_old,
             size_s = size_s,
             size_s2 = size_s2)

# Inits
logit_model <- glm(formula = y ~ day_s * (age_young + age_old) 
                   + size_s + size_s2,
                   family = binomial(link = "logit"))

logit_model$coefficients
summary(logit_model)

inits <- list(psi = 0.5,
              beta = logit_model$coefficients,
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

## Check the assigned samplers
glmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))
glmIndicatorConf$printSamplers(c("z[5]", "beta[5]"))


# ~~~ d. Build and run the RJMCMC ----------------------------------------------

mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(reversible_jump <- runMCMC(CMCMCIndicatorRJ, 
                                        niter = 150000, 
                                        nburnin = 10000,
                                        nchains = 2))

save(reversible_jump,
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/reversible_jump.RData")


# ~~~ e. Individual inclusion probabilities ------------------------------------

load("07_results/01_interim_results/model_outputs/covariate_model_selection/reversible_jump.RData")

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
                              "day*age old")) +
  scale_y_continuous(limits = c(0, 1))

ggsave("07_results/01_interim_results/model_outputs/covariate_model_selection/graphs/incl_prob_reversible_jump.png",
       height = 4, width = 5)

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

# ~ 2. With 11 predictors -------------------------------------------------------

# ~~~ a. Write the model -------------------------------------------------------

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0, 1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- zbeta[1] +
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
  }
})


# ~~~ b. Constants and data ----------------------------------------------------

## Set up the model.
constants <- list(N = length(y), 
                  numVars = 11)
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
logit_model <- glm(formula = y ~ day_s * (age_young + age_old) + 
                   size_s + size_s2 + 
                   ice_free_days_previous_s + winter_AO_s + prior_spring_AO_s,
                   family = binomial(link = "logit"))

logit_model$coefficients
summary(logit_model)

inits <- list(psi = 0.5,
              beta = logit_model$coefficients,
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

## Check the assigned samplers
glmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))
glmIndicatorConf$printSamplers(c("z[5]", "beta[5]"))


# ~~~ d. Build and run the RJMCMC ----------------------------------------------

mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(reversible_jump <- runMCMC(CMCMCIndicatorRJ, 
                                       niter = 150000, 
                                       nburnin = 5000,
                                       nchains = 2))

save(reversible_jump,
     file = "07_results/01_interim_results/model_outputs/covariate_model_selection/reversible_jump_2.RData")


# ~~~ e. Individual inclusion probabilities ------------------------------------

load("07_results/01_interim_results/model_outputs/covariate_model_selection/reversible_jump_2.RData")

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("intercept", "day number", "age young", "age old",
                                                     "size", "size2", "day*age young",
                                                     "day*age old", "ice free days t-1",
                                                     "winter AO", "spring AO t-1"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = c("intercept", "day number", "age young", "age old",
                              "size", "size2", "day*age young",
                              "day*age old", "ice free days t-1",
                              "winter AO", "spring AO t-1"))

ggsave("07_results/01_interim_results/model_outputs/covariate_model_selection/graphs/incl_prob_reversible_jump_2.png",
       height = 4, width = 5)


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
                        "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]", "beta[11]",
                        "psi"),
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)

#______________________----
# ~ 3. With 8 predictors & no indicator --------------------------------------

# This is the other method to use the reversible jump method, i.e. without indicator

# ~~~ a. Write the model -------------------------------------------------------

glmNoIndicatorCode <- nimbleCode({

  for(i in 1:numVars) {
    beta[i] ~ dnorm(0, sd = 100)
  }
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta[1] +
      beta[2] * day_s[i] +
      beta[3] * age_young[i] +
      beta[4] * age_old[i] +
      beta[5] * size_s[i] +
      beta[6] * size_s2[i] +
      beta[7] * day_s[i] * age_young[i] +
      beta[8] * day_s[i] * age_old[i]
  }
})



# ~~~ b. Constants and data ----------------------------------------------------

## Set up the model.
constants <- list(N = length(y), 
                  numVars = 8)
dat  <- list(y = y, 
             day_s = day_s,
             age_young = age_young,
             age_old = age_old,
             size_s = size_s,
             size_s2 = size_s2)

# Inits
logit_model <- glm(formula = y ~ day_s * (age_young + age_old) 
                   + size_s + size_s2,
                   family = binomial(link = "logit"))

logit_model$coefficients
summary(logit_model)

inits <- list(beta = logit_model$coefficients)

glmNoIndicatorModel <- nimbleModel(code = glmNoIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)


# ~~~ c. Configuring RJMCMC  ---------------------------------------------------
glmNoIndicatorConf <- configureMCMC(glmNoIndicatorModel)
configureRJ(glmNoIndicatorConf,
            targetNodes = 'beta',
            priorProb = 0.5,
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
glmNoIndicatorConf$printSamplers(c("beta[1]"))


# ~~~ d. Build and run the RJMCMC ----------------------------------------------

mcmcNoIndicatorRJ <- buildMCMC(glmNoIndicatorConf)

cNoIndicatorModel <- compileNimble(glmNoIndicatorModel)

CMCMCNoIndicatorRJ <- compileNimble(mcmcNoIndicatorRJ, 
                                  project = glmNoIndicatorModel)

set.seed(1)
system.time(reversible_jump <- runMCMC(CMCMCNoIndicatorRJ, 
                                       niter = 150000, 
                                       nburnin = 10000,
                                       nchains = 2))



# ~~~ e. Individual inclusion probabilities ------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])


betaCols <- grep("beta\\[", colnames(allchains))
posterior_inclusion_proportions <- colMeans(apply(allchains[, betaCols],
                                                  2, function(x) x != 0))
posterior_inclusion_proportions <- data.frame(inclusion_prob = colMeans(apply(X = allchains[, betaCols],
                                                                            MARGIN = 2, 
                                                                            FUN = function(x) x != 0)),
                                              parameter = c("intercept", "day number", "age young", "age old",
                                                            "size", "size2", "day*age young",
                                                            "day*age old"))
                                                             

ggplot(data = posterior_inclusion_proportions, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("intercept", "day number", "age young", "age old",
                              "size", "size2", "day*age young",
                              "day*age old"))


# ~~~ f. Parameter distributions -----------------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(reversible_jump[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(reversible_jump[["chain1"]])[1]),
                   rep(2, times = dim(reversible_jump[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
                        "beta[6]", "beta[7]", "beta[8]"),
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump_3_(no_indicator).png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)






#______________________----
# ~ 2. With 6 predictors -------------------------------------------------------





