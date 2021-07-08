#==============================================================================#
#                                                                              #
#            Variable selection using the reversible jump MCMC method          #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")

# Load the environmental data --------------------------------------------------

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


# ~ 1. With 5 predictors -------------------------------------------------------

# ~~~ a. Generate the data ------------------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

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

true_betas <- c(-1.5, -2, -1.75, -1.75, 0)

z <- 1 - 1.5*day_number_s - 2*ice_free_days_previous_s - 
  1.75*age_young - 1.75*age_old + 0*size_s

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

logit_model <- glm(formula = y ~ day_number_s + ice_free_days_previous_s + 
                     age_young + age_old + size_s,
                   family = binomial(link = "logit"))
logit_model$coefficients
summary(logit_model)
# ~~~ b. Reversible jump with indicator variables -------------------------------------

# Next we set up the model. In this case we explicitly include indicator variables 
# that include or exclude the corresponding predictor variable. For this example 
# we assume the indicator variables are exchangeable and we include the inclusion 
# probability in the inference.

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0 +
      zbeta[1] * day_number_s[i] +
      zbeta[2] * ice_free_days_previous_s[i] +
      zbeta[3] * age_young[i] +
      zbeta[4] * age_old[i] +
      zbeta[5] * size_s[i]
    
  }
})

## Set up the model.
constants <- list(N = length(y), 
                  numVars = 5)
inits <- list(psi = 0.5,
              beta0 = logit_model$coefficients[1],
              beta = logit_model$coefficients[1:constants$numVars + 1],
              z = sample(0:1, constants$numVars, 0.5))



dat  <- list(y = y, 
             day_number_s = day_number_s,
             ice_free_days_previous_s = ice_free_days_previous_s,
             age_young = age_young,
             age_old = age_old,
             size_s = size_s)

glmIndicatorModel <- nimbleModel(code = glmIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)

# The above model code can potentially be used to set up variable selection in 
# NIMBLE without using RJMCMC, since the indicator variables can turn the regression 
# parameters off and on. However, in that case the MCMC sampling can be inefficient 
# because a given regression parameter can wander widely in the parameter space 
# when the corresponding variable is not in the model. This can make it difficult 
# for the variable to be put back into the model, unless the prior for that 
# variable is (perhaps artificially) made somewhat informative. Configuring 
# RJMCMC sampling via our NIMBLE function configureRJ results in the MCMC not 
# sampling the regression coefficients for variables for iterations where the 
# variables are not in the model.

# Configuring RJMCMC  ----------------------------------------------------------

# The RJMCMC sampler can be added to the MCMC configuration by calling the function 
# configureRJ(). In the example considered we introduced z as indicator variables 
# associated with the regression coefficients beta. We can pass these, respectively,
# to configureRJ using the arguments indicatorNodes and targetNodes. The control 
# arguments allow one to specify the mean and the scale of the normal proposal 
# distribution used when proposing to put a coefficient back into the model.

glmIndicatorConf <- configureMCMC(glmIndicatorModel)
glmIndicatorConf$addMonitors('z')

configureRJ(conf = glmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

# Checking the assigned samplers we see that the indicator variables are each 
# assigned an RJ_indicator sampler whose targetNode is the corresponding coefficient, 
# while the beta parameters have a RJ_toggled sampler. The latter sampler is a
# modified version of the original sampler to the targetNode that is invoked only 
# when the variable is currently in the model.

## Check the assigned samplers
glmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))
glmIndicatorConf$printSamplers(c("z[5]", "beta[5]"))

# Build and run the RJMCMC -----------------------------------------------------
mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, 
                                        niter = 50000, 
                                        nburnin = 5000,
                                        nchains = 2))



# ~~~ c. Looking at the results -------------------------------------------------------
# We can look at the sampled values of the indicator and corresponding coefficient 
# for some of the variables.

par(mfrow = c(2, 2))
plot(samplesIndicator[["chain1"]][,'beta[2]'], pch = 16, cex = 0.4, main = "beta[2] traceplot")
plot(samplesIndicator[["chain1"]][,'z[2]'], pch = 16, cex = 0.4, main = "z[2] traceplot")
plot(samplesIndicator[["chain1"]][,'beta[5]'], pch = 16, cex = 0.4, main = "beta[5] traceplot")
plot(samplesIndicator[["chain1"]][,'z[5]'], pch = 16, cex = 0.4, main = "z[5] traceplot")

# ~~~ d. Individual inclusion probabilities -------------------------------------------

# Now let’s look at the inference on the variable selection problem. We see that 
# the fourth and fifth predictors are almost always included (these are the ones 
# with the largest true coefficient values), while the others, including some 
# variables that are truly associated with the outcome but have smaller true 
# coefficient values, are almost never included.

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("day number", "ice free days",
                                                     "age young", "age old", "size"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("day number", "ice free days",
                              "age young", "age old", "size"))


# ~~~ e. Parameter distributions ------------------------------------------------------

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])
betaCols <- grep("beta\\[", colnames(allchains))

chains_l <- as.data.frame(allchains[, c(betaCols, max(betaCols) + 1)]) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "psi"),
               names_to = "parameter")

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump_dummy_data.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)


#______________________----
# ~ 2. With 6 predictors -------------------------------------------------------

# ~~~ a. Generate the data ------------------------------------------------------------

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

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

# Spring AO t-1 
prior_spring_AO <- data_model$prior_spring_AO
prior_spring_AO_s <- as.vector(scale(prior_spring_AO))

true_betas <- c(-1.5, -2, -1.75, -1.75, 0, 3)
z <- 1 - 1.5*day_number_s - 2*ice_free_days_previous_s - 
  1.75*age_young - 1.75*age_old + 0*size_s + 3*prior_spring_AO_s

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

logit_model <- glm(formula = y ~ day_number_s + ice_free_days_previous_s + 
                     age_young + age_old + size_s + prior_spring_AO_s,
                   family = binomial(link = "logit"))
logit_model$coefficients
summary(logit_model)


# ~~~ b. Reversible jump with indicator variables -------------------------------------

# Next we set up the model. In this case we explicitly include indicator variables 
# that include or exclude the corresponding predictor variable. For this example 
# we assume the indicator variables are exchangeable and we include the inclusion 
# probability in the inference.

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0 +
      zbeta[1] * day_number_s[i] +
      zbeta[2] * ice_free_days_previous_s[i] +
      zbeta[3] * age_young[i] +
      zbeta[4] * age_old[i] +
      zbeta[5] * size_s[i] +
      zbeta[6] * prior_spring_AO_s[i]
  }
})

z <- 1 - 1.5*day_number_s - 2*ice_free_days_previous_s - 
  1.75*age_young - 1.75*age_old + 0*size_s + 3*prior_spring_AO_s


## Set up the model.
constants <- list(N = length(y), 
                  numVars = 6)
inits <- list(psi = 0.5,
              beta0 = logit_model$coefficients[1],
              beta = logit_model$coefficients[1:constants$numVars + 1],
              z = sample(0:1, constants$numVars, 0.5))

dat  <- list(y = y, 
             day_number_s = day_number_s,
             ice_free_days_previous_s = ice_free_days_previous_s,
             age_young = age_young,
             age_old = age_old,
             prior_spring_AO_s = prior_spring_AO_s,
             size_s = size_s)


glmIndicatorModel <- nimbleModel(code = glmIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)

# The above model code can potentially be used to set up variable selection in 
# NIMBLE without using RJMCMC, since the indicator variables can turn the regression 
# parameters off and on. However, in that case the MCMC sampling can be inefficient 
# because a given regression parameter can wander widely in the parameter space 
# when the corresponding variable is not in the model. This can make it difficult 
# for the variable to be put back into the model, unless the prior for that 
# variable is (perhaps artificially) made somewhat informative. Configuring 
# RJMCMC sampling via our NIMBLE function configureRJ results in the MCMC not 
# sampling the regression coefficients for variables for iterations where the 
# variables are not in the model.

# Configuring RJMCMC  ----------------------------------------------------------

# The RJMCMC sampler can be added to the MCMC configuration by calling the function 
# configureRJ(). In the example considered we introduced z as indicator variables 
# associated with the regression coefficients beta. We can pass these, respectively,
# to configureRJ using the arguments indicatorNodes and targetNodes. The control 
# arguments allow one to specify the mean and the scale of the normal proposal 
# distribution used when proposing to put a coefficient back into the model.

glmIndicatorConf <- configureMCMC(glmIndicatorModel)
glmIndicatorConf$addMonitors('z')

configureRJ(conf = glmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

# Checking the assigned samplers we see that the indicator variables are each 
# assigned an RJ_indicator sampler whose targetNode is the corresponding coefficient, 
# while the beta parameters have a RJ_toggled sampler. The latter sampler is a
# modified version of the original sampler to the targetNode that is invoked only 
# when the variable is currently in the model.

## Check the assigned samplers
glmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))
glmIndicatorConf$printSamplers(c("z[5]", "beta[5]"))

# Build and run the RJMCMC -----------------------------------------------------
mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, 
                                        niter = 50000, 
                                        nburnin = 5000,
                                        nchains = 2))

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])

round(apply(X = allchains, MARGIN = 2, FUN = mean), 2)

# ~~~ c. Looking at the results -------------------------------------------------------
# We can look at the sampled values of the indicator and corresponding coefficient 
# for some of the variables.

par(mfrow = c(2, 2))
plot(samplesIndicator[["chain1"]][,'beta[2]'], pch = 16, cex = 0.4, main = "beta[2] traceplot")
plot(samplesIndicator[["chain1"]][,'z[2]'], pch = 16, cex = 0.4, main = "z[2] traceplot")
plot(samplesIndicator[["chain1"]][,'beta[5]'], pch = 16, cex = 0.4, main = "beta[5] traceplot")
plot(samplesIndicator[["chain1"]][,'z[5]'], pch = 16, cex = 0.4, main = "z[5] traceplot")

# ~~~ d. Individual inclusion probabilities -------------------------------------------

# Now let’s look at the inference on the variable selection problem. We see that 
# the fourth and fifth predictors are almost always included (these are the ones 
# with the largest true coefficient values), while the others, including some 
# variables that are truly associated with the outcome but have smaller true 
# coefficient values, are almost never included.

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("day number", "ice free days",
                                                     "age young", "age old", "size",
                                                     "prior spring AO"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("day number", "ice free days",
                              "age young", "age old", "size", "prior spring AO"),
                   labels = c("day number(-1.5)", "ice free days(-2)",
                              "age young(-1.75)", "age old(-1.75)", "size(0)",
                              "prior spring AO(3)")) +
  labs(x = "covariate(true beta)",
       y = "inclusion prbability")

ggsave(filename = 
         "07_results/01_interim_results/model_outputs/covariate_model_selection/graphs/reversible_jump_dummy_data_2_inclusion_prob.png",
       height = 4, width = 7)


# ~~~ e. Parameter distributions ------------------------------------------------------

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])
betaCols <- grep("beta\\[", colnames(allchains))

chains_l <- as.data.frame(allchains[, c(betaCols, max(betaCols) + 1)]) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "psi"),
               names_to = "parameter")

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump_dummy_data_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)



#______________________----
# ~ 3. With 6 pred + random effect-------------------------------------------------------

# ~~~ a. Generate the data ------------------------------------------------------------

# Date capture
day_number <- data_model$day_number
day_number_s <- as.vector(scale(day_number))

# Number of ice-free days t-1
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- as.vector(scale(ice_free_days_previous))

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

# Spring AO t-1 
prior_spring_AO <- data_model$prior_spring_AO
prior_spring_AO_s <- as.vector(scale(prior_spring_AO))

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


# random effect
random_temp <- rnorm(mean = 0, sd = 2, n = length(unique(year)))
year_list <- unique(data_model$year)
random <- c(rep(NA, times = nrow(data_model)))

for (k in 1:length(year_list)) {
  x <- which(data_model$year == year_list[k])
  random[x] <- random_temp[k]
}

true_betas <- c(-1.5, -2, -1.75, -1.75, 0, 3)
z <- 1 - 1.5*day_number_s - 2*ice_free_days_previous_s - 
  1.75*age_young - 1.75*age_old + 0*size_s + 3*prior_spring_AO_s + random

p <- plogis(z)              # inverse logit
y <- rbinom(244, 1, p)      # bernoulli response variable

logit_model <- glm(formula = y ~ day_number_s + ice_free_days_previous_s + 
                     age_young + age_old + size_s + prior_spring_AO_s,
                   family = binomial(link = "logit"))
logit_model$coefficients
summary(logit_model)


# ~~~ b. Reversible jump with indicator variables -------------------------------------

# Next we set up the model. In this case we explicitly include indicator variables 
# that include or exclude the corresponding predictor variable. For this example 
# we assume the indicator variables are exchangeable and we include the inclusion 
# probability in the inference.

glmIndicatorCode <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 1.5)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:N) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0 +
      zbeta[1] * day_number_s[i] +
      zbeta[2] * ice_free_days_previous_s[i] +
      zbeta[3] * age_young[i] +
      zbeta[4] * age_old[i] +
      zbeta[5] * size_s[i] +
      zbeta[6] * prior_spring_AO_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
})

## Set up the model.
constants <- list(N = length(y), 
                  numVars = 6,
                  year = as.numeric(year))

inits <- list(psi = 0.5,
              sigma = runif(1),
              beta0 = logit_model$coefficients[1],
              beta = logit_model$coefficients[1:constants$numVars + 1],
              z = sample(0:1, constants$numVars, 0.5))



dat  <- list(y = y, 
             nbyear = length(levels(year)),
             day_number_s = day_number_s,
             ice_free_days_previous_s = ice_free_days_previous_s,
             age_young = age_young,
             age_old = age_old,
             prior_spring_AO_s = prior_spring_AO_s,
             size_s = size_s)

glmIndicatorModel <- nimbleModel(code = glmIndicatorCode, 
                                 constants = constants,
                                 inits = inits, 
                                 data = dat)

# The above model code can potentially be used to set up variable selection in 
# NIMBLE without using RJMCMC, since the indicator variables can turn the regression 
# parameters off and on. However, in that case the MCMC sampling can be inefficient 
# because a given regression parameter can wander widely in the parameter space 
# when the corresponding variable is not in the model. This can make it difficult 
# for the variable to be put back into the model, unless the prior for that 
# variable is (perhaps artificially) made somewhat informative. Configuring 
# RJMCMC sampling via our NIMBLE function configureRJ results in the MCMC not 
# sampling the regression coefficients for variables for iterations where the 
# variables are not in the model.

# Configuring RJMCMC  ----------------------------------------------------------

# The RJMCMC sampler can be added to the MCMC configuration by calling the function 
# configureRJ(). In the example considered we introduced z as indicator variables 
# associated with the regression coefficients beta. We can pass these, respectively,
# to configureRJ using the arguments indicatorNodes and targetNodes. The control 
# arguments allow one to specify the mean and the scale of the normal proposal 
# distribution used when proposing to put a coefficient back into the model.

glmIndicatorConf <- configureMCMC(glmIndicatorModel)
glmIndicatorConf$addMonitors('z')

configureRJ(conf = glmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

# Checking the assigned samplers we see that the indicator variables are each 
# assigned an RJ_indicator sampler whose targetNode is the corresponding coefficient, 
# while the beta parameters have a RJ_toggled sampler. The latter sampler is a
# modified version of the original sampler to the targetNode that is invoked only 
# when the variable is currently in the model.

## Check the assigned samplers
glmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))
glmIndicatorConf$printSamplers(c("z[5]", "beta[5]"))

# Build and run the RJMCMC -----------------------------------------------------
mcmcIndicatorRJ <- buildMCMC(glmIndicatorConf)

cIndicatorModel <- compileNimble(glmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, 
                                  project = glmIndicatorModel)

set.seed(1)
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, 
                                        niter = 50000, 
                                        nburnin = 5000,
                                        nchains = 2))



# ~~~ c. Looking at the results -------------------------------------------------------
# We can look at the sampled values of the indicator and corresponding coefficient 
# for some of the variables.

par(mfrow = c(2, 2))
plot(samplesIndicator[["chain1"]][,'beta[2]'], pch = 16, cex = 0.4, main = "beta[2] traceplot")
plot(samplesIndicator[["chain1"]][,'z[2]'], pch = 16, cex = 0.4, main = "z[2] traceplot")
plot(samplesIndicator[["chain1"]][,'beta[5]'], pch = 16, cex = 0.4, main = "beta[5] traceplot")
plot(samplesIndicator[["chain1"]][,'z[5]'], pch = 16, cex = 0.4, main = "z[5] traceplot")

# ~~~ d. Individual inclusion probabilities -------------------------------------------

# Now let’s look at the inference on the variable selection problem. We see that 
# the fourth and fifth predictors are almost always included (these are the ones 
# with the largest true coefficient values), while the others, including some 
# variables that are truly associated with the outcome but have smaller true 
# coefficient values, are almost never included.

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = )


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("day number", "ice free days",
                              "age young", "age old", "size", "prior spring AO"),
                   labels = c("day number(-1.5)", "ice free days(-2)",
                              "age young(-1.75)", "age old(-1.75)", "size(0)",
                              "prior spring AO(3)")) +
  labs(x = "covariate(true beta)",
       y = "inclusion prbability")

ggsave(filename = 
         "07_results/01_interim_results/model_outputs/covariate_model_selection/graphs/reversible_jump_dummy_data_2_inclusion_prob.png",
       height = 4, width = 7)


# ~~~ e. Parameter distributions ------------------------------------------------------

allchains <- rbind(samplesIndicator[["chain1"]],
                   samplesIndicator[["chain2"]])
betaCols <- grep("beta\\[", colnames(allchains))

chains_l <- as.data.frame(allchains[, c(betaCols, max(betaCols) + 1)]) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
  mutate(chain = as.factor(chain)) %>% 
  pivot_longer(cols = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "psi"),
               names_to = "parameter")

chains_l <- as.data.frame(allchains) %>%
  mutate(iteration = rep(seq(1:dim(samplesIndicator[["chain1"]])[1]), 
                         times = 2),
         chain = c(rep(1, times = dim(samplesIndicator[["chain1"]])[1]),
                   rep(2, times = dim(samplesIndicator[["chain2"]])[1]))) %>%
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
            "07_results/01_interim_results/model_outputs/covariate_model_selection/diagnostic_plots/reversible_jump_dummy_data_2.png",
          plot = diagnostic_plot,
          ncol = 3,
          nrow = 1,
          base_width = 3,
          base_height = nrows * 1.1)
