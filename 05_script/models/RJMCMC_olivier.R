#==============================================================================#
#                                                                              #
#            Variable selection using the reversible jump MCMC method          #
#                             (multinomial)                                    #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(mlogit)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")


data_model <- read_csv("06_processed_data/CR_data/CR_f_clean.csv")

# Response variable. Cub number = 0, 1 or 4 (i.e 2 or 3)
y <- data_model %>%
  pull(cub_number_2)

# Date of capture
day <- as.numeric(data_model$day_number)
day_s <- as.vector(scale(day))


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

# Inits, use the output of the glm function
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
                                       niter = 50000, 
                                       nburnin = 10000,
                                       nchains = 2))



# ~~~ e. Individual inclusion probabilities ------------------------------------

allchains <- rbind(reversible_jump[["chain1"]],
                   reversible_jump[["chain2"]])

zCols <- grep("z\\[", colnames(allchains))
posterior_inclusion_prob <- data.frame(inclusion_prob = colMeans(allchains[, zCols]),
                                       parameter = c("intercept", "day number"))


ggplot(data = posterior_inclusion_prob, aes(x = parameter, y = inclusion_prob)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = c("intercept", "day number"))


# ~~~ f. Parameter distributions, mixing, convergence --------------------------

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
(diagnostic_plot <- plot_grid(trace.plots,
                             density.plots, 
                             running.mean.plot,
                             ncol = 3, nrow = 1)
)
