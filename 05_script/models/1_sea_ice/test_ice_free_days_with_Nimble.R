#==============================================================================#
#                                                                              #
#             Test model with ice free days only, using Nimble                 #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(mlogit)
library(R2jags)
library(viridis)
library(ggmcmc)
library(gridExtra)
library(nimble)
Sys.setenv(LANG = "en")

writeLines(
  'PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"',
  con = file("~/.Renviron", open = "a")
)


# A. Ready the data ------------------------------------------------------------

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
  filter(year != 1992) # Remove the captures from 1992 since I can't calculate the ice free days in 1991



# build dataset

# The variables will be stored in lists because JAGS requires lists
y <- factor(as.numeric(data_model$cub_number))
summary(y)

# Renumérotation des années
{year <- data_model$year
year <- factor(year) 
year2 <- NULL
for (i in 1:length(year)){
  year2 <- c(year2, which(year[i] == levels(year)))
}
year <- factor(year2)
nbyear <- length(levels(year))
year
}

# Predictor
ice_free_days_previous <- data_model$ice_free_days_previous
ice_free_days_previous_s <- (ice_free_days_previous - mean(ice_free_days_previous))/sd(ice_free_days_previous) 

N <- length(y) # nb of reproductive events
J <- length(levels(y)) # number of categories


# B. Models --------------------------------------------------------------------

# ~ 1. Model without random effect ---------------------------------------------

# Put all the data in a list
# dat <- list(y, N, J, year, nbyear, ice_free_days_previous_s)
# names(dat) <- c("y", "N", "J","year", "nbyear", "ice_free_days_previous_s")
# dat <- list(y, N, J, ice_free_days_previous_s)
dat <- list(as.numeric(y), as.numeric(ice_free_days_previous_s))

names(dat) <- c("y", "ice_free_days_previous_s")

# Define the parameters to estimate
params <- c("a0", "a1", "b0", "c0") #, "sigma1" ,"eps1")


# Generate starting values
temp.dat <- data.frame(y = y, var_scaled = ice_free_days_previous_s)
temp.dat$yfac <- as.factor(temp.dat$y)   # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| var_scaled, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)
# inits1 <- list(a0 = coefs[1], b0 = coefs[2], c0 = coefs[3],
#                a1 = coefs[4],
#                sigma1 = runif(1))
# inits2 <- list(a0 = coefs[1] + 0.1, b0 = coefs[2] - 0.1, c0 = coefs[3] + 0.1,
#                a1 = coefs[4] -0.1,
#                sigma1 = runif(1))


inits <- function() list(inits1 = list(a0 = coefs[1], b0 = coefs[2], c0 = coefs[3],
                                       a1 = coefs[4]),
                         inits2 = list(a0 = coefs[1] + 0.1, b0 = coefs[2] - 0.1, c0 = coefs[3] + 0.1,
                                       a1 = coefs[4] -0.1))


# Define constants
my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y))) 



# ~~~ a. Nimble model --------------------------------------------------------------

Nimble_ie_free_days_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] #+ sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_previous_s[i] #+ sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_previous_s[i] #+ sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  # for (i in 1:nbyear) {
  #   eps1[i] <- sigma1 * eta1[i]
  #   eta1[i] ~ dnorm(0.00000E+00, 1)
  # }
  # 
  # sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
})

writeLines(
  'PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"',
  con = file("~/.Renviron", open = "a")
)

hmm.phiflowREp <- nimbleCode({
  for (t in 1:(T-1)){
    logit(phi[t]) <- beta[1] + beta[2] * flow[t] + eps[t] 
    eps[t] ~ dnorm(0, sd = sdeps) 
    ...  
  }
  sdeps ~ dunif(0,10) 
  ...
})  
  

# ~~~ b. Run the model ---------------------------------------------------------

mcmc.output <- nimbleMCMC(code = Nimble_ie_free_days_common,     # model code  
                          data = dat,                  # data
                          constants = my.constants,        # constants
                          inits = inits,          # initial values
                          monitors = params,   # parameters to monitor
                          niter = 30000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2)

  

# 
processed_output <- ggs(as.mcmc(mcmc.output)) %>%
  filter(Parameter %in% c("a0", "a1", "a2", "b0", "b1", "c0", "c1", "deviance"))

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
       filename = paste0("07_results/01_interim_results/model_outputs/graphs/diagnostic_plots/model_",
                         model_code, "_", slope, toupper(mode), ".png"),
       width = 17.5, height = 2.5*nbr_rows) 



# ~ 2. Model with random effect ------------------------------------------------

# Put all the data in a list
# dat <- list(y, N, J, year, nbyear, ice_free_days_previous_s)
# names(dat) <- c("y", "N", "J","year", "nbyear", "ice_free_days_previous_s")
# dat <- list(y, N, J, ice_free_days_previous_s)

# dat <- list(as.numeric(y), as.numeric(year), as.numeric(nbyear), as.numeric(ice_free_days_previous_s))
# names(dat) <- c("y", "year", "nbyear", "ice_free_days_previous_s")

dat <- list(as.numeric(y), as.numeric(ice_free_days_previous_s))
names(dat) <- c("y", "ice_free_days_previous_s")

# Define the parameters to estimate
params <- c("a0", "a1", "b0", "c0", "sigma1", "eps1") #, "sigma1" ,"eps1")


# Generate starting values
temp.dat <- data.frame(y = y, var_scaled = ice_free_days_previous_s)
temp.dat$yfac <- as.factor(temp.dat$y)   # Ajout de Y en facteur
mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
mlogit.mod <- mlogit(yfac ~ 1| var_scaled, 
                     data = mnl.dat, 
                     reflevel = "0")

coefs <- as.vector(summary(mlogit.mod)$coefficients)

inits <- function() list(inits1 = list(a0 = coefs[1], b0 = coefs[2], c0 = coefs[3],
                                       a1 = coefs[4], sdeps = runif(1)),
                         inits2 = list(a0 = coefs[1] + 0.1, b0 = coefs[2] - 0.1, c0 = coefs[3] + 0.1,
                                       a1 = coefs[4] -0.1, sdeps = runif(1)))


# Define constants
my.constants <- list(N = length(y), # nb of females captured
                     J = length(levels(y)),
                     year = as.numeric(year),
                     nbyear = nbyear) 



# ~~~ a. Nimble model --------------------------------------------------------------

library(nimble)
Nimble_ice_free_days_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] + eps1[year[i]] # eps[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_previous_s[i] + eps1[year[i]] # eps[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_previous_s[i] + eps1[year[i]] # eps[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1) 
  }
   
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
})



# ~~~ b. Run the model ---------------------------------------------------------
start <- Sys.time()
mcmc.output <- nimbleMCMC(code = Nimble_ice_free_days_common,     # model code  
                          data = dat,                  # data
                          constants = my.constants,        # constants
                          inits = inits,          # initial values
                          monitors = params,   # parameters to monitor
                          niter = 25000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2)
end <- Sys.time()
end - start


# 
processed_output <- ggs(as.mcmc(mcmc.output)) %>%
  filter(Parameter %in% c("a0", "a1", "a2", "b0", "b1", "c0", "c1", "deviance"))

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
       filename = paste0("07_results/01_interim_results/model_outputs/graphs/diagnostic_plots/model_",
                         model_code, "_", slope, toupper(mode), ".png"),
       width = 17.5, height = 2.5*nbr_rows)



# ~~~ c. Plot the model --------------------------------------------------------

res <- rbind(mcmc.output$chain1, mcmc.output$chain2)

# If females without cubs are included

b1cub <- res[, c(1, 2)]
b2cub <- res[, c(3, 2)] 
b3cub <- res[, c(4, 2)] 


range <- range(ice_free_days_previous_s)

lengthgrid <- 100
grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                   to = range[2] + 0.1*(range[2] - range[1]), 
                   length = lengthgrid)

grid <- grid_scaled * sd(ice_free_days_previous) + mean(ice_free_days_previous)


q0cub <- q1cub <- q2cub <- q3cub <- matrix(data = NA, 
                                           nrow = dim(b2cub)[1], 
                                           ncol = lengthgrid)
for (i in 1:lengthgrid) {
  for (j in 1:dim(b2cub)[1]) {
    q0cub[j, i] <- exp(0)
    q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
    q2cub[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid_scaled[i])	
    q3cub[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid_scaled[i])		
  }
}
# Backtransform
p0cub <- p1cub <- p2cub <- p3cub <- matrix(NA, dim(b2cub)[1], lengthgrid)
for (i in 1:lengthgrid){
  for (j in 1:dim(b2cub)[1]){
    norm <- (q0cub[j, i] + q1cub[j, i] + q2cub[j, i] + q3cub[j,i])
    p0cub[j, i] <- q0cub[j, i]/norm
    p1cub[j, i] <- q1cub[j, i]/norm
    p2cub[j, i] <- q2cub[j, i]/norm
    p3cub[j, i] <- q3cub[j, i]/norm
  }
}
df.for.plot <- data.frame(var = grid,
                          mean_p_0_cub = apply(p0cub, 2, mean),
                          mean_p_1_cub = apply(p1cub, 2, mean),
                          mean_p_2_cub = apply(p2cub, 2, mean),
                          mean_p_3_cub = apply(p3cub, 2, mean),
                          ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                          ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                          ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                          ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                          ci_p_2_cub_2.5 = apply(p2cub, 2, quantile, probs = 0.025),
                          ci_p_2_cub_97.5 = apply(p2cub, 2, quantile, probs = 0.975),
                          ci_p_3_cub_2.5 = apply(p3cub, 2, quantile, probs = 0.025),
                          ci_p_3_cub_97.5 = apply(p3cub, 2, quantile, probs = 0.975)) %>%
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


ggplot(data = df.for.plot, 
       aes(x = var, y = value, group = name, linetype = type, color = as.factor(cub_number))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE,                       
                      labels = color_labels) +
  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = "Ice free days during the previous year",
       y = "Probability", 
       color = "",
       linetype = "") 
