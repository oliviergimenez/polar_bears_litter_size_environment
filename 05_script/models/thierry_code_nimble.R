#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#                       Try Thierry"s code with Nimble                         #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(nimble)
library(cowplot)
Sys.setenv(LANG = "en")


# A. Thierry's code as it is ---------------------------------------------------
# ~ 1. Load and process the data -----------------------------------------------
data.location <- paste(getwd(),"/03_methods_protocols/Model selection Bayesian Kuo & Mallick/Model-Recr-v2/1.orig.data", sep="")


data <- read.csv(file = paste0(data.location, "/tbl-data.csv"), sep = ";")

# Number of rows
nind <- dim(data)[1]

# #-----------------------------------------------------------------------#
# ### COVARIATE: MomAge
# #-------------------------------------------#
# # GET THE MomAge Values - NON-STANDARDIZED
# MomAge.ns <- data$MomAge
# 
# # STANDARDIZE THE VALUE of MomAge
# min.age <- min(MomAge.ns, na.rm = TRUE)
# max.age <- max(MomAge.ns, na.rm = TRUE)
# dist.age <- min.age:max.age
# mean.age <- mean(dist.age)
# sd.age <- sd(dist.age)
# 
# MomAge <- (MomAge.ns - mean.age) / sd.age
# head(MomAge); length(MomAge)


#-----------------------------------------------------------------------#
### COVARIATE: Cohort Size (= nb of pups for each year/cohort)
#-------------------------------------------#
# Get the Table with Cohort Size - NON-STANDARDIZED
APC <- read.csv( file = paste(data.location, "/APC.csv", sep = ""), sep = ";")

coh.size.ns <- NA
for (i in 1:dim(data)[1]){
  coh.size.ns[i] <- APC[which(APC$season == data$ChildCohort[i]), 2]
}
head(coh.size.ns); length(coh.size.ns)

# STANDARDIZE THE VALUE of Cohort Size
mean.coh <- mean(coh.size.ns)
sd.coh <- sd(coh.size.ns)
coh.size <- (coh.size.ns - mean.coh) / sd.coh
head(coh.size) ; mean(coh.size) ; sd(coh.size) ; length(coh.size)


# #-----------------------------------------------------------------------#
# ### COVARIATE: MomExp standardized by age  
# #   =>   Beta.exp represents the "cost" of having a first-time mom
# #-------------------------------------------#
# 
# MomExp <- data$MomExpSTD
# head(MomExp) ; length(MomExp)

#-----------------------------------------------------------------------#




### 	ALPHA EMPIRICAL DISTRIBUTION

#### MATRIX OF alpha"s POSTERIOR DISTRIBUTION SAMPLES
load(file = paste(data.location, "/alphas.RData", sep = ""))
alphas[1:10,1:10] ; dim(alphas)

spl <- seq(from = 1, to = dim(alphas)[1], by = 10)
head(spl) ; length(spl) ; max(spl)

alpha.matrix <- as.matrix(t(alphas[spl,]))
rownames(alpha.matrix) <- NULL
alpha.matrix[1:10,1:10]
dim(alpha.matrix)
rm(alphas)

alpha <- as.matrix(alpha.matrix[data$MomAlphaID,])
rownames(alpha) <- NULL
alpha[1:20,1:10] ; dim(alpha)
rm(alpha.matrix)

### Equal Proba for each Posterior draw of alpha
K.tot <- dim(alpha)[2]
p.k <- rep(1/K.tot, times = K.tot)
head(p.k)
length(p.k)




### ********************************************************** ###
y <- data$Recr
head(y) ; length(y)


### ********************************************************** ###
# Number of Parameter weight for each Model
num.par.wt <- as.matrix(read.csv("03_methods_protocols/Model selection Bayesian Kuo & Mallick/Model-Recr-v2/Model-Matrix.csv", sep = ";", header = F))
num.par.wt

dim.M <- dim(num.par.wt)[1] 
p.M <- rep(1/dim.M, dim.M)


### ********************************************************** ###
# Bundle data
my.constants <- list(nind = nind,
                     coh.size = coh.size,
                     p.k = p.k, 
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            alpha = alpha,
            num.par.wt = num.par.wt)

### **********************************************************




# ~ 1bis. Reduce the sample size -----------------------------------------------

nind.bis <- 300

set.seed(seed = 1111)
index <- sample(1:nind, size = nind.bis,
               replace = FALSE)

y.bis <- y[index]
coh.size.bis <- coh.size[index]
p.k.bis <- p.k[1:500]

alpha.bis <- alpha[index, 1:500]
K.tot <- dim(alpha.bis)[2]

### ********************************************************** ###
# Bundle data
my.constants <- list(nind = nind.bis,
                     coh.size = coh.size.bis,
                     p.k = p.k.bis, 
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = y.bis,
            alpha = alpha.bis,
            num.par.wt = num.par.wt)

### **********************************************************



# ~ 2. Code of the model -------------------------------------------------------

model_test_recr <- nimbleCode({
  k ~ dcat(p.k[])
  
  for (i in 1:nind){
    logit(pi[i]) <- mu.pi + 
      beta.coh*coh.size[i] + 
      w.alpha*beta.alpha*alpha[i,k]
  } # "i"
  
  for (i in 1:nind){
    y[i] ~ dbern(pi[i])
  } # "i"
  
  # ------------#
  #   PRIORS	  #
  # ------------#
  
  mu.pi ~ dnorm(0, V.mu)					
  beta.coh ~ dnorm(0, V.beta.coh)		
  beta.alpha ~ dnorm(0, V.beta.alpha)	
  
  V ~ dunif(0,100)
  
  ## Equations ##
  V.mu <- V*num.par.wt[M,1]
  V.beta.coh <- V*num.par.wt[M,2]
  V.beta.alpha <- V*num.par.wt[M,3]
  
  
  #######################################################+
  # Models and associated parameter inclusion indicators #
  #######################################################+
  
  M ~ dcat(p.M[])
  
  w.alpha <- equals(M, 2)
  
  for (i in 1:dim.M){
    Indicator.Model[i] <- equals(M, i)
  }
  
})

# ~ 3.  initial values ---------------------------------------------------------
inits <- function() list(mu.pi = rnorm(1, 0, 10), 
                         beta.alpha = rnorm(1, 0, 10),
                         beta.coh = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1),
                         k = sample(K.tot, size = 1)) 

# ~ 4.  Parameters to monitor --------------------------------------------------
params <- c("mu.pi", "beta.coh", "beta.alpha", "V", "Indicator.Model", "w.alpha",
            "k", "M") 


# ~ 5. Run the model -----------------------------------------------------------

start <- Sys.time()
fit_test_recr <-  nimbleMCMC(code = model_test_recr,     # model code  
                               data = dat,                                   
                               constants = my.constants,        
                               inits = inits,          
                               monitors = params,   # parameters to monitor
                               niter = 100000,                  # nb iterations
                               nburnin = 10000,              # length of the burn-in
                               nchains = 2,
                               summary = TRUE,
                               WAIC = TRUE)
end <- Sys.time()
end - start

save(fit_test_recr, file = "07_results/01_interim_results/model_outputs/model_recr_model_selection.RData")

load("07_results/01_interim_results/model_outputs/model_recr_model_selection.RData")

fit_test_recr$summary



# ~6. Check convergence --------------------------------------------------------

params.plot <- c("mu.pi", "beta.coh", "beta.alpha", "V", 
                 "w.alpha",
                 "k", "M") 

chain1 <- data.frame(fit_test_recr[["samples"]][["chain1"]]) %>%
  dplyr::select(params.plot) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_test_recr[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_test_recr[["samples"]][["chain2"]]) %>%
  dplyr::select(params.plot) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_test_recr[["samples"]][["chain2"]])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

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
  geom_density(alpha = 0.25) +
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
diagnostic_plot




#___________________________________________________________________________----
# B. Thierry's code with only one alpha col ------------------------------------

# ~ 1. Load and process the data -----------------------------------------------
data.location <- paste(getwd(),"/03_methods_protocols/Model selection Bayesian Kuo & Mallick/Model-Recr-v2/1.orig.data", sep="")

data <- read.csv(file = paste0(data.location, "/tbl-data.csv"), sep = ";")

# Number of rows
nind <- dim(data)[1]

# #-----------------------------------------------------------------------#
# ### COVARIATE: MomAge
# #-------------------------------------------#
# # GET THE MomAge Values - NON-STANDARDIZED
# MomAge.ns <- data$MomAge
# 
# # STANDARDIZE THE VALUE of MomAge
# min.age <- min(MomAge.ns, na.rm = TRUE)
# max.age <- max(MomAge.ns, na.rm = TRUE)
# dist.age <- min.age:max.age
# mean.age <- mean(dist.age)
# sd.age <- sd(dist.age)
# 
# MomAge <- (MomAge.ns - mean.age) / sd.age
# head(MomAge); length(MomAge)


#-----------------------------------------------------------------------#
### COVARIATE: Cohort Size (= nb of pups for each year/cohort)
#-------------------------------------------#
# Get the Table with Cohort Size - NON-STANDARDIZED
APC <- read.csv( file = paste(data.location, "/APC.csv", sep = ""), sep = ";")

coh.size.ns <- NA
for (i in 1:dim(data)[1]){
  coh.size.ns[i] <- APC[which(APC$season == data$ChildCohort[i]), 2]
}
head(coh.size.ns); length(coh.size.ns)

# STANDARDIZE THE VALUE of Cohort Size
mean.coh <- mean(coh.size.ns)
sd.coh <- sd(coh.size.ns)
coh.size <- (coh.size.ns - mean.coh) / sd.coh
head(coh.size) ; mean(coh.size) ; sd(coh.size) ; length(coh.size)


# #-----------------------------------------------------------------------#
# ### COVARIATE: MomExp standardized by age  
# #   =>   Beta.exp represents the "cost" of having a first-time mom
# #-------------------------------------------#
# 
# MomExp <- data$MomExpSTD
# head(MomExp) ; length(MomExp)

#-----------------------------------------------------------------------#

### 	ALPHA EMPIRICAL DISTRIBUTION

#### MATRIX OF alpha"s POSTERIOR DISTRIBUTION SAMPLES
load(file = paste(data.location, "/alphas.RData", sep = ""))
alphas[1:10,1:10] ; dim(alphas)

spl <- seq(from = 1, to = dim(alphas)[1], by = 10)
head(spl) ; length(spl) ; max(spl)

alpha.matrix <- as.matrix(t(alphas[spl,]))
rownames(alpha.matrix) <- NULL
alpha.matrix[1:10,1:10]
dim(alpha.matrix)
rm(alphas)

alpha <- as.matrix(alpha.matrix[data$MomAlphaID,])
rownames(alpha) <- NULL
alpha[1:20,1:10] ; dim(alpha)
rm(alpha.matrix)




### ********************************************************** ###
y <- data$Recr
head(y) ; length(y)


### ********************************************************** ###
# Number of Parameter weight for each Model
num.par.wt <- as.matrix(read.csv("03_methods_protocols/Model selection Bayesian Kuo & Mallick/Model-Recr-v2/Model-Matrix.csv", sep = ";", header = F))
num.par.wt

dim.M <- dim(num.par.wt)[1] 
p.M <- rep(1/dim.M, dim.M)


### ********************************************************** ###
# Bundle data
my.constants <- list(nind = nind,
                     coh.size = coh.size,
                     alpha = alpha[, 1], 
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = as.numeric(y),
            alpha = alpha,
            num.par.wt = num.par.wt)

### **********************************************************




# ~ 1bis. Reduce the sample size -----------------------------------------------

nind.bis <- 300

set.seed(seed = 1111)
index <- sample(1:nind, size = nind.bis,
                replace = FALSE)

y.bis <- y[index]
coh.size.bis <- coh.size[index]

alpha.bis <- alpha[index, 1]

### ********************************************************** ###
# Bundle data
my.constants <- list(nind = nind.bis,
                     coh.size = coh.size.bis,
                     alpha = alpha.bis, 
                     dim.M = dim.M,  
                     p.M = p.M)

dat <- list(y = y.bis,
            num.par.wt = num.par.wt)

### **********************************************************



# ~ 2. Code of the model -------------------------------------------------------

model_test_recr_one_alpha_col <- nimbleCode({
  for (i in 1:nind){
    logit(pi[i]) <- mu.pi + 
      beta.coh*coh.size[i] + 
      w.alpha*beta.alpha*alpha[i]
  } # "i"
  
  for (i in 1:nind){
    y[i] ~ dbern(pi[i])
  } # "i"
  
  # ------------#
  #   PRIORS	  #
  # ------------#
  
  mu.pi ~ dnorm(0, V.mu)					
  beta.coh ~ dnorm(0, V.beta.coh)		
  beta.alpha ~ dnorm(0, V.beta.alpha)	
  
  V ~ dunif(0,100)
  
  ## Equations ##
  V.mu <- V*num.par.wt[M,1]
  V.beta.coh <- V*num.par.wt[M,2]
  V.beta.alpha <- V*num.par.wt[M,3]
  
  
  #######################################################+
  # Models and associated parameter inclusion indicators #
  #######################################################+
  
  M ~ dcat(p.M[])
  
  w.alpha <- equals(M, 2)
  
  for (i in 1:dim.M){
    Indicator.Model[i] <- equals(M, i)
  }
  
})

# ~ 3.  initial values ---------------------------------------------------------
inits <- function() list(mu.pi = rnorm(1, 0, 10), 
                         beta.alpha = rnorm(1, 0, 10),
                         beta.coh = rnorm(1, 0, 10),
                         V = runif(1, 1, 100),
                         M = sample(dim.M, size = 1))

# ~ 4.  Parameters to monitor --------------------------------------------------
params <- c("mu.pi", "beta.coh", "beta.alpha", "V", "Indicator.Model", "w.alpha",
            "M") 


# ~ 5. Run the model -----------------------------------------------------------

start <- Sys.time()
fit_test_recr_one_alpha_col <-  nimbleMCMC(code = model_test_recr_one_alpha_col,     # model code  
                                           data = dat,                                   
                                           constants = my.constants,        
                                           inits = inits,          
                                           monitors = params,   # parameters to monitor
                                           niter = 100000,                  # nb iterations
                                           nburnin = 10000,              # length of the burn-in
                                           nchains = 2,
                                           summary = TRUE,
                                           WAIC = TRUE)
end <- Sys.time()
end - start

save(fit_test_recr_one_alpha_col, file = "07_results/01_interim_results/model_outputs/model_recr_with_only_one_alpha_model_selection.RData")

load("07_results/01_interim_results/model_outputs/model_recr_with_only_one_alpha_model_selection.RData")

fit_test_recr_one_alpha_col$summary



# ~6. Check convergence --------------------------------------------------------

params.plot <- c("mu.pi", "beta.coh", "beta.alpha", "V", 
                 "w.alpha", "M") 

chain1 <- data.frame(fit_test_recr_one_alpha_col[["samples"]][["chain1"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "1",
         iteration = seq(1, dim(fit_test_recr_one_alpha_col[["samples"]][["chain1"]])[1], by = 1))
chain2 <- data.frame(fit_test_recr_one_alpha_col[["samples"]][["chain2"]]) %>%
  dplyr::select(all_of(params.plot)) %>%
  mutate(chain = "2",
         iteration = seq(1, dim(fit_test_recr_one_alpha_col[["samples"]][["chain2"]])[1], by = 1))
chains <- rbind(chain1, chain2) 

chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 

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
  geom_density(alpha = 0.25) +
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
diagnostic_plot
