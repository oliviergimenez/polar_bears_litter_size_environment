#------------------------------------------------------------------------------#
#
#                          Understand Thierry's code
#
#------------------------------------------------------------------------------#

library(tidyverse)
library(mlogit)
library(viridis)
library(ggmcmc)
library(gridExtra)
library(nimble)
Sys.setenv(LANG = "en")

# Code O'Hara et Sillanpaa -----------------------------------------------------

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



model {
  for(i in 1:N) {
    for(pp in 1:MarkerN)   {    # This is a long-handed way of doing inprod(), but is quicker
      SumMark[pp, i] <- X[i,10*pp-9]*beta[10*pp-9] + 
        X[i,10*pp-8]*beta[10*pp-8] +
        X[i,10*pp-7]*beta[10*pp-7] + 
        X[i,10*pp-6]*beta[10*pp-6] +
        X[i,10*pp-5]*beta[10*pp-5] + X[i,10*pp-4]*beta[10*pp-4] +
        X[i,10*pp-3]*beta[10*pp-3] + X[i,10*pp-2]*beta[10*pp-2] +
        X[i,10*pp-1]*beta[10*pp-1] + X[i,10*pp]*beta[10*pp]
    }
    mu[i] <- alpha + sum(SumMark[, i]) + X[i,10*MarkerN+1]*beta[10*MarkerN+1] +
      X[i,10*MarkerN+2]*beta[10*MarkerN+2] +X[i,10*MarkerN+3]*beta[10*MarkerN+3] +
      X[i,10*MarkerN+4]*beta[10*MarkerN+4] +X[i,10*MarkerN+5]*beta[10*MarkerN+5] +
      X[i,10*MarkerN+6]*beta[10*MarkerN+6] +X[i,10*MarkerN+7]*beta[10*MarkerN+7]
    y[i] ~ dnorm(mu[i], tau)             # Likelihood
    for(j in 1:NMarkers)   {         X[i,j] ~ dbern(0.5)      }    # Missing data
  }
  for(j in 1:NMarkers)   {
    beta [j] ~ dnorm(0,1.0E-6)                # Regression coefficient
    Ind [j] <- step(abs(beta[j]) - 0.05)       # Indicator
  }
  alpha ~ dnorm(0,1.0E-6)                   # Intercept
  tau ~ dgamma(1.0E-4,1.0E-4)             # Residual precision
}
Initial values
list(alpha=58.3)   list(alpha=57)   list(alpha=58)   list(alpha=59)   list(alpha=60)
list(tau=1, beta=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )


# Other code -------------------------------------------------------------------

list(x = c(1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,
           4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
           5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
           6,6,6,6,6,6,6,6,6,6,6,6,6,6,
           7,7,7,7,7,7,7,7,7,8,8,8,8,8))

list(T = 25,
     Model = 0)

model <- nimbleCode({
  for (i in 1:100){
    x[i] ~ dbin(p[i],10)
    logit(p[i]) <- theta[i]
    theta[i] ~ dnorm(mu,tau.theta)
  }
  Model ~ dbern(0.15)
  T ~ dunif(0,1000)
  tau.theta <- T*(1-Model) + Model*1.0E3
  tau.mu <- T/(1+Model)
  mu ~ dnorm(0.0,tau.mu)
  sd <- 1/sqrt(T)
})



start <- Sys.time()
assign(x = paste0("fit_", model_code),
       value = nimbleMCMC(code = get(paste0("model_", model_code)),     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits_null,          
                          monitors = params,   # parameters to monitor
                          niter = 20000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE))
end <- Sys.time()
end - start


# mean			sd		MC_error	val2.5pc		median	val97.5pc	start	sample
# Model		0.1722	0.3775	0.005287	0.0					0.0			1.0				1		422000
# 
# 
# BF_{01} = (Posterior odds)/(Prior Odds) =  (0.1722/0.8278) /(0.15/0.85) = 1.18 in 
# favor of Random Effects model 0 (note that sample variance only 6% larger than nominal under homogeneity).




# Thierry's code ---------------------------------------------------------------


# Generate data

model{				
  # ---------------------------------------------------------------------------#
  ### MODEL ###
  k ~ dcat(p.k[])
  
  for (i in 1:nind){
    y[i] ~ dbern(pi[i])
    
    logit(pi[i]) <- beta0 +
      beta.coh*coh.size[i] + 
      w.alpha*beta.alpha*alpha[i,k]
  } # 'i'
  
  # ------------#
  #   PRIORS	  #
  # ------------#
  
  beta0 ~ dnorm(0, tau.mu)					
  beta.coh ~ dnorm(0, tau.beta.coh)		
  beta.alpha ~ dnorm(0, tau.beta.alpha)	
  
  V ~ dunif(0,100)
  tau.V <- pow(V,-1)    # tau.V = 1/V
  
  ## Equations ##
  tau.mu <- tau.V/num.par.wt[M,1]
  tau.beta.coh <- tau.V/num.par.wt[M,2]
  tau.beta.alpha <- tau.V/num.par.wt[M,3]
  
  
  ########################################################
  # Models and associated parameter inclusion indicators #
  ########################################################
  
  M ~ dcat(p.M[])
  
  w.alpha <- equals(M, 2) # w.alpha <- 1 if M = 2, otherwise, w.alpha = 0
  
  for (i in 1:dim.M){
    Indicator.Model[i] <- equals(M,i)
  }
  
  
  # ---------------------#
  #  Derived Parameters  #
  # ---------------------#
  
  mu.pi.orig <- 1/(1+exp(-mu.pi))
  
  # 1.COHORT EFFECT#
  l.Coh.0 <- mu.pi							
  l.Coh.1 <- mu.pi + beta.coh		  			
  ## ORIGINAL SCALE : Logit BACK transformation (expit)
  Coh.0 <- 1/(1+exp(-l.Coh.0))
  Coh.1 <- 1/(1+exp(-l.Coh.1))
  beta.coh.orig <- Coh.1 - Coh.0
  
  # 2.MOM ALPHA EFFECT#				
  l.Al.0 <- mu.pi							
  l.Al.1 <- mu.pi + beta.alpha * 0.6794	  		
  ## ORIGINAL SCALE : Logit BACK transformation (expit)
  Al.0 <- 1/(1+exp(-l.Al.0))
  Al.1 <- 1/(1+exp(-l.Al.1))
  beta.alpha.orig.1SD <- Al.1 - Al.0
  
  
  # -------------------------------------------------------------------------- #
  ### MODEL ENDS ###
} 			

















