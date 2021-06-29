#==============================================================================#
#                                                                              #
#                        Nimble models for ice-free days                       #
#                                                                              #
#==============================================================================#


# A. March-May ice t-1 ==========================================================
# ~ 1. March-May ice t-1: effect 1c_VS_0c (1.4.2_D) ============================

model_1.4.2_D_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_previous_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 2. March-May ice t-1: effect 2-3c_VS_0c (1.4.2_D) =============================

model_1.4.2_D_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*mar_may_ice_previous_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 3. March-May ice t-1: effect common (1.4.2_D) =============================

model_1.4.2_D_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_previous_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*mar_may_ice_previous_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 4. March-May ice t-1: effect distinct (1.4.2_D) =============================

model_1.4.2_D_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_previous_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*mar_may_ice_previous_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})






# ~ 5. March-May ice t-1: binomial =============================================

model_1.4.2_D_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1 * mar_may_ice_previous_s[i] + eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})





# B. March-May ice t-2 ==========================================================
# ~ 1. March-May ice t-2: effect 1c_VS_0c (1.4.3_D) ==============================

model_1.4.3_D_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 2. March-May ice t-2: effect 2-3c_VS_0c (1.4.3_D) =============================

model_1.4.3_D_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 3. March-May ice t-2: effect common (1.4.3_D) =============================

model_1.4.3_D_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 4. March-May ice t-2: effect distinct (1.4.3_D) =============================

model_1.4.3_D_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*mar_may_ice_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})


#