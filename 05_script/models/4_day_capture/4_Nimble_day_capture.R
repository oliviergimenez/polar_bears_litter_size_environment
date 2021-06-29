#==============================================================================#
#                                                                              #
#                        Nimble models for ice-free days                       #
#                                                                              #
#==============================================================================#



# ~ 1. Day of capture: effect 1c_VS_0c ============================

model_4_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*day_s[i] + eps1[year[i]]
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


# ~ 2. Day of capture: effect 2-3c_VS_0c =======================================

model_4_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*day_s[i] + eps1[year[i]]
    
    
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

# ~ 3. Day of capture: effect common ===========================================

model_4_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*day_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*day_s[i] + eps1[year[i]]
    
    
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

# ~ 4. Day of capture: effect distinct =========================================

model_4_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*day_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*day_s[i] + eps1[year[i]]
    
    
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






# ~ 5. Day of capture: binomial =============================================

model_4_effect_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b0 + b1 * day_s[i] + eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})






