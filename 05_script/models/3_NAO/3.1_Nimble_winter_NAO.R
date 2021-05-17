#==============================================================================#
#                                                                              #
#                        Nimble models for winter NAO                          #
#                                                                              #
#==============================================================================#

# A. Winter NAO t  ==============================================================

# ~ 1. Winter NAO t: effect 1c_VS_0c (3.1.1.1) ==================================

model_3.1.1.1_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_s[i] + eps1[year[i]]
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


# ~ 2. Winter NAO t: effect 2-3c_VS_0c (3.1.1.1) ================================

model_3.1.1.1_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 3. Winter NAO t: effect common (3.1.1.1) =============================

model_3.1.1.1_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 4. Winter NAO t: effect distinct (3.1.1.1) =================================

model_3.1.1.1_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})



# ~ 5. Winter NAO t FACTOR : effect 1c_VS_0c (3.1.1.2) ==========================

model_3.1.1.2_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_pos[i] + eps1[year[i]]
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


# ~ 6. Winter NAO t FACTOR: effect 2-3c_VS_0c (3.1.1.2) =========================

model_3.1.1.2_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 7. Winter NAO t FACTOR: effect common (3.1.1.2) =============================

model_3.1.1.2_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 8. Winter NAO t FACTOR: effect distinct (3.1.1.2) ===========================

model_3.1.1.2_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})


# B. Winter NAO t-1  ============================================================

# ~ 1. Winter NAO t-1: effect 1c_VS_0c (3.1.2.1) ====================================

model_3.1.2.1_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_s[i] + eps1[year[i]]
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


# ~ 2. Winter NAO t-1: effect 2-3c_VS_0c (3.1.2.1) ================================

model_3.1.2.1_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 3. Winter NAO t-1: effect common (3.1.2.1) =============================

model_3.1.2.1_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 4. Winter NAO t-1: effect distinct (3.1.2.1) =================================

model_3.1.2.1_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})



# ~ 5. Winter NAO t-1 FACTOR : effect 1c_VS_0c (3.1.2.2) ==========================

model_3.1.2.2_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_pos[i] + eps1[year[i]]
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


# ~ 6. Winter NAO t-1 FACTOR: effect 2-3c_VS_0c (3.1.2.2) =========================

model_3.1.2.2_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 7. Winter NAO t-1 FACTOR: effect common (3.1.2.2) =============================

model_3.1.2.2_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 8. Winter NAO t-1 FACTOR: effect distinct (3.1.2.2) ===========================

model_3.1.2.2_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_prior_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})


# C. Winter NAO t-2  ============================================================

# ~ 1. Winter NAO t-2: effect 1c_VS_0c (3.1.3.1) ====================================

model_3.1.3.1_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_s[i] + eps1[year[i]]
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


# ~ 2. Winter NAO t-2: effect 2-3c_VS_0c (3.1.3.1) ================================

model_3.1.3.1_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 3. Winter NAO t-2: effect common (3.1.3.1) =============================

model_3.1.3.1_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 4. Winter NAO t-2: effect distinct (3.1.3.1) =================================

model_3.1.3.1_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_s[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_2y_prior_s[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})



# ~ 5. Winter NAO t-2 FACTOR : effect 1c_VS_0c (3.1.3.2) ==========================

model_3.1.3.2_effect_1c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
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


# ~ 6. Winter NAO t-2 FACTOR: effect 2-3c_VS_0c (3.1.3.2) =========================

model_3.1.3.2_effect_2_3c_VS_0c <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
})


# ~ 7. Winter NAO t-2 FACTOR: effect common (3.1.3.2) =============================

model_3.1.3.2_effect_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
})

# ~ 8. Winter NAO t-2 FACTOR: effect distinct (3.1.3.2) ===========================

model_3.1.3.2_effect_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_NAO_2y_prior_pos[i] + eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  a0 ~ dnorm(0.00000E+00, sd = 1.5)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  a1 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
  
})









# 10. Winter NAO: common slope (3.1.1.1_common_bis) =================================

model_3.1.1.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_NAO_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
}



# 11. Winter NAO w/o females w/ no cubs: distinct slopes (3.1.1.1_distinct_bis) ==============================

model_3.1.1.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
  
}



# 12. Previous Winter NAO w/o females w/ no cubs: common slope (3.1.2.1_common_bis) =============================

model_3.1.2.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_NAO_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
}



# 13. Previous Winter NAO w/o females w/ no cubs: distinct slopes  (3.1.2.1_distinct_bis) ========================

model_3.1.2.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
  
}



# 14. 2y prior Winter NAO w/o females w/ no cubs: common slope (3.1.3.1_common_bis) ==========================

model_3.1.3.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_NAO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
}




# 15. 2y prior Winter NAO w/o females w/ no cubs: distinct slopes  (3.1.3.1_distinct_bis) =====================

model_3.1.3.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
  
}






# 16. Winter NAO FACTOR w/o females w/ no cubs: distinct slopes (3.1.1.2_distinct_bis_bis) ====

model_3.1.1.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_pos[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
}


# 17. Prior Winter NAO FACTOR w/o females w/ no cubs: distinct slopes (3.1.2.2_distinct_bis) ====

model_3.1.2.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_prior_pos[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
}


# 18. 2y prior Winter NAO FACTOR w/o females w/ no cubs: distinct slopes (3.1.3.2_distinct_bis) ====

model_3.1.3.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_NAO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_NAO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
}



