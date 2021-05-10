#==============================================================================#
#                                                                              #
#                         Jags models for Winter AO                           #
#                                                                              #
#==============================================================================#



# 1. Ice-free days previous year: common slope (1.1.2_D_common) =============================

model_1.1.2_D_common_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
}



# 2. Ice-free days previous year: distinct slopes  (1.1.2_D_distinct) ========================

model_1.1.2_D_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
  
}



# 3. Ice-free days 2y-prior: common slope (1.1.3_D_common) ==========================

model_1.1.3_D_common_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
}




# 4. Ice-free days 2y-prior: distinct slopes  (1.1.3_D_distinct) =====================

model_1.1.3_D_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  b1 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
  c1 ~ dnorm(0.00000E+00, 0.1)
  
}



# 5. Ice-free days previous year w/o females w/ no cubs: common slope (1.1.2_D_common_bis) =============================

model_1.1.2_D_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
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



# model_1.1.2_D_common_slope_bis <- function()  {
#   
#   for (i in 1:N) {
#     y[i] ~ dcat(p[i, 1:J])
#     log(q[i, 1]) <- 0
#     log(q[i, 2]) <- b0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
#     log(q[i, 3]) <- c0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
#     
#     for (j in 1:J) {
#       p[i, j] <- q[i, j]/sum(q[i, 1:J])
#     }
#   }
#   
#   for (i in 1:nbyear) {
#     eps1[i] <- sigma1 * eta1[i]
#     eta1[i] ~ dnorm(0.00000E+00, 1)
#   }
#   
#   sigma1 ~ dunif(0.00000E+00, 10)
#   
#   b0 ~ dnorm(0.00000E+00, 0.1)
#   b1 ~ dnorm(0.00000E+00, 0.1)
#   c0 ~ dnorm(0.00000E+00, 0.1)
# }



# 6. Ice-free days previous year w/o females w/ no cubs: distinct slopes (1.1.2_D_distinct_bis) ========================

model_1.1.2_D_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
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



# 7. Ice-free days 2y-prior w/o females w/ no cubs: common slope (1.1.3_D_common_bis) ==========================

model_1.1.3_D_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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




# 8. Ice-free days 2y-prior w/o females w/ no cubs: distinct slopes  (1.1.3_D_distinct_bis) =====================

model_1.1.3_D_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*ice_free_days_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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



# 9. NIMBLE test, Ice-free days previous year: common slope (1.1.2_D_common) =============================
library(nimble)
Nimble_ie_free_days_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*ice_free_days_previous_s[i] + sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  a0 ~ dnorm(0.00000E+00, 0.1)
  a1 ~ dnorm(0.00000E+00, 0.1)
  b0 ~ dnorm(0.00000E+00, 0.1)
  c0 ~ dnorm(0.00000E+00, 0.1)
})




