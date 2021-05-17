#==============================================================================#
#                                                                              #
#                     Nimble models with several predictors                    #
#                                                                              #
#==============================================================================#


# 1. Individual traits + ice-free days (common effect) =========================

model_1 <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days
      b[3]*size_s[i] + b[4]*size_s2[i] +                       # Size
      b[5]*age_young[i] + b[6]*age_old[i] +                    # Age 
      b[7]*day_s[i] +                                          # Day of capture
      b[8]*age_young[i]*day_s[i] + b[9]*age_old[i]*day_s[i] +  # Interaction day of capture * age
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[10] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[3]*size_s[i] + b[4]*size_s2[i] + 
      b[5]*age_young[i] + b[6]*age_old[i] +
      b[7]*day_s[i] +
      b[8]*age_young[i]*day_s[i] + b[9]*age_old[i]*day_s[i] + 
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:10){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})



# 2. Ice-free days previous year: distinct slopes  (1.1.2_D_distinct) ========================

model_1.1.2_D_common_5A <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b[1] + b[2]*ice_free_days_previous_s[i] + b[3]*age_s[i] + b[4]*age2_s[i] + 
      b[5]*size_s[i] + b[6]*size2_s[i] + b[7]*day_s[i] + sigma1 * eta1[year[i]]
    
    log(q[i, 3]) <- b[8] +  b[2]*ice_free_days_previous_s[i] + b[3]*age_s[i] + b[4]*age2_s[i] + 
      b[5]*size_s[i] + b[6]*size2_s[i] + b[7]*day_s[i] + sigma1 * eta1[year[i]]
    
    log(q[i, 4]) <- b[9] + b[2]*ice_free_days_previous_s[i] + b[5]*size_s[i] + 
      sigma1 * eta1[year[i]]
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] <- sigma1 * eta1[i]
    eta1[i] ~ dnorm(0.00000E+00, 1)
  }
  
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for(i in 1:9){
    b[i] ~ dnorm(0.00000E+00, 0.1)
  }  
}
