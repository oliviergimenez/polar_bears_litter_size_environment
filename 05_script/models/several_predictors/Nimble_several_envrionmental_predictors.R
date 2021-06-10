#==============================================================================#
#                                                                              #
#              Nimble models with several environmental predictors             #
#                                                                              #
#==============================================================================#


# 1.Ice-free days + spring AO (both common effect) =============================

model_1.1.2.D_2.2.2.1_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*prior_spring_AO_s[i] +                              # Spring AO
      b[4]*day_s[i] +                                          # Day of capture
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[5] + 
      b[2]*ice_free_days_previous_s[i] +                                     
      b[3]*prior_spring_AO_s[i] +                               
      b[4]*day_s[i] +                                          
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:5){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})



