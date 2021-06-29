#==============================================================================#
#                                                                              #
#                  Nimble models with environmental covariates                 #
#                                                                              #
#==============================================================================#


# A. Multinomial ===============================================================

# ~ 1. Ice-free days + spring AO (both common effect) =============================

model_1.1.2.D_2.2.2.1_common <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*prior_spring_AO_s[i] +                              # Spring AO
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[4] + 
      b[2]*ice_free_days_previous_s[i] +                                     
      b[3]*prior_spring_AO_s[i] +                               
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:4){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 2. Ice-free days + spring AO (both common effect) + day capture (distinct) ====

model_1.1.2.D_2.2.2.1_common_4_distinct <- nimbleCode({
  
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
      b[6]*day_s[i] +                                          
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:6){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 3. Ice-free days * day capture (distinct) + spring AO (both common effect)  =====

# Here the coefficient for the interaction is common

model_1.1.2.D_2.2.2.1_common_4_distinct_2 <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*prior_spring_AO_s[i] +                              # Spring AO
      b[4]*day_s[i] +                                          # Day of capture
      b[5]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[6] + 
      b[2]*ice_free_days_previous_s[i] +                                     
      b[3]*prior_spring_AO_s[i] +                               
      b[7]*day_s[i] + 
      b[5]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:7){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 4. Ice-free days * day capture (distinct) + spring AO (both common effect)  =====

# Here the coefficient for the interaction is distinct

model_1.1.2.D_2.2.2.1_common_4_distinct_3 <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*prior_spring_AO_s[i] +                              # Spring AO
      b[4]*day_s[i] +                                          # Day of capture
      b[5]*ice_free_days_previous_s[i]*day_s[i] +
    eps1[year[i]]                                              # Random effect of year
    
    log(q[i, 3]) <- b[6] + 
      b[2]*ice_free_days_previous_s[i] +                                     
      b[3]*prior_spring_AO_s[i] +                               
      b[7]*day_s[i] + 
      b[8]*ice_free_days_previous_s[i]*day_s[i] +                # Distinct interaction
    eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:8){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 5. Ice-free days t-1 + spring AO t-1 + winter AO + (all common effect) + day capture (distinct) ====

model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*winter_AO_s[i] +                                    # Winter AO t-1
      b[4]*prior_spring_AO_s[i] +                              # Spring AO t-1 
      b[5]*day_s[i] +                                          # Day of capture
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[6] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[3]*winter_AO_s[i] +
      b[4]*prior_spring_AO_s[i] +    
      b[7]*day_s[i] +                                          
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:7){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})



# ~ 6. Ice-free days t-1 * day capture (distinct) + spring AO t-1 + winter AO + (all common effect)  ====

# In this model, there is an interaction 
model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_2 <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*winter_AO_s[i] +                                    # Winter AO t-1
      b[4]*prior_spring_AO_s[i] +                              # Spring AO t-1 
      b[5]*day_s[i] +                                          # Day of capture
      b[6]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[7] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[3]*winter_AO_s[i] +
      b[4]*prior_spring_AO_s[i] +    
      b[8]*day_s[i] +                                          
      b[6]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:8){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 7. Ice-free days t-1 * day capture (distinct) + spring AO t-1 + winter AO + (all common effect)  ====

# In this model, there is an interaction. The coefficients for the interaction are distinct

model_1.1.2.D_2.1.1.1_2.2.2.1_common_4_distinct_3 <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*winter_AO_s[i] +                                    # Winter AO t-1
      b[4]*prior_spring_AO_s[i] +                              # Spring AO t-1 
      b[5]*day_s[i] +                                          # Day of capture
      b[6]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[7] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[3]*winter_AO_s[i] +
      b[4]*prior_spring_AO_s[i] +    
      b[8]*day_s[i] +                                          
      b[9]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:9){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 8. Ice-free days t-1 * day capture (distinct) ==============================

# In this model, there is an interaction. The coefficient for the interaction is common

model_1.1.2.D_common_4_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*day_s[i] +                                          # Day of capture
      b[4]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[5] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[6]*day_s[i] +                                          
      b[4]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:6){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})


# ~ 9. Ice-free days t-1 * day capture (distinct) + winter AO ==================

# In this model, there is an interaction. Thr coefficient is common
model_1.1.2.D_2.1.1.1_common_4_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*ice_free_days_previous_s[i] +                       # Ice free days              
      b[3]*winter_AO_s[i] +                                    # Winter AO t-1
      b[4]*day_s[i] +                                          # Day of capture
      b[5]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[6] + 
      b[2]*ice_free_days_previous_s[i] + 
      b[3]*winter_AO_s[i] +
      b[7]*day_s[i] +                                          
      b[5]*ice_free_days_previous_s[i]*day_s[i] +
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:7){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})



# ~ 10. Ice-free days t-1 * day capture (distinct) + spring AO t-1 + winter AO + (all common effect)  ====

# In this model, there is an interaction 
model_2.1.1.1_2.2.2.1_common_4_distinct <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    
    log(q[i, 2]) <- b[1] + 
      b[2]*winter_AO_s[i] +                                    # Winter AO t-1
      b[3]*prior_spring_AO_s[i] +                              # Spring AO t-1 
      b[4]*day_s[i] +                                          # Day of capture
      eps1[year[i]]                                            # Random effect of year
    
    log(q[i, 3]) <- b[5] + 
      b[2]*winter_AO_s[i] +
      b[3]*prior_spring_AO_s[i] +    
      b[6]*day_s[i] +                                          
      eps1[year[i]]
    
    
    for (j in 1:J) {
      p[i, j] <- q[i, j]/sum(q[i, 1:J])
    }
  }
  
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  
  sigma1 ~ dunif(0, 10)
  
  for(i in 1:6){
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  } 
})





#==============================================================================#

# B. Binomial ==================================================================

# ~ 1. Ice-free days + winter AO + spring AO + day capture =====================

model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * winter_AO_s[i] +
      b[4] * prior_spring_AO_s[i] +
      b[5] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:5) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})


# ~ 2. Ice-free days * day capture + winter AO + spring AO =====================

# Add an interaction between ice-free days and day of capture
model_1.1.2_D_2.1.1.1_2.2.2.1_4_binomial_2 <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * winter_AO_s[i] +
      b[4] * prior_spring_AO_s[i] +
      b[5] * day_s[i] +
      b[6] * ice_free_days_previous_s[i] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:6) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})


# ~ 3. Ice-free days + winter AO + spring AO + day capture =====================

# Remove winter AO

model_1.1.2_D_2.2.2.1_4_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * prior_spring_AO_s[i] +
      b[4] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:4) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})



# ~ 4. Ice-free days + winter AO + spring AO + day capture =====================

# Remove spring AO

model_1.1.2_D_2.1.1.1_4_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * winter_AO_s[i] +
      b[4] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:4) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})


# ~ 5. Ice-free days + winter AO + spring AO ===================================

model_1.1.2_D_2.1.1.1_2.2.2.1_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * winter_AO_s[i] +
      b[4] * prior_spring_AO_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:4) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})



# ~ 6. Winter AO + spring AO + day capture =====================================

# Remove ice-free days

model_2.1.1.1_2.2.2.1_4_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * winter_AO_s[i] +
      b[3] * prior_spring_AO_s[i] +
      b[4] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:4) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})


# ~ 7. Ice-free days + day capture =============================================

# remove winter AO and spring AO t-1

model_1.1.2_D_4_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:3) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})




# ~ 8. Ice-free days + day capture =============================================

# remove winter AO and spring AO t-1. Add the interaction between ice-free days and 
# date of capture

model_1.1.2_D_4_binomial_2 <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- b[1] + 
      b[2] * ice_free_days_previous_s[i] +
      b[3] * day_s[i] +
      b[4] * ice_free_days_previous_s[i] * day_s[i] +
      eps1[year[i]]
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0.00000E+00, 10)
  
  for (i in 1:4) {
    b[i] ~ dnorm(0.00000E+00, sd = 1.5)
  }
})







