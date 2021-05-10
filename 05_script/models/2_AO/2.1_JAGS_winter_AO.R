#==============================================================================#
#                                                                              #
#                         Jags models for Winter AO                           #
#                                                                              #
#==============================================================================#


# 1. Winter AO: common slope (2.1.1.1_common) =================================

model_2.1.1.1_common_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*w_AO_s[i] + sigma1 * eta1[year[i]]
    
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



# 2. Winter AO: distinct slopes (2.1.1.1_distinct) ==============================

model_2.1.1.1_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_s[i] + sigma1 * eta1[year[i]]
    
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



# 3. Previous Winter AO: common slope (2.1.2.1_common) =============================

model_2.1.2.1_common_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    
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



# 4. Previous Winter AO: distinct slopes  (2.1.2.1_distinct) ========================

model_2.1.2.1_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    
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



# 5. 2y prior Winter AO: common slope (2.1.3.1_common) ==========================

model_2.1.3.1_common_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + a1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + a1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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




# 6. 2y prior Winter AO: distinct slopes  (2.1.3.1_distinct) =====================

model_2.1.3.1_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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

# 7. Winter AO: distinct slopes FACTOR (2.1.1.2_distinct) =======================

model_2.1.1.2_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_pos[i] + sigma1 * eta1[year[i]]
    
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


# 8. Previous Winter AO: distinct slopes FACTOR (2.1.2.2_distinct) =======================

model_2.1.2.2_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_prior_pos[i] + sigma1 * eta1[year[i]]
    
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


# 9. 2y priro Winter AO: distinct slopes FACTOR (2.1.3.2_distinct) ============

model_2.1.3.2_distinct_slope <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + a1*w_AO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- b0 + b1*w_AO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 4]) <- c0 + c1*w_AO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    
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





# 10. Winter AO: common slope (2.1.1.1_common_bis) =================================

model_2.1.1.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_AO_s[i] + sigma1 * eta1[year[i]]
    
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



# 11. Winter AO w/o females w/ no cubs: distinct slopes (2.1.1.1_distinct_bis) ==============================

model_2.1.1.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_s[i] + sigma1 * eta1[year[i]]
    
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



# 12. Previous Winter AO w/o females w/ no cubs: common slope (2.1.2.1_common_bis) =============================

model_2.1.2.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    
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



# 13. Previous Winter AO w/o females w/ no cubs: distinct slopes  (2.1.2.1_distinct_bis) ========================

model_2.1.2.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_prior_s[i] + sigma1 * eta1[year[i]]
    
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



# 14. 2y prior Winter AO w/o females w/ no cubs: common slope (2.1.3.1_common_bis) ==========================

model_2.1.3.1_common_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + b1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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




# 15. 2y prior Winter AO w/o females w/ no cubs: distinct slopes  (2.1.3.1_distinct_bis) =====================

model_2.1.3.1_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_2y_prior_s[i] + sigma1 * eta1[year[i]]
    
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
  
}






# 16. Winter AO FACTOR w/o females w/ no cubs: distinct slopes (2.1.1.2_distinct_bis) ====

model_2.1.1.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_pos[i] + sigma1 * eta1[year[i]]
    
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


# 17. Prior Winter AO FACTOR w/o females w/ no cubs: distinct slopes (2.1.2.2_distinct_bis) ====

model_2.1.2.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_prior_pos[i] + sigma1 * eta1[year[i]]
    
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


# 18. 2y prior Winter AO FACTOR w/o females w/ no cubs: distinct slopes (2.1.3.2_distinct_bis) ====

model_2.1.3.2_distinct_slope_bis <- function()  {
  
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- b0 + b1*w_AO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    log(q[i, 3]) <- c0 + c1*w_AO_2y_prior_pos[i] + sigma1 * eta1[year[i]]
    
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



