#==============================================================================#
#                                                                              #
#                      Function to run and interpret models                    #
#                                                                              #
#==============================================================================#


null_model <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dcat(p[i, 1:J])
    log(q[i, 1]) <- 0
    log(q[i, 2]) <- a0 + eps1[year[i]] # eps[year[i]]
    log(q[i, 3]) <- b0 + eps1[year[i]] # eps[year[i]]
    
    
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
})



get_coefs_and_params <- function(y, var_scaled, effect, mode) {
  temp.dat <- data.frame(y = y, var_scaled = var_scaled)
  temp.dat$yfac <- as.factor(temp.dat$y)   # Ajout de Y en facteur
  mnl.dat <- mlogit.data(temp.dat, varying = NULL, choice = "yfac", shape = "wide") 
  mlogit.mod <- mlogit(yfac ~ 1| var_scaled, 
                       data = mnl.dat, 
                       reflevel = ifelse(mode == "_bis", "1", "0"))
  
  all_coefs <- as.vector(summary(mlogit.mod)$coefficients)
  
  if (mode == "") {
    
    if (effect == "1c_VS_0c") {
      coefs <- c(all_coefs[1:3])
      params <- c("a0", "b0", "a1", "sigma1" ,"eps1")
    }
    if (effect == "2_3c_VS_0c") {
      coefs <- c(all_coefs[c(1, 2, 4)])
      params <- c("a0", "b0", "b1", "sigma1" ,"eps1")
    }
    if (effect == "common") { # Same as first possibility
      coefs <- c(all_coefs[c(1:3, 3)])
      params <- c("a0", "b0", "a1", "sigma1" ,"eps1")
    }
    if (effect == "distinct") {
      coefs <- c(all_coefs[c(1:4)])
      params <- c("a0", "b0", "a1", "b1", "sigma1" ,"eps1")
    }
    
    
  } else {
    print("not done")
  }
  return(list(coefs = coefs,
              params = params))
}










check_convergence <- function(jags_output, model_code, slope, mode) {
  processed_output <- ggs(as.mcmc(jags_output)) %>%
    filter(Parameter %in% c("a0", "a1", "a2", "b0", "b1", "c0", "c1", "deviance"))
  
  f1 <- ggs_traceplot(processed_output) + 
    theme_bw() +
    theme(legend.position = "none")
  f2 <- ggs_density(processed_output) + 
    theme_bw()
  f3 <- ggs_running(processed_output) + 
    theme_bw() +
    theme(legend.position = "none")
  x <- grid.arrange(f1, f2, f3, ncol = 3, nrow = 1)
  
  nbr_rows <- length(unique(processed_output$Parameter))
  
  ggsave(x,
         filename = paste0("07_results/01_interim_results/model_outputs/graphs/diagnostic_plots/model_",
                           model_code, "_", slope, toupper(mode), ".png"),
         width = 17.5, height = 2.5*nbr_rows)
  
}

get_probabilities <- function(model_code, effect, mode, var_scaled, var) { 
  res <- rbind(get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain1,
               get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain2)
  
  # Create grid of x values
  range <- range(var_scaled)
  
  lengthgrid <- 100
  grid_scaled <- seq(from = range[1] - 0.1*(range[2] - range[1]), 
                     to = range[2] + 0.1*(range[2] - range[1]), 
                     length = lengthgrid)
  
  grid <- grid_scaled * sd(var) + mean(var)
  
  
  q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                      nrow = dim(res)[1], 
                                      ncol = lengthgrid)
  
  if (effect == "1c_VS_0c") {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, 3] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
        q2_3cub[j, i] <- exp(b2_3cub[j])	
      }
    }
  }
  if (effect == "2_3c_VS_0c") {
    b1cub <- res[, 1]
    b2_3cub <- res[, c(2,3)] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j])	
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid_scaled[i])	
      }
    }
  }
  if (effect %in% c("common")) {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, c(3, 2)] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid_scaled[i])
      }
    }
  }
  if (effect %in% c("distinct")) {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, c(3, 4)] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid_scaled[i])	
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid_scaled[i])
      }
    }
  }
  # Backtransform
  p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
  for (i in 1:lengthgrid){
    for (j in 1:dim(res)[1]){
      norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
      p0cub[j, i] <- q0cub[j, i]/norm
      p1cub[j, i] <- q1cub[j, i]/norm
      p2_3cub[j, i] <- q2_3cub[j, i]/norm
    }
  }
  
  df.for.plot <- data.frame(var = grid,
                            mean_p_0_cub = apply(p0cub, 2, mean),
                            mean_p_1_cub = apply(p1cub, 2, mean),
                            mean_p_2_3_cub = apply(p2_3cub, 2, mean),
                            ci_p_0_cub_2.5 = apply(p0cub, 2, quantile, probs = 0.025),
                            ci_p_0_cub_97.5 = apply(p0cub, 2, quantile, probs = 0.975),
                            ci_p_1_cub_2.5 = apply(p1cub, 2, quantile, probs = 0.025),
                            ci_p_1_cub_97.5 = apply(p1cub, 2, quantile, probs = 0.975),
                            ci_p_2_3_cub_2.5 = apply(p2_3cub, 2, quantile, probs = 0.025),
                            ci_p_2_3_cub_97.5 = apply(p2_3cub, 2, quantile, probs = 0.975)) %>%
    pivot_longer(cols = c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub",
                          "ci_p_0_cub_2.5", "ci_p_0_cub_97.5", 
                          "ci_p_1_cub_2.5", "ci_p_1_cub_97.5", 
                          "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5")) %>%
    mutate(cub_number = ifelse(name %in% c("mean_p_0_cub", "ci_p_0_cub_2.5", "ci_p_0_cub_97.5"), 0, 
                               ifelse(name %in% c("mean_p_1_cub", "ci_p_1_cub_2.5", "ci_p_1_cub_97.5"), 1, 
                                      ifelse(name %in% c("mean_p_2_3_cub", "ci_p_2_3_cub_2.5", "ci_p_2_3_cub_97.5"), 2, 10))),
           type = ifelse(name %in% c("mean_p_0_cub", "mean_p_1_cub", "mean_p_2_3_cub"), "mean", "credible_interval"))
  
  color_labels <- c("no cubs", "1 cub", "2-3 cubs")
  return(list(df.for.plot, color_labels))
}




get_probabilities_factor <- function(model_code, mode, slope, var_scaled, var) {
  res <- get(paste0("fit_", model_code, "_", slope, mode))$BUGSoutput$sims.matrix
  
  # If females without cubs are included
  if(mode == "") {
    if(slope == "common") {
      b1cub <- res[, c(1, 2)]
      b2cub <- res[, c(3, 2)] 
      b3cub <- res[, c(4, 2)] 
    } else {
      b1cub <- res[, c(1, 2)]
      b2cub <- res[, c(3, 4)] 
      b3cub <- res[, c(5, 6)] 
    }
    
    lengthgrid <- 2
    grid <- seq(from = 0,
                to = 1, 
                length = lengthgrid)
    
    q0cub <- q1cub <- q2cub <- q3cub <- matrix(data = NA, 
                                               nrow = dim(b2cub)[1], 
                                               ncol = lengthgrid)
    for (i in 1:lengthgrid){
      for (j in 1:dim(b2cub)[1]){
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid[i])	
        q2cub[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid[i])	
        q3cub[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid[i])		
      }}
    # backtransform
    p0cub <- p1cub <- p2cub <- p3cub <- matrix(NA, dim(b2cub)[1], lengthgrid)
    for (i in 1:lengthgrid){
      for (j in 1:dim(b2cub)[1]){
        norm <- (q0cub[j, i] + q1cub[j, i] + q2cub[j, i] + q3cub[j,i])
        p0cub[j, i] <- q0cub[j, i]/norm
        p1cub[j, i] <- q1cub[j, i]/norm
        p2cub[j, i] <- q2cub[j, i]/norm
        p3cub[j, i] <- q3cub[j, i]/norm
      }
    }
    
    df.for.plot <- data.frame(var = c(rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2),
                                      rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2),
                                      rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2),
                                      rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2)),
                              probability = c(p0cub[, 1], p0cub[, 2],
                                              p1cub[, 1], p1cub[, 2],
                                              p2cub[, 1], p2cub[, 2],
                                              p3cub[, 1], p3cub[, 2]),
                              nbr_cub = c(rep("no cubs", times = length(p1cub)),
                                          rep("1 cub", times = length(p1cub)),
                                          rep("2 cubs", times = length(p1cub)),
                                          rep("3 cubs", times = length(p1cub)))) %>%
      mutate(nbr_cub = factor(nbr_cub,      # Reordering group factor levels
                              levels = c("no cubs", "1 cub", "2 cubs", "3 cubs")))
    color_labels <- c("no cubs", "1 cub", "2 cubs", "3 cubs")
    
    
    
    # If females without cubs are excluded
  } else {
    if(slope == "common") {
      b2cub <- res[, c(1, 2)]
      b3cub <- res[, c(3, 2)] 
    } else {
      b2cub <- res[, c(1, 2)] 
      b3cub <- res[, c(3, 4)]
    }
    
    range <- range(var_scaled)
    
    lengthgrid <- 2
    grid <- seq(from = 0,
                to = 1, 
                length = lengthgrid)
    
    q1cub <- q2cub <- q3cub <- matrix(data = NA, 
                                      nrow = dim(b2cub)[1], 
                                      ncol = lengthgrid)
    
    for (i in 1:lengthgrid){
      for (j in 1:dim(b2cub)[1]){
        q1cub[j, i] <- exp(0)
        q2cub[j, i] <- exp(b2cub[j, 1] + b2cub[j, 2] * grid[i])	
        q3cub[j, i] <- exp(b3cub[j, 1] + b3cub[j, 2] * grid[i])		
      }
    }
    
    # backtransform
    p1cub <- matrix(NA, dim(b2cub)[1], lengthgrid)
    p2cub <- p1cub
    p3cub <- p1cub
    for (i in 1:lengthgrid){
      for (j in 1:dim(b2cub)[1]){
        norm <- (q1cub[j, i] + q2cub[j, i] + q3cub[j,i])
        p1cub[j, i] <- q1cub[j, i]/norm
        p2cub[j, i] <- q2cub[j, i]/norm
        p3cub[j, i] <- q3cub[j, i]/norm
      }
    }
    df.for.plot <- data.frame(var = c(rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2),
                                      rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2),
                                      rep(0, times = length(p1cub)/2),
                                      rep(1, times = length(p1cub)/2)),
                              probability = c(p1cub[, 1], p1cub[, 2],
                                              p2cub[, 1], p2cub[, 2],
                                              p3cub[, 1], p3cub[, 2]),
                              nbr_cub = c(rep("1 cub", times = length(p1cub)),
                                          rep("2 cubs", times = length(p1cub)),
                                          rep("3 cubs", times = length(p1cub)))) %>%
      mutate(nbr_cub = factor(nbr_cub,      # Reordering group factor levels
                              levels = c("1 cub", "2 cubs", "3 cubs")))
    
    color_labels <- c("1 cub", "2 cubs", "3 cubs")
  }
  return(list(df.for.plot, color_labels))
}

