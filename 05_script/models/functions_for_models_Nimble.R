#==============================================================================#
#                                                                              #
#                      Function to run and interpret models                    #
#                                                                              #
#==============================================================================#


# Null model -------------------------------------------------------------------
model_0.0.0.0 <- nimbleCode({
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

# Null model binomial ----------------------------------------------------------

model_0.0.0.0_binomial <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0
  }
  for (i in 1:nbyear) {
    eps1[i] ~ dnorm(0, sd = sigma1)
  }
  sigma1 ~ dunif(0, 10)
  b0 ~ dnorm(0.00000E+00, sd = 1.5)
  b1 ~ dnorm(0.00000E+00, sd = 1.5)
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


# Create the wAIC table 
# /!\ do not run it again or it will erase all the saved wAIC values !! 
# Changed the name of the file just in case
# wAIC_table <- data.frame(model_code = NA, 
#                          variable = NA,
#                          effect = NA,
#                          wAIC_1 = NA,
#                          wAIC_2 = NA,
#                          wAIC_3 = NA,
#                          wAIC_4 = NA,
#                          wAIC_5 = NA) %>%
#   filter(!is.na(variable))
# write_csv(wAIC_table, "07_results/01_interim_results/model_outputs/wAIC_table_new.csv")

save_wAIC_null_model <- function(model_code) {
  wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
                         col_types = cols(wAIC_1 = col_character(),
                                          wAIC_2 = col_character(), 
                                          wAIC_3 = col_character(), 
                                          wAIC_4 = col_character(), 
                                          wAIC_5 = col_character()))
  wAIC <- get(paste0("fit_", model_code))$WAIC
  row <- which(wAIC_table$model_code == model_code)
  if (length(row) == 1) {
    column.wAIC.to.fill <- which(is.na(wAIC_table[row, ]))
    if (length(column.wAIC.to.fill) > 0) {
      wAIC_table[row, min(column.wAIC.to.fill)] <- as.character(round(wAIC, 3))
    }
  }
  if (length(row) == 0) {
    row.to.add <- c(model_code, "none", "none (multinomial)", as.character(round(wAIC, 3)), NA, NA, NA, NA)
    wAIC_table <- rbind(wAIC_table,
                        row.to.add) %>%
      arrange(model_code) 
    colnames(wAIC_table) <- c("model_code", "variable", "effect", "wAIC_1", 
                              "wAIC_2", "wAIC_3", "wAIC_4", "wAIC_5")
    
  }
  # Remove potential duplicated values
  for (i in 1:nrow(wAIC_table)) {
    x <- as.numeric(as.vector(wAIC_table[i, c("wAIC_1",            # Keep only the columns with the wAIC
                                              "wAIC_2", "wAIC_3", "wAIC_4", "wAIC_5")]))
    
    y <- x[-which(is.na(x))]         # Remove the columns with NAs among the wAIC columns
    
    new.wAIC.columns <- unique(y)
    new.k.row <- c(new.wAIC.columns, rep(NA, times = 5 - length(new.wAIC.columns)))
    for (j in 1:5) {
      wAIC_table[i, j + 3] <- as.character(new.k.row[j])
    }
  }
  # print(paste0("wAIC: ", round(wAIC, 3)))
  print(round(wAIC, 3))
  return(wAIC_table)
}



save_wAIC <- function(model_code, var_short_name, effect) {
  # wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
  #                        col_types = cols(wAIC_1 = col_double(),
  #                                         wAIC_2 = col_double(), 
  #                                         wAIC_3 = col_double(), 
  #                                         wAIC_4 = col_double(), 
  #                                         wAIC_5 = col_double()))
  wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
                         col_types = cols(wAIC_1 = col_character(),
                                          wAIC_2 = col_character(), 
                                          wAIC_3 = col_character(), 
                                          wAIC_4 = col_character(), 
                                          wAIC_5 = col_character()))
  wAIC <- get(paste0("fit_", model_code, "_effect_", effect))$WAIC
  row <- which(wAIC_table$model_code == model_code & wAIC_table$effect == effect)
  if (length(row) == 1) {
    column.wAIC.to.fill <- which(is.na(wAIC_table[row, ]))
    if (length(column.wAIC.to.fill) > 0) {
      wAIC_table[row, min(column.wAIC.to.fill)] <- as.character(round(wAIC, 3))
    }
  }
  if (length(row) == 0) {
    row.to.add <- c(model_code, substr(var_short_name, 0, nchar(var_short_name)-2), 
                    effect, as.character(round(wAIC, 3)), NA, NA, NA, NA)
    wAIC_table <- rbind(wAIC_table,
                        row.to.add) %>%
      arrange(model_code) 
    colnames(wAIC_table) <- c("model_code", "variable", "effect", "wAIC_1", 
                              "wAIC_2", "wAIC_3", "wAIC_4", "wAIC_5")
    
  }
  # Remove potential duplicated values
  for (i in 1:nrow(wAIC_table)) {
    x <- as.numeric(as.vector(wAIC_table[i, c("wAIC_1",            # Keep only the columns with the wAIC
                                              "wAIC_2", "wAIC_3", "wAIC_4", "wAIC_5")]))
    
    y <- x[-which(is.na(x))]         # Remove the columns with NAs among the wAIC columns
    
    new.wAIC.columns <- unique(y)
    new.k.row <- c(new.wAIC.columns, rep(NA, times = 5 - length(new.wAIC.columns)))
    for (j in 1:5) {
      wAIC_table[i, j + 3] <- as.character(new.k.row[j])
    }
  }
  # print(paste0("wAIC: ", round(wAIC, 3)))
  print(round(wAIC, 3))
  return(wAIC_table)
}





# check_convergence <- function(params, effect, model_code) {
#   # Load file
#   load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
#                      model_code, "_effect_", effect, ".RData"))
#   nimble_output <- get(paste0("fit_", model_code, "_effect_", effect))
#   # rm(list(get(paste0("fit_", model_code, "_effect_", effect))))
#   print(paste0("check convergence of model_", model_code, "_effect_", effect))
#   
#   # Process Nimble output into dataframe
#   chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
#     select(params[-length(params)]) %>%
#     mutate(chain = "chain 1",
#            iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
#   chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
#     select(params[-length(params)]) %>%
#     mutate(chain = "chain 2",
#            iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
#   chains <- rbind(chain1, chain2)
#   
#   
#   # Plot
#   if (effect %in% c("1c_VS_0c", "common")) {
#     nrows <- 4
#     
#     trace_a0 <- ggplot(data = chains, aes(x = iteration, y = a0, color = chain)) +
#       geom_line() +
#       labs(y = "a0")
#     density_a0 <- ggplot(data = chains, 
#                          aes(x = a0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "a0") +
#       theme(legend.position = "none")
#     
#     trace_a1 <- ggplot(data = chains, aes(x = iteration, y = a1, color = chain)) +
#       geom_line() +
#       labs(y = "a1")
#     density_a1 <- ggplot(data = chains, 
#                          aes(x = a1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "a1") +
#       theme(legend.position = "none")
#     
#     trace_b0 <- ggplot(data = chains, aes(x = iteration, y = b0, color = chain)) +
#       geom_line() +
#       labs(y = "b0")
#     density_b0 <- ggplot(data = chains, 
#                          aes(x = b0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "b0") +
#       theme(legend.position = "none")
#     
#     trace_sigma1 <- ggplot(data = chains, aes(x = iteration, y = sigma1, color = chain)) +
#       geom_line() +
#       labs(y = "sigma1")
#     density_sigma1 <- ggplot(data = chains, 
#                              aes(x = sigma1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "sigma1") +
#       theme(legend.position = "none")
#     
#     diagnostic_plot <- plot_grid(trace_a0, density_a0,
#               trace_a1, density_a1,
#               trace_b0, density_b0,
#               trace_sigma1, density_sigma1,
#               ncol = 2, nrow = nrows)
#   } 
#   if (effect == "2_3c_VS_0c") {
#     nrows <- 4
#     
#     trace_a0 <- ggplot(data = chains, aes(x = iteration, y = a0, color = chain)) +
#       geom_line() +
#       labs(y = "a0")
#     density_a0 <- ggplot(data = chains, 
#                          aes(x = a0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "a0") +
#       theme(legend.position = "none")
#     
#     trace_b0 <- ggplot(data = chains, aes(x = iteration, y = b0, color = chain)) +
#       geom_line() +
#       labs(y = "b0")
#     density_b0 <- ggplot(data = chains, 
#                          aes(x = b0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "b0") +
#       theme(legend.position = "none")
#     
#     trace_b1 <- ggplot(data = chains, aes(x = iteration, y = b1, color = chain)) +
#       geom_line() +
#       labs(y = "b1")
#     density_b1 <- ggplot(data = chains, 
#                          aes(x = b1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "b1") +
#       theme(legend.position = "none")
#     
#     trace_sigma1 <- ggplot(data = chains, aes(x = iteration, y = sigma1, color = chain)) +
#       geom_line() +
#       labs(y = "sigma1")
#     density_sigma1 <- ggplot(data = chains, 
#                              aes(x = sigma1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "sigma1") +
#       theme(legend.position = "none")
#     
#     diagnostic_plot <- plot_grid(trace_a0, density_a0,
#               trace_b0, density_b0,
#               trace_b1, density_b1,
#               trace_sigma1, density_sigma1,
#               ncol = 2, nrow = nrows)
#   }
#   if (effect == "distinct") {
#     nrows <- 5
#     
#     trace_a0 <- ggplot(data = chains, aes(x = iteration, y = a0, color = chain)) +
#       geom_line() +
#       labs(y = "a0")
#     density_a0 <- ggplot(data = chains, 
#                          aes(x = a0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "a0") +
#       theme(legend.position = "none")
#     
#     trace_a1 <- ggplot(data = chains, aes(x = iteration, y = a1, color = chain)) +
#       geom_line() +
#       labs(y = "a1")
#     density_a1 <- ggplot(data = chains, 
#                          aes(x = a1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "a1") +
#       theme(legend.position = "none")
#     
#     trace_b0 <- ggplot(data = chains, aes(x = iteration, y = b0, color = chain)) +
#       geom_line() +
#       labs(y = "b0")
#     density_b0 <- ggplot(data = chains, 
#                          aes(x = b0, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "b0") +
#       theme(legend.position = "none")
#     
#     trace_b1 <- ggplot(data = chains, aes(x = iteration, y = b1, color = chain)) +
#       geom_line() +
#       labs(y = "b1")
#     density_b1 <- ggplot(data = chains, 
#                          aes(x = b1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "b1") +
#       theme(legend.position = "none")
#     
#     trace_sigma1 <- ggplot(data = chains, aes(x = iteration, y = sigma1, color = chain)) +
#       geom_line() +
#       labs(y = "sigma1")
#     density_sigma1 <- ggplot(data = chains, 
#                              aes(x = sigma1, color = chain, fill = chain)) +
#       geom_density(alpha = 0.25) +
#       labs(x = "sigma1") +
#       theme(legend.position = "none")
#     
#     diagnostic_plot <- plot_grid(trace_a0, density_a0,
#               trace_a1, density_a1,
#               trace_b0, density_b0,
#               trace_b1, density_b1,
#               trace_sigma1, density_sigma1,
#               ncol = 2, nrow = nrows)
#     
#   }
#   save_plot(filename = paste0("07_results/01_interim_results/model_outputs/graph/diagnostic_plots/fit_",
#                               model_code, "_effect_", effect, ".png"), 
#             plot = diagnostic_plot,
#             ncol = 2,
#             nrow = nrows)
#             # width = 10,
#             # height = nrows * 3)
#   return(diagnostic_plot)
# }




check_convergence <- function(params, effect, model_code) {
  # Load file
  load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                     model_code, "_effect_", effect, ".RData"))
  nimble_output <- get(paste0("fit_", model_code, "_effect_", effect))
  # rm(list(get(paste0("fit_", model_code, "_effect_", effect))))
  print(paste0("check convergence of model_", model_code, "_effect_", effect))
  
  # Process Nimble output into dataframe
  chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
    dplyr::select(params[-length(params)]) %>%
    mutate(chain = "1",
           iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
  chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
    dplyr::select(params[-length(params)]) %>%
    mutate(chain = "2",
           iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
  chains <- rbind(chain1, chain2) 
  
  chains_l <- pivot_longer(chains, cols = params[-length(params)], names_to = "parameter") 
  
  param.mean <- chains_l %>%
    group_by(parameter, chain) %>%
    summarize(m = mean(value))
  
  param.running.mean <- chains_l %>%
    arrange(parameter, iteration) %>%
    group_by(parameter, chain) %>%
    mutate(rm = cumsum(value)/iteration)
  
  trace.plots <- ggplot(data = chains_l, 
                        aes(x = iteration, y = value, color = chain)) +
    geom_line() +
    labs(y = "trace") +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free",
                ncol = 1)
  
  density.plots <- ggplot(data = chains_l, 
                          aes(x = value, color = chain, fill = chain)) +
    geom_density(alpha = 0.25) +
    labs(x = "density") +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free_y",
                ncol = 1)
  
  running.mean.plot <- ggplot(param.running.mean, 
                              aes(x = iteration, y = rm, color = chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept = m), param.mean,
               colour = "black", alpha = 0.5) + 
    ylab("running Mean") +
    facet_grid(parameter ~ chain, scales = "free")
    
  # Plot all the plots together
  diagnostic_plot <- plot_grid(trace.plots,
                               density.plots, 
                               running.mean.plot,
                               ncol = 3, nrow = 1)
 
  nrows = length(params)
  save_plot(filename = paste0("07_results/01_interim_results/model_outputs/graph/diagnostic_plots/fit_",
                              model_code, "_effect_", effect, ".png"), 
            plot = diagnostic_plot,
            ncol = 3,
            nrow = 1,
            base_width = 10,
            base_height = nrows * 3)
  return(diagnostic_plot)
}



check_convergence_several_predictors <- function(params.plot, nimble_output) {
  # Process Nimble output into dataframe
  chain1 <- data.frame(nimble_output[["samples"]][["chain1"]]) %>%
    dplyr::select(params.plot) %>%
    mutate(chain = "1",
           iteration = seq(1, dim(nimble_output[["samples"]][["chain1"]])[1], by = 1))
  chain2 <- data.frame(nimble_output[["samples"]][["chain2"]]) %>%
    dplyr::select(params.plot) %>%
    mutate(chain = "2",
           iteration = seq(1, dim(nimble_output[["samples"]][["chain2"]])[1], by = 1))
  chains <- rbind(chain1, chain2) 
  
  chains_l <- pivot_longer(chains, cols = params.plot, names_to = "parameter") 
  
  param.mean <- chains_l %>%
    group_by(parameter, chain) %>%
    summarise(m = mean(value))
  
  param.running.mean <- chains_l %>%
    arrange(parameter, iteration) %>%
    group_by(parameter, chain) %>%
    mutate(rm = cumsum(value)/iteration)
  
  trace.plots <- ggplot(data = chains_l, 
                        aes(x = iteration, y = value, color = chain)) +
    geom_line() +
    labs(y = "trace") +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free",
                ncol = 1)
  
  density.plots <- ggplot(data = chains_l, 
                          aes(x = value, color = chain, fill = chain)) +
    geom_density(alpha = 0.25) +
    labs(x = "density") +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free_y",
                ncol = 1)
  
  running.mean.plot <- ggplot(param.running.mean, 
                              aes(x = iteration, y = rm, color = chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept = m), param.mean,
               colour = "black", alpha = 0.5) + 
    ylab("running Mean") +
    facet_grid(parameter ~ chain, scales = "free")
  
  # Plot all the plots together
  diagnostic_plot <- plot_grid(trace.plots,
                               density.plots, 
                               running.mean.plot,
                               ncol = 3, nrow = 1)
  return(diagnostic_plot)
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




get_probabilities_factor <- function(model_code, effect, mode, var_scaled, var) {
  res <- rbind(get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain1,
               get(paste0("fit_", model_code, "_effect_", effect, mode))$samples$chain2)
  
  # Create grid of x values
  lengthgrid <- 2
  grid <- seq(from = 0,
              to = 1, 
              length = lengthgrid)
  
  q0cub <- q1cub <- q2_3cub <- matrix(data = NA, 
                                      nrow = dim(res)[1], 
                                      ncol = lengthgrid)
  
  if (effect == "1c_VS_0c") {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, 3] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid[i])	
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
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid[i])	
      }
    }
  }
  if (effect %in% c("common")) {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, c(3, 2)] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid[i])	
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid[i])
      }
    }
  }
  if (effect %in% c("distinct")) {
    b1cub <- res[, c(1, 2)]
    b2_3cub <- res[, c(3, 4)] 
    
    for (i in 1:lengthgrid) {
      for (j in 1:dim(res)[1]) {
        q0cub[j, i] <- exp(0)
        q1cub[j, i] <- exp(b1cub[j, 1] + b1cub[j, 2] * grid[i])	
        q2_3cub[j, i] <- exp(b2_3cub[j, 1] + b2_3cub[j, 2] * grid[i])
      }
    }
  }
  # backtransform
  p0cub <- p1cub <- p2_3cub <- matrix(NA, dim(res)[1], lengthgrid)
  for (i in 1:lengthgrid){
    for (j in 1:dim(res)[1]){
      norm <- (q0cub[j, i] + q1cub[j, i] + q2_3cub[j, i])
      p0cub[j, i] <- q0cub[j, i]/norm
      p1cub[j, i] <- q1cub[j, i]/norm
      p2_3cub[j, i] <- q2_3cub[j, i]/norm
    }
  }
  df.for.plot <- data.frame(var = c(rep(0, times = length(p1cub)/2),
                                    rep(1, times = length(p1cub)/2),
                                    rep(0, times = length(p1cub)/2),
                                    rep(1, times = length(p1cub)/2),
                                    rep(0, times = length(p1cub)/2),
                                    rep(1, times = length(p1cub)/2)),
                            probability = c(p0cub[, 1], p0cub[, 2],
                                            p1cub[, 1], p1cub[, 2],
                                            p2_3cub[, 1], p2_3cub[, 2]),
                            nbr_cub = c(rep("no cubs", times = length(p1cub)),
                                        rep("1 cub", times = length(p1cub)),
                                        rep("2-3 cubs", times = length(p1cub)))) %>%
    mutate(nbr_cub = factor(nbr_cub,      # Reordering group factor levels
                            levels = c("no cubs", "1 cub", "2-3 cubs")))
  color_labels <- c("no cubs", "1 cub", "2-3 cubs")
  return(list(df.for.plot, color_labels))
}

