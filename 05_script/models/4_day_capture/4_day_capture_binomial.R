#==============================================================================#
#                                                                              #
#                   Models with only date of capture (binomial)                #
#                                                                              #
#==============================================================================#


library(tidyverse)
library(viridis)
library(ggmcmc)
library(cowplot)
library(nimble)
Sys.setenv(LANG = "en")

# READY VARIABLES ===========================================================
# CR data
data_model <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv") %>%
  mutate(success = ifelse(cub_number == 1, 0, 1))

y <- data_model %>%
  pull(success)

# Renumérotation des années
{year <- data_model$year
  year <- factor(year) # laisse tomber les modalites inutiles 
  year2 <- NULL
  for (i in 1:length(year)){
    year2 <- c(year2, which(year[i] == levels(year)))
  }
  year <- factor(year2)
  nbyear <- length(levels(year))
  year
}


dat <- list(y = as.numeric(y))

# Load the Nimble models + the ancillary functions
source("05_script/models/4_day_capture/4_Nimble_day_capture.R")
source("05_script/models/functions_for_models_Nimble.R")



# RUN MODELS ===================================================================

# A. Null model ================================================================




# B. Model with only ice-free days (Bernoulli) =================================

# ~~~ a. Run the model ---------------------------------------------------------

{model_code <- "4"
effect <- "binomial"

var <- data_model$day_number
var_scaled <- scale(var)
var_short_name <- "day_s"
var_full_name <- "Day of capture"

# Are females without cubs taken into account ?

my.constants <- list(N = length(y),            # nb of females captured
                     year = as.numeric(year),
                     nbyear = nbyear,
                     as.numeric(var_scaled)) 
names(my.constants)[4] <- var_short_name
}

# Define the parameters to estimate
params <- c("b0", "b1", "sigma1", "eps1")

temp.data <- data.frame(y = y, 
                        var_scaled = var_scaled)
mylogit <- glm(as.factor(y) ~ var_scaled, 
               data = temp.data, 
               family = "binomial")

inits <- function() list(b0 = mylogit$coefficients[1] + round(runif(n = 1, -1, 1))/10, 
                         b1 = mylogit$coefficients[2] + round(runif(n = 1, -1, 1))/10,
                         sigma1 = runif(1))


# Run the model
start <- Sys.time()
assign(x = paste0("fit_", model_code, "_effect_", effect),
       value = nimbleMCMC(code = get(paste0("model_", model_code, "_effect_", effect)),     # model code  
                          data = dat,                                   
                          constants = my.constants,        
                          inits = inits,          
                          monitors = params,   # parameters to monitor
                          thin = 10,
                          niter = 25000,                  # nb iterations
                          nburnin = 5000,              # length of the burn-in
                          nchains = 2,
                          summary = TRUE,
                          WAIC = TRUE))
end <- Sys.time()
end - start

get(paste0("fit_", model_code, "_effect_", effect))$WAIC
# 313.0506

save(list = paste0("fit_", model_code, "_effect_", effect), 
     file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))

# ~~~ b. Check convergence -----------------------------------------------------

check_convergence(params = params,
                  effect = effect,
                  model_code = model_code)


# ~~~ c. Plot ------------------------------------------------------------------

load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                   model_code, "_effect_", effect, ".RData"))


# Get df for ggplot
assign(x = paste0("fit_", model_code, "_effect_", effect, "_for_plot"),
       value = get_probabilities_binomial(model_code, effect, var_scaled, var))


ggplot(data = get(paste0("fit_", model_code, "_effect_", effect, "_for_plot")), 
       aes(x = var, y = value, linetype = type, group = name)) +
  geom_line() +

  scale_linetype_manual(limits = c("mean", "credible_interval"),
                        values = c("solid", "dotted"),
                        labels = c("Mean", "CI")) +
  theme_bw() +
  labs(x = var_full_name,
       y = "Probability of having 2-3 cubs", 
       linetype = "") 

ggsave(filename = paste0("D:/polar_bears_litter_size_environment/07_results/01_interim_results/model_outputs/graph/model_", 
                         model_code,  "_effect_", effect, ".png"),
       width = 6, height = 3)


test <- fit_4_effect_binomial_for_plot %>%
  dplyr::select(-type) %>%
  pivot_wider(names_from = name, values_from = value)

ggplot(data = test) +
  geom_line(aes(x = var, y = mean_p_2_3_cub)) +
  geom_ribbon(aes(x = var, ymin = ci_p_2_3_cub_2.5,
                  ymax = ci_p_2_3_cub_97.5),
              alpha = 0.15) +
  theme_bw() +
  labs(x = var_full_name,
       y = "Probability of having 2-3 cubs", 
       linetype = "") 

