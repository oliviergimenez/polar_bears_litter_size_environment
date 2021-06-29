#==============================================================================#
#                                                                              #
#                          Effect structure analysis                           #
#                                                                              #
#==============================================================================#

library(tidyverse)
Sys.setenv(LANG = "en")


# ~ 1. Plot the wAIC of each single-covariate model ---------------------------- 
wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
                       col_types = cols(wAIC_1 = col_double(),
                                        wAIC_2 = col_double(), 
                                        wAIC_3 = col_double(), 
                                        wAIC_4 = col_double(), 
                                        wAIC_5 = col_double())) %>%
  filter(!is.na(wAIC_1))

wAIC_table <- wAIC_table %>%
  mutate(mean_wAIC = apply(X = wAIC_table[, c("wAIC_1", "wAIC_2", "wAIC_3",
                                              "wAIC_4", "wAIC_5")],
                           MARGIN = 1, FUN = mean, na.rm = TRUE),
         variable = factor(variable, levels = unique(wAIC_table$variable)))

levels(wAIC_table$variable) <- c("null_model", "ice free days t-1", "ice-free days t-2", 
                                       "sea ice retreat date t-1", "sea ice retreat date t-2",
                                       "sea ice advance date t-1", "sea ice advance date t-2",
                                       "winter AO t", "winter AO t-1", "winter AO t-2", 
                                       "spring AO t", "spring AO t-1", "spring AO t-2",
                                       "winter NAO t", "winter NAO t-1", "winter NAO t-2", 
                                       "spring NAO t", "spring NAO t-1", "spring NAO t-2", 
                                 "day capture")


ggplot(data = wAIC_table[-1,], aes(x = effect, y = mean_wAIC)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ variable) +
  geom_text(aes(label = mean_wAIC), nudge_y = 0.3, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("07_results/01_interim_results/model_outputs/graph/wAIC_effect_structure.png",
       width = 10, height = 8.5)



# ~ 2. Plot paremeter estimates with 80% CI ------------------------------------

df.model.outputs <- data.frame(value = NA,
                               parameter = NA,
                               model_code = NA,
                               effect = NA,
                               variable = NA) %>%
  filter(!is.na(value))

for (k in 2:nrow(wAIC_table)) {
  
  load(file = paste0("07_results/01_interim_results/model_outputs/model_", 
                     wAIC_table$model_code[k], "_effect_", wAIC_table$effect[k], ".RData"))
  nimble.output <- get(paste0("fit_", wAIC_table$model_code[k], 
                              "_effect_", wAIC_table$effect[k]))
  
  if (wAIC_table$effect[k] %in% c("1c_VS_0c", "common")) {
    value <-  c(nimble.output$samples$chain1[, "a1"], 
              nimble.output$samples$chain2[, "a1"])
    df.model_outputs_k <- data.frame(value = value,
                                     parameter = rep("a1", times = length(value)),
                                     model_code = rep(wAIC_table$model_code[k], times = length(value)),
                                     effect = rep(wAIC_table$effect[k], times = length(value)),
                                     variable = rep(wAIC_table$variable[k], times = length(value)))
    # value <- c(nimble.output$samples$chain1[1:dim(nimble.output$samples$chain1)[1]], 
    #             nimble.output$samples$chain2[1:dim(nimble.output$samples$chain2)[1]])
    # parameter <- rep("a1", times = length(value))
    
    df.model.outputs <- rbind(df.model.outputs,
                              df.model_outputs_k)
  }
  if (wAIC_table$effect[k] == "2_3c_VS_0c") {
    value <- c(nimble.output$samples$chain1[, "b1"], 
              nimble.output$samples$chain2[, "b1"])
    df.model_outputs_k <- data.frame(value = value,
                                     parameter = rep("b1", times = length(value)),
                                     model_code = rep(wAIC_table$model_code[k], times = length(value)),
                                     effect = rep(wAIC_table$effect[k], times = length(value)),
                                     variable = rep(wAIC_table$variable[k], times = length(value)))
    # value <- c(nimble.output$samples$chain1[1:dim(nimble.output$samples$chain1)[1]], 
    #             nimble.output$samples$chain2[1:dim(nimble.output$samples$chain2)[1]])
    # parameter <- rep("a1", times = length(value))
    
    df.model.outputs <- rbind(df.model.outputs,
                              df.model_outputs_k)
  }
  if (wAIC_table$effect[k] == "distinct") {
    value <- c(nimble.output$samples$chain1[, "a1"], 
              nimble.output$samples$chain2[, "a1"],
              nimble.output$samples$chain1[, "b1"], 
              nimble.output$samples$chain2[, "b1"],
              nimble.output$samples$chain1[, "a1"] - nimble.output$samples$chain1[, "b1"],
              nimble.output$samples$chain2[, "a1"] - nimble.output$samples$chain2[, "b1"])
    df.model_outputs_k <- data.frame(value = value,
                                     parameter = c(rep("a1", times = length(value)/3),
                                                   rep("b1", times = length(value)/3),
                                                   rep("a1-b1", times = length(value)/3)),
                                     model_code = rep(wAIC_table$model_code[k], times = length(value)),
                                     effect = rep(wAIC_table$effect[k], times = length(value)),
                                     variable = rep(wAIC_table$variable[k], times = length(value)))
    # value <- c(nimble.output$samples$chain1[1:dim(nimble.output$samples$chain1)[1]], 
    #             nimble.output$samples$chain2[1:dim(nimble.output$samples$chain2)[1]])
    # parameter <- rep("a1", times = length(value))
    
    df.model.outputs <- rbind(df.model.outputs,
                              df.model_outputs_k)
  }
  rm(list = c(paste0("fit_", wAIC_table$model_code[k], "_effect_", wAIC_table$effect[k])))
}
df.model.outputs <- df.model.outputs %>% 
  mutate(param_effect = paste0(effect, "_", parameter),
         variable = factor(variable, levels = unique(df.model.outputs$variable)))

levels(df.model.outputs$variable) <- c("ice free days t-1", "ice-free days t-2", 
                                       "sea ice retreat date t-1", "sea ice retreat date t-2",
                                       "sea ice advance date t-1", "sea ice advance date t-2",
                                       "winter AO t", "winter AO t-1", "winter AO t-2", 
                                       "spring AO t", "spring AO t-1", "spring AO t-2",
                                       "winter NAO t", "winter NAO t-1", "winter NAO t-2", 
                                       "spring NAO t", "spring NAO t-1", "spring NAO t-2",
                                       "day capture")
unique(df.model.outputs$variable)


get_mean_and_CI <- function(x, lower, upper) {
  output <- data.frame(y = mean(x), 
                       ymin = quantile(x, probs = lower),
                       ymax = quantile(x, probs = upper))
  return(output)
}


ggplot(data = df.model.outputs, aes(x = param_effect , y = value, color = parameter)) + 
  stat_summary(fun.data = get_mean_and_CI, 
               fun.args = list(lower = 0.20, upper = 0.80),
               geom = "pointrange") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # stat_summary(fun.data = mean_sdl, 
  #              geom = "pointrange") +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) + 
                     #values = c("#FFAE03", "#FE4E00", "#E9190F")) +
  scale_x_discrete(limits = c("1c_VS_0c_a1", "2_3c_VS_0c_b1", "common_a1", 
                                "distinct_a1", "distinct_b1", "distinct_a1-b1"),
                   # labels = c("1c VS 0c: a1", "2-3c VS 0c: b1", "common: a1",
                   #            "distinct: a1", "distinct: b1", "distinct: a1-b1"),
                   labels = c("a1 (1c VS 0c)", "b1 (2-3c VS 0c)", "a1 (common)",
                              "a1 (distinct)", "b1 (distinct)", "a1-b1 (distinct)")) +
  facet_wrap(~ variable) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "",
       color = "")

ggsave("07_results/01_interim_results/model_outputs/graph/estimates_effect_structure.png",
       width = 10, height = 8.5)




# ~ 3. Select the right effect structure for obvious cases ---------------------

wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
                       col_types = cols(wAIC_1 = col_double(),
                                        wAIC_2 = col_double(), 
                                        wAIC_3 = col_double(), 
                                        wAIC_4 = col_double(), 
                                        wAIC_5 = col_double())) %>%
  filter(!is.na(wAIC_1),
         variable != "none")

wAIC_table <- wAIC_table %>%
  mutate(mean_wAIC = apply(X = wAIC_table[, c("wAIC_1", "wAIC_2", "wAIC_3",
                                              "wAIC_4", "wAIC_5")],
                           MARGIN = 1, FUN = mean, na.rm = TRUE))

wAIC_table_2 <- c()
covariates <- unique(wAIC_table$variable)
for (k in 1:length(covariates)) {
  wAIC_table_k <- wAIC_table %>%
    filter(variable == covariates[k]) 
  
  min_wAIC_k <- min(wAIC_table_k$mean_wAIC)
  
  wAIC_table_k <- wAIC_table_k %>%
    mutate(delta_wAIC = mean_wAIC - min_wAIC_k) %>%
    filter(delta_wAIC < 2)
  
  wAIC_table_2 <- rbind(wAIC_table_2, 
                        wAIC_table_k)
  
}

wAIC_table_2 <- wAIC_table_2 %>%
  mutate(variable = factor(variable, levels = unique(wAIC_table$variable)))

levels(wAIC_table_2$variable) <- c("ice free days t-1", "ice-free days t-2", 
                                 "sea ice retreat date t-1", "sea ice retreat date t-2",
                                 "sea ice advance date t-1", "sea ice advance date t-2",
                                 "winter AO t", "winter AO t-1", "winter AO t-2", 
                                 "spring AO t", "spring AO t-1", "spring AO t-2",
                                 "winter NAO t", "winter NAO t-1", "winter NAO t-2", 
                                 "spring NAO t", "spring NAO t-1", "spring NAO t-2",
                                 "day capture")


ggplot(data = wAIC_table_2, aes(x = effect, y = mean_wAIC)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ variable) +
  geom_text(aes(label = mean_wAIC), nudge_y = 0.3, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("07_results/01_interim_results/model_outputs/graph/wAIC_effect_structure_2.png",
       width = 10, height = 8.5)



# ~ 4. Visually select the best model structure for ambiguous cases ------------- 

# For covariates for which the best effect structure was not 2 wAIC points away 
# from the second best, I will assess which effect structure is the best. To do so
# I will use the 80% CI plots of the parameters estimates. 
# If for instance the distinct effect model appears slightly better than the
# common effect model, but that the 80% CI of the two parameters overlap, then 
# I'll select the common effect model, since it has less parameters.


wAIC_table_3 <- c()

wAIC_table_3 <- rbind(wAIC_table_3,
                      wAIC_table_2[which(wAIC_table_2$model_code == "1.1.2_D" &    # ice free days t-1
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "1.2.2_D" &    # sea ice retreat t-1
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "1.2.3_D" &    # sea ice retreat t-2
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "1.3.2_D" &    # sea ice advance t-1
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "2.2.2.1"), ], # Spring AO t-1
                      
                      wAIC_table_2[which(wAIC_table_2$model_code == "2.2.3.1"&     # Spring AO t-2
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.1.1.1"&     # Winter NAO t
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.1.2.1"&     # Winter NAO t-1
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.1.3.1"&     # Winter NAO t-2
                                          wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.2.1.1"&     # Spring NAO t
                                     wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.2.2.1"&     # Spring NAO t-1
                                           wAIC_table_2$effect == "common"), ],
                      wAIC_table_2[which(wAIC_table_2$model_code == "3.2.3.1"&     # Spring NAO t-2
                                           wAIC_table_2$effect == "common"), ])
