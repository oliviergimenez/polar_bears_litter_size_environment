#==============================================================================#
#                                                                              #
#                               wAIC analysis                                  #
#                                                                              #
#==============================================================================#

library(tidyverse)
Sys.setenv(LANG = "en")



wAIC_table <- read_csv("07_results/01_interim_results/model_outputs/wAIC_table.csv", 
                       col_types = cols(wAIC_1 = col_double(),
                                        wAIC_2 = col_double(), 
                                        wAIC_3 = col_double(), 
                                        wAIC_4 = col_double(), 
                                        wAIC_5 = col_double()))

wAIC_table <- wAIC_table %>%
  mutate(mean_wAIC = apply(X = wAIC_table[, c("wAIC_1", "wAIC_2", "wAIC_3",
                                              "wAIC_4", "wAIC_5")],
                           MARGIN = 1, FUN = mean, na.rm = TRUE))

wAIC_null_model <- wAIC_table$mean_wAIC[1]

wAIC_table <- wAIC_table %>% 
  mutate(delta_wAIC = mean_wAIC - wAIC_table$mean_wAIC[1])


ggplot(data = wAIC_table[-1,], aes(x = effect, y = delta_wAIC)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ variable) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("07_results/01_interim_results/model_outputs/graph/delata_wAIC.png",
       width = 10, height = 8.5)




covariates <- unique(wAIC_table$variable)[-1]
for (k in 1:length(covariates)) {
  wAIC_table_k <- wAIC_table %>%
    filter(variable == covariates[k])
  best_effect_k
}


