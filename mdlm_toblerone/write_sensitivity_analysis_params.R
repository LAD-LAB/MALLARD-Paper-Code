library(dplyr)

setwd("/data/davidlab/users/jds/mdlm_publication_2/mdlm_toblerone/")


# Set up base parameters (the ones for the main analysis)
base <- within(list(), {
  c0 <- 25
  W_scale_mean <- 0
  W_scale_var <- 2
  V_scale_mean <- 1
  V_scale_var <- 2
  W_lkj <- 1
  V_lkj <- 1
})

# Set up Basic Structure
params.base <- data.frame(base)
params <- data.frame(base)

add_combo <- function(d, param, value){
  for (v in seq_along(value)){
    p <- params.base
    p[[param]] <- value[v]
    d <- rbind(d, p)
  }
  d
}

# Set up Combinations
params <- params %>% 
  add_combo("c0", c(10, 50)) %>% 
  add_combo("W_scale_mean", c(-.5, 2)) %>% 
  add_combo("W_scale_var", c(1, 3)) %>% 
  add_combo("V_scale_mean", c(-.5, 2)) %>% 
  add_combo("V_scale_var", c(1, 3)) %>% 
  add_combo("W_lkj", c(2,3)) %>% 
  add_combo("V_lkj", c(2,3))

write.table(params, file="sensitivity_analysis_params.txt")
