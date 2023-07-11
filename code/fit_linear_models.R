# this function fits a linear model to study the relationship
# between either indirect effects or trait matching and 
# each network descriptor. the functions returns a data frame
# containing the estimates of each model

model_fit <- function(data_results, flip_sequence, response_variable){
  
  # facilitate the splitting of data
  data_results <- data_results %>% ungroup()
  
  # fit linear models: effect of each structure measure on indirect effects or trait matching
  if(response_variable=="network_indirect_effects"){
    fitlmlist_size <- lmList(network_indirect_effects~log(size)|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_conn <- lmList(network_indirect_effects~connectance|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_mod <- lmList(network_indirect_effects~modularity_z|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_nest <- lmList(network_indirect_effects~nestedness_z|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
  }
  if(response_variable=="network_matching_all_species"){
    fitlmlist_size <- lmList(network_matching_all_species~log(size)|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_conn <- lmList(network_matching_all_species~connectance|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_mod <- lmList(network_matching_all_species~modularity_z|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
    fitlmlist_nest <- lmList(network_matching_all_species~nestedness_z|anta_ratio, subset(data_results, sequence==flip_sequence), pool = FALSE)
  }
  
  # extract slope estimates
  slope_est <- rbind(summary(fitlmlist_size)$coefficients[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                     summary(fitlmlist_conn)$coefficients[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                     summary(fitlmlist_mod)$coefficients[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                     summary(fitlmlist_nest)$coefficients[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)))
  
  # etract confidence intervals
  CI_low <- rbind(confint(fitlmlist_size)[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                  confint(fitlmlist_conn)[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                  confint(fitlmlist_mod)[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                  confint(fitlmlist_nest)[,1,2] %>% melt() %>% mutate(anta_ratio=rownames(.)))
  CI_high <- rbind(confint(fitlmlist_size)[,2,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                   confint(fitlmlist_conn)[,2,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                   confint(fitlmlist_mod)[,2,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                   confint(fitlmlist_nest)[,2,2] %>% melt() %>% mutate(anta_ratio=rownames(.)))
  
  # extract p-values
  p_val <- rbind(summary(fitlmlist_size)$coefficients[,4,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                 summary(fitlmlist_conn)$coefficients[,4,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                 summary(fitlmlist_mod)$coefficients[,4,2] %>% melt() %>% mutate(anta_ratio=rownames(.)),
                 summary(fitlmlist_nest)$coefficients[,4,2] %>% melt() %>% mutate(anta_ratio=rownames(.)))
  
  # extract R2
  R2_val <- c(as.numeric(summary(fitlmlist_size)$adj.r.squared),
              as.numeric(summary(fitlmlist_conn)$adj.r.squared),
              as.numeric(summary(fitlmlist_mod)$adj.r.squared),
              as.numeric(summary(fitlmlist_nest)$adj.r.squared))
  
  # output dataframe
  df_out <- data.frame(expand.grid(anta_ratio=as.numeric(unique(slope_est$anta_ratio)),
                                   explanatory_var=c("size(log)", "connectance", "modularity", "nestedness")),
                       slope_estimate=slope_est$value,
                       confint_low=CI_low$value,
                       confint_high=CI_high$value,
                       p_value=p_val$value,
                       R2=R2_val) %>%
    mutate(sequence=flip_sequence)
  
  return(df_out)
}
