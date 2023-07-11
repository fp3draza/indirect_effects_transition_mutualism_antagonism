# clean up
rm(list=ls())

# Requires
library(data.table)
library(dplyr)

# Load data
complete_data <- fread('~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
complete_data <- complete_data[,2:ncol(complete_data)]

# Summarise network scale
complete_data_network_scale_summary <- complete_data %>% 
                                                        group_by(network_name, network_type, anta_ratio, sequence) %>%
                                                        mutate(is_na = if_else(is.na(network_indirect_effects), 1, 0)) %>% 
                                                        summarise(time_to_equilibrium = mean(time_to_equilibrium, na.rm = TRUE),
                                                                  size = mean(n_sp, na.rm = TRUE),
                                                                  modularity_z = mean(modularity_z, na.rm = TRUE),
                                                                  nestedness_z = mean(nestedness_z, na.rm = TRUE),
                                                                  connectance = mean(connectance, na.rm = TRUE),
                                                                  rate_of_adaptive_change = mean(rate_of_adaptive_change, na.rm = TRUE),
                                                                  network_indirect_effects = mean(network_indirect_effects, na.rm = TRUE),
                                                                  network_trait_change = mean(network_trait_change, na.rm = TRUE),
                                                                  network_theta_match = mean(network_theta_change, na.rm = TRUE),
                                                                  network_matching_all_species = mean(network_matching_all_species_no_diagonal, na.rm = TRUE),
                                                                  network_matching_interacting_species = mean(network_matching_only_interactions, na.rm = TRUE),
                                                                  is.na = sum(is_na)/length(is_na),
                                                                  z_final_mean = mean(z_final, na.rm = TRUE),
                                                                  z_final_sd = sd(z_final, na.rm = TRUE))


# Summarise species scale
complete_data_species_scale_summary <- complete_data %>% 
                                                group_by(network_name, network_type, anta_ratio, sequence, replica, sp_names) %>%
                                                mutate(species_name = paste0(sp_names,1:length(sp_names)),
                                                       is_na = if_else(is.na(network_indirect_effects), 1, 0)) 

complete_data_species_scale_summary_final <- complete_data_species_scale_summary %>% 
        group_by(network_name, network_type, anta_ratio, sequence, species_name) %>%
        mutate(is_na = if_else(is.na(network_indirect_effects), 1, 0)) %>% 
        summarise(is_na = mean(is_na, na.rm = TRUE),
                  degree = mean(degree, na.rm = TRUE),
                  species_matching_all_species_no_diagonal = mean(species_matching_all_species_no_diagonal, na.rm = TRUE),
                  species_matching_only_interactions = mean(species_matching_only_interactions, na.rm = TRUE),
                  species_trait_change = mean(species_trait_change, na.rm = TRUE),
                  species_indirect_effects_col_sum = mean(species_indirect_effects_col_sum, na.rm = TRUE),
                  species_theta_matching = mean(species_theta_matching, na.rm = TRUE),
                  z_final_mean = mean(z_final, na.rm = TRUE),
                  z_final_sd = sd(z_final, na.rm = TRUE))

# Store output
write.csv(complete_data_network_scale_summary, file = '~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
write.csv(complete_data_species_scale_summary_final, file = '~/indirect_effects_transition_mutualism_antagonism/results/species_scale_results_summary.csv')
