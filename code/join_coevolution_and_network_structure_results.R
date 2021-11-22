# Requires
require(dplyr)
require(data.table)

# Load processed data
specialist_antagonistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/antagonistic_results_specialist.csv')
specialist_antagonistic_data <- specialist_antagonistic_data[,2:ncol(specialist_antagonistic_data)]
generalist_antagonistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/antagonistic_results_generalist.csv')
generalist_antagonistic_data <- generalist_antagonistic_data[,2:ncol(generalist_antagonistic_data)]
random_antagonistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/antagonistic_results_random.csv')
random_antagonistic_data <- random_antagonistic_data[,2:ncol(random_antagonistic_data)]
specialist_mutualistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/mutualistic_results_specialist.csv')
specialist_mutualistic_data <- specialist_mutualistic_data[,2:ncol(specialist_mutualistic_data)]
generalist_mutualistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/mutualistic_results_generalist.csv')
generalist_mutualistic_data <- generalist_mutualistic_data[,2:ncol(generalist_mutualistic_data)]
random_mutualistic_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/mutualistic_results_random.csv')
random_mutualistic_data <- random_mutualistic_data[,2:ncol(random_mutualistic_data)]

# Bind data
complete_data <- rbind(specialist_antagonistic_data, generalist_antagonistic_data, random_antagonistic_data,
                       specialist_mutualistic_data, generalist_mutualistic_data, random_mutualistic_data)

# Change column names
colnames(complete_data) <- c('network_name', 'network_type', 'anta_ratio', 'sequence', 
                             'replica', 'specialist_or_interacts_with_one', 'is_specialist',
                             'interacts_with_specialist', 'network_indirect_effects', 'species_indirect_effects_row_sum',
                             'species_indirect_effects_row_sum_over_ind', 'species_indirect_effects_col_sum',
                             'species_indirect_effects_col_sum_over_ind', 'time_to_equilibrium',
                             'rate_of_adaptive_change', 'degree', 'sp_names', 'species_matching_all_species_include_diagonal',
                             'species_matching_all_species_no_diagonal', 'species_matching_only_interactions', 
                             'species_trait_change', 'species_theta_matching', 'network_trait_change', 'network_theta_change',
                             'network_matching_all_species_include_diagonal', 'network_matching_all_species_no_diagonal',
                             'network_matching_only_interactions','z_final','z_initial','z_theta')


# Load network structure data
network_structure_data <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_structure.csv')
network_structure_data <- network_structure_data[,2:ncol(network_structure_data)]

# JOIN STATEMENT
complete_data <- complete_data %>% left_join(network_structure_data)

# Write file
write.csv(complete_data, '~/mutualistic_antagonistic_indirect_effects/results/coevolution_network_structure_results.csv')
