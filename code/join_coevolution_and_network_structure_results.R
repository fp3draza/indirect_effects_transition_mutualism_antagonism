# This scripts load all post-processed results and merges them into a single file.

# clean up
rm(list=ls())

# Load packages
require(dplyr)
require(data.table)
require(tidyverse)

# Define directory where data is stored 
dir_data <- '~/indirect_effects_transition_mutualism_antagonism/results/'

# read all files corresponding to the output of the coevolution simulations
list_mutualistic <- list.files(path = dir_data, pattern = "mutualistic*", full.names = TRUE)
list_antagonistic <- list.files(path = dir_data, pattern = "antagonistic*", full.names = TRUE)
list_full <- c(list_mutualistic, list_antagonistic)
complete_data <- list_full %>% 
  map_df(~fread(.))

# remove first column
complete_data <- complete_data[ , 2:ncol(complete_data)]

# Change column names
colnames(complete_data) <- c('network_name', 'network_type', 'anta_ratio', 'sequence', 
                             'replica', 'network_indirect_effects', 'time_to_equilibrium',
                             'rate_of_adaptive_change', 'degree', 'sp_names',
                             'species_matching_all_species_no_diagonal', 'species_matching_only_interactions',
                             'species_trait_change', 'species_theta_matching', 
                             'species_indirect_effects_col_sum', 'network_trait_change',
                             'network_theta_change', 'network_matching_all_species_no_diagonal',
                             'network_matching_only_interactions', 'z_final','z_initial','z_theta')


# Load network structure data
network_structure_data <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_structure.csv')
network_structure_data <- network_structure_data[,2:ncol(network_structure_data)]

# JOIN STATEMENT
complete_data <- complete_data %>% left_join(network_structure_data)

# Create directory where data will be stored
dir.create('~/indirect_effects_transition_mutualism_antagonism/results/', showWarnings = FALSE)

# Write file
write.csv(complete_data, '~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
