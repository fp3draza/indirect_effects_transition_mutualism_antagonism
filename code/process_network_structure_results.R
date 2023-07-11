# this script joins the files containing the structural metrics of mutualistic and antagonistic interactions
# into a single file, which is stored in the .csv format

# load packages
library(dplyr)

# clean up
rm(list=ls())

# load data
load('~/indirect_effects_transition_mutualism_antagonism/output/network_structure/structuretable.mut.RData')
load('~/indirect_effects_transition_mutualism_antagonism/output/network_structure/structuretable.anta.RData')

# Bind data
network_structure <- rbind(structure_data_mutualistic, structure_data_antagonistic)

# Remove rownames
rownames(network_structure) <- NULL

# Tidy network name column and rename
network_structure <- network_structure %>% rename(network_name = File)
network_structure <- network_structure %>% mutate(network_name = tools::file_path_sans_ext(network_name))

# Tidy column names
network_structure <- network_structure %>% mutate(modularity_z = z_score_modularity,
                                                  nestedness_z = z_score_nestedness)

# Create directory where data will be stored
dir.create('~/indirect_effects_transition_mutualism_antagonism/results/', showWarnings = FALSE)

# Save data to file
write.csv(network_structure, '~/indirect_effects_transition_mutualism_antagonism/results/network_structure.csv')
