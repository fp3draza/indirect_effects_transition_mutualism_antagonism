# requires
require(dplyr)
require(FactoMineR)

# load data
load('~/indirect_effects_transition_mutualism_antagonism/output/network_structure/structuretable.mut.RData')
load('~/indirect_effects_transition_mutualism_antagonism/output/network_structure/structuretable.anta.RData')
load('~/indirect_effects_transition_mutualism_antagonism/output/network_structure/degreetab.RData')

# Bind data
network_structure <- rbind(strtab, strtab.anta)

# Remove rownames
rownames(network_structure) <- NULL
rownames(degreetab) <- NULL

# Tidy network name column and rename
network_structure <- network_structure %>% rename(network_name = File)
network_structure <- network_structure %>% mutate(network_name = tools::file_path_sans_ext(network_name))

# Join with degree data
degreetab <- degreetab %>% select(meandegree, vardegree)
network_structure <- cbind(network_structure, degreetab)

# Standardize nestedness and modularity
network_structure <- network_structure %>% mutate(modularity_s = (Modularity - mmean) / (mmean),
                                                  nestedness_s = ((obs.NODF - nmean) / (nmean)) * 0.01,
                                                  modularity_z = (Modularity - mmean) / msd,
                                                  nestedness_z = (obs.NODF - nmean) / nsd)


# Fit PCA 
pca_observed_values <- PCA(network_structure %>% select(nsp, connectance, obs.NODF, Modularity), graph = FALSE)
pca_observed_values_updated <- PCA(network_structure %>% select(nsp, connectance, obs_nestedness_fortuna, Modularity), graph = FALSE)
pca_mean_values <- PCA(network_structure %>% select(nsp, connectance, nmean, mmean), graph = FALSE)
pca_standardised_values <- PCA(network_structure %>% select(nsp, connectance, nestedness_s, modularity_s), graph = FALSE)
pca_zscore_values <- PCA(network_structure %>% select(nsp, connectance, nestedness_z, modularity_z), graph = FALSE)
pca_size_connectance <- PCA(network_structure %>% select(nsp, connectance), graph = FALSE)
pca_vardegree <- PCA(network_structure %>% select(nsp, connectance, obs.NODF, Modularity, vardegree), graph = FALSE)
pca_observerd_values_mutualistic <- PCA(network_structure[1:34,] %>% select(nsp, connectance, obs.NODF, Modularity), graph = FALSE)
pca_observerd_values_antagonistic <- PCA(network_structure[35:nrow(network_structure),] %>% select(nsp, connectance, obs.NODF, Modularity), graph = FALSE)
pca_vardegree_mutualistic <- PCA(network_structure[1:34,] %>% select(nsp, connectance, obs.NODF, Modularity, vardegree), graph = FALSE)
pca_vardegree_antagonistic <- PCA(network_structure[35:nrow(network_structure),] %>% select(nsp, connectance, obs.NODF, Modularity, vardegree), graph = FALSE)

# Add network coordinates on PCAs
network_structure <- network_structure %>% mutate(pca_observed_values = as.vector(pca_observed_values$ind$coord[,1]),
                             pca_observed_values_updated = as.vector(pca_observed_values_updated$ind$coord[,1]),
                             pca_mean_values = as.vector(pca_mean_values$ind$coord[,1]),
                             pca_standardised_values = as.vector(pca_standardised_values$ind$coord[,1]),
                             pca_zscore_values = as.vector(pca_zscore_values$ind$coord[,1]),
                             pca_size_connectance = as.vector(pca_size_connectance$ind$coord[,1]),
                             pca_observed_by_network_type = c(as.vector(pca_observerd_values_mutualistic$ind$coord[,1]), as.vector(pca_observerd_values_antagonistic$ind$coord[,1])),
                             pca_vardegree = as.vector(pca_vardegree$ind$coord[,1]),
                             pca_vardegree_by_network_type = c(as.vector(pca_vardegree_mutualistic$ind$coord[,1]), as.vector(pca_vardegree_antagonistic$ind$coord[,1])))
# Save data to file
write.csv(network_structure, '~/indirect_effects_transition_mutualism_antagonism/results/network_structure.csv')
