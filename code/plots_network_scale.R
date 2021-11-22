### Network scale plots

# Requires
require(data.table)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

# Load data
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]

##############################################
# ANTA RATIO
##############################################

####### TRAIT CHANGE ####### 
# Trait change across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_trait_change,
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('Trait change in network') + scale_fill_brewer(palette = "Oranges")

# Trait change across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_trait_change,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Trait change in network') + scale_fill_brewer(palette = "Accent")

# Trait change across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_trait_change,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Trait change in network') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)

####### TRAIT MATCHING ####### 
# Trait matching across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_matching_all_species,
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('Trait matching in network (all species)') + scale_fill_brewer(palette = "Oranges")

# Trait matching across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_matching_all_species,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Trait matching in network (all species)') + scale_fill_brewer(palette = "Accent")

# Trait matching across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_matching_all_species,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Trait matching in network (all species)') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)

####### INDIRECT EFFECTS ####### 
# Indirect effects across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_indirect_effects,
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('Network indirect effects') + scale_fill_brewer(palette = "Oranges")

# Indirect effects across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_indirect_effects,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Network indirect effects') + scale_fill_brewer(palette = "Accent")

# Indirect effects across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = network_indirect_effects,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Network indirect effects') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)


####### Time to equilibirum ####### 
# Time to equilibrium across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = log(time_to_equilibrium),
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('log Time to equilibrium') + scale_fill_brewer(palette = "Oranges")

# Time to equilibrium across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = log(time_to_equilibrium),
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('log Time to equilibrium') + scale_fill_brewer(palette = "Accent")

# Time to equilibrium across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = log(time_to_equilibrium),
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('log Time to equilibrium') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)


####### Rate of adaptive change ####### 
# Rate of adaptive change across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = rate_of_adaptive_change,
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('Rate of adaptive change') + scale_fill_brewer(palette = "Oranges")

# Rate of adaptive change across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = rate_of_adaptive_change,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Rate of adaptive change') + scale_fill_brewer(palette = "Accent")

# Rate of adaptive change across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = rate_of_adaptive_change,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Rate of adaptive change') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)

####### IS NA ####### 
# IS NA change across antagonism ratio
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = is.na,
           fill = as.factor(anta_ratio))) +
  geom_boxplot(alpha = 0.6) + geom_jitter(alpha = 0.5) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position = 'none', text = element_text(size=18)) +
  xlab('Fraction of antagonistic interactions') + ylab('Fraction of runs that did not converge') + scale_fill_brewer(palette = "Oranges")

# IS NA across antagonism ratio by network type
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = is.na,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Fraction of runs that did not converge') + scale_fill_brewer(palette = "Accent")

# IS NA change across antagonism ratio by network type and sequence
ggplot(data = network_scale_results,
       aes(x = as.factor(anta_ratio),
           y = is.na,
           fill = network_type)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Fraction of antagonistic interactions') + ylab('Fraction of runs that did not converge') + scale_fill_brewer(palette = "Accent") + facet_grid(~sequence)


##############################################
# NETWORK STRUCTURE
##############################################

# Indirect effects

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_indirect_effects)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Indirect effects') 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_indirect_effects,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Indirect effects') + 
  scale_colour_brewer(palette = "Accent") 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_indirect_effects,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Indirect effects') + 
  scale_colour_brewer(palette = "Accent") + facet_grid(sequence~as.factor(anta_ratio))


# Trait matching

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_matching_all_species)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Trait matching') 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_matching_all_species,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Trait matching') + 
  scale_colour_brewer(palette = "Accent") 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = network_matching_all_species,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Trait matching') + 
  scale_colour_brewer(palette = "Accent") + facet_grid(sequence~as.factor(anta_ratio))

# Rate of adaptive change

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = rate_of_adaptive_change)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Rate of adaptive change') 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = rate_of_adaptive_change,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Rate of adaptive change') + 
  scale_colour_brewer(palette = "Accent") 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = rate_of_adaptive_change,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Rate of adaptive change') + 
  scale_colour_brewer(palette = "Accent") + facet_grid(sequence~as.factor(anta_ratio))

# Is na

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = is.na)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Fraction of runs that did not converge') 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = is.na,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Fraction of runs that did not converge') + 
  scale_colour_brewer(palette = "Accent") 

ggplot(data = network_scale_results, 
       aes(x = pca_observed_values,
           y = is.na,
           col = network_type)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(aspect.ratio = 1, text = element_text(size=18), legend.position = 'bottom') +
  xlab('Network structure') + ylab('Fraction of runs that did not converge') + 
  scale_colour_brewer(palette = "Accent") + facet_grid(sequence~as.factor(anta_ratio))
