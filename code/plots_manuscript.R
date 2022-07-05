# manuscript figures

# requires
require(ggplot2)
require(data.table)
require(dplyr)
require(ggpubr)
library(lme4)
library(lmerTest)

###################
# figure 2

# load network data
network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
                                    filter(sequence %in% c('generalist','specialist'),
                                    network_type == 'mutualistic_network') %>% 
                                    mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# define color palette
color_pal <- c("#0072B2", "#D55E00")

# Summarise data
trait_distribution_summary <- network_scale_results_for_plots %>% 
  group_by(anta_ratio, sequence) %>% 
  summarise(z_final_mean_summary = mean(z_final_mean),
            z_final_sd = mean(z_final_sd)) 

full_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_mutualistic <- full_results %>% na.omit() %>% 
                                        filter(network_type == 'mutualistic_network', sequence  %in% c('generalist','specialist')) %>% 
                                        mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))
  
# figure 2a
figure_2a <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_indirect_effects,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Contribution of indirect effects \nto coevolution\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)

# figure 2b
figure_2b <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = rate_of_adaptive_change,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Rate of adaptive change\n') + scale_fill_manual(values=color_pal)  + theme(aspect.ratio = 1)

  
# figure 2c
figure_2c <- ggplot(data = full_results_mutualistic, aes(x = as.factor(anta_ratio), 
                                            y = z_final, fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7)  + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Trait values after coevolution\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)

# figure 2d
figure_2d <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_matching_all_species,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Network trait matching\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)


# merge figures
figure_2 <- ggarrange(figure_2a, figure_2b, figure_2c, figure_2d, # list of plots
                      labels = "AUTO", # labels
                      common.legend = T, # COMMON LEGEND
                      legend = "bottom", # legend position
                      nrow = 2,
                      ncol = 2,
                      hjust = -2)

###################
# figure 3
# Load data
full_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_mutualistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'mutualistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"),
         guild =  if_else(sp_names == 'C', 'consumer', 'resource'))

full_results_antagonistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'antagonistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"),
         guild =  if_else(sp_names == 'C', 'consumer', 'resource'))

full_results_mutualistic_summary <- full_results_mutualistic %>% group_by(anta_ratio, sequence, guild) %>%
  summarise(iqr = IQR(z_final), sd = sd(z_final))

full_results_antagonistic_summary <- full_results_antagonistic %>% group_by(anta_ratio, sequence, guild) %>%
  summarise(iqr = IQR(z_final), sd = sd(z_final))

figure_3 <- ggplot(data = full_results_mutualistic_summary, aes(x = as.factor(anta_ratio), y = sd, col = as.factor(sequence), group = interaction(sequence,guild))) + 
  geom_point() + geom_line(aes(linetype = guild)) +
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=16), legend.position = 'bottom')  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1) + ylab('Standard deviation of trait values\n after coevoluiton\n') + xlab('\nFraction of antagonistic interactions')


###################
# figure 4 and 5
# Load data
species_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/species_scale_results_summary.csv')
species_scale_results <- species_scale_results[,2:ncol(species_scale_results)]
network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_slim <- network_scale_results %>% select(network_name, network_matching_all_species,
                                                               network_type, anta_ratio, sequence)

# Join data
species_scale_results <- species_scale_results %>% left_join(network_scale_results_slim)
species_scale_results_for_plots <- species_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'mutualistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first")) %>%
  rename(`indirect effects` = species_indirect_effects_col_sum) 
species_scale_results_for_plots <- species_scale_results_for_plots %>% mutate(guild = if_else(grepl('C',species_name), 'consumer', 'resource'))

# figure 4 and 5
figure_4 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = species_matching_all_species_no_diagonal, col = sequence))  + 
  scale_x_continuous(trans='log10')  +  geom_point(alpha = 0.3, shape = 21) + geom_smooth(method = 'lm') + theme_minimal() + theme(legend.position = 'bottom', text = element_text(size=15), plot.subtitle = element_text(hjust = 0.5), panel.margin.x=unit(1, "lines"), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species trait matching'), subtitle = 'Fraction of antagonistic interactions') + scale_color_manual(values=color_pal) 

figure_5 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = `indirect effects`, col = sequence))  + 
  scale_x_continuous(trans='log10')  +  geom_point(alpha = 0.3, shape = 21) + geom_smooth(method = 'lm') + theme_minimal() + theme(legend.position = 'bottom', text = element_text(size=15), panel.margin.x=unit(1, "lines"), plot.subtitle = element_text(hjust = 0.5), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species contribution to indirect effects'), subtitle = 'Fraction of antagonistic interactions') + scale_color_manual(values=color_pal) 

# new figure
new_figure <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = . %>% filter(m == 0.7, network == 'M_SD_034'), aes(x = as.numeric(trait_matching), y = as.numeric(indirect_effects), col = as.factor(anta_ratio))) + 
  geom_point() + theme_minimal() + facet_grid(~strategy) + theme(aspect.ratio = 1, text = element_text(size=18),  legend.title=element_blank()) + xlab('\nNetwork trait matching') + ylab('Contribution of indirect effects\n to coevolution\n')
