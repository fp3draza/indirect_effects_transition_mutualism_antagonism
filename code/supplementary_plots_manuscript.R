# supp manuscript figures

# requires
require(ggplot2)
require(data.table)
require(dplyr)
require(ggpubr)
library(lme4)
library(lmerTest)

# define color palette
color_pal <- c("#0072B2", "#D55E00")


###################
# figure S1
# load network data
network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'antagonistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# summarise data
full_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_mutualistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'antagonistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# figure S1a
figure_S1a <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_indirect_effects, 
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Contribution of indirect effects\n to coevolution\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)


# figure S1b
figure_S1b <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = rate_of_adaptive_change,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Rate of adaptive change\n') + scale_fill_manual(values=color_pal)  + theme(aspect.ratio = 1)

# figure S1c
figure_S1c <- ggplot(data = full_results_mutualistic, aes(x = as.factor(anta_ratio), 
                                                          y = z_final, fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7)  + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=16)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Trait values after coevolution\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)


# figure S1d
figure_S1d <- ggplot(data = network_scale_results_for_plots,
                     aes(x = as.factor(anta_ratio),
                         y = network_matching_all_species,
                         fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=15)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Network trait matching\n') + scale_fill_manual(values=color_pal) + theme(aspect.ratio = 1)


# merge figures
figure_S1 <- ggarrange(figure_S1a, figure_S1b, figure_S1c,figure_S1d,  # list of plots
                      labels = "AUTO", # labels
                      common.legend = T, # COMMON LEGEND
                      legend = "bottom", # legend position
                      nrow = 2,
                      ncol = 2,
                      hjust = -2)



###################
# figure S2

network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'mutualistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


figure_S2a <- ggplot(data = network_scale_results_for_plots, 
                     aes(x = pca_observed_values_by_network_type,
                         y = network_indirect_effects,
                         col = sequence)) + 
  geom_point() + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(text = element_text(size=16), legend.position = 'bottom', plot.subtitle = element_text(hjust = 0.5)) +
  xlab('\nNetwork structure') + ylab('Contribution of indirect effects\n to coevolution\n') + geom_smooth(method = "lm", se = FALSE) + scale_color_manual(values=color_pal) + 
  theme(legend.title=element_blank()) + labs(subtitle = 'Fraction of antagonistic interactions')

figure_S2b <- ggplot(data = network_scale_results_for_plots, 
                     aes(x = pca_observed_values_by_network_type,
                         y = network_matching_all_species,
                         col = sequence)) + 
  geom_point() + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(text = element_text(size=16), legend.position = 'bottom', plot.subtitle = element_text(hjust = 0.5)) +
  xlab('\nNetwork structure') + ylab('Network trait matching\n') + geom_smooth(method = "lm", se = FALSE) + scale_color_manual(values=color_pal) + 
  theme(legend.title=element_blank()) + labs(subtitle = 'Fraction of antagonistic interactions')


figure_S2 <- ggarrange(figure_S2a, figure_S2b, # list of plots
                       labels = "AUTO", # labels
                       common.legend = T, # COMMON LEGEND
                       legend = "bottom", # legend position
                       nrow = 2)


###################
# figure S3

# Read data and wrangle
network_scale_results <- read.csv("~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv",header = T,row.names = 1)
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'mutualistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

# Fit models
fitlmlist <- lmList(network_indirect_effects~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="specialist first"))
fitlmlist2 <- lmList(network_indirect_effects~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="generalist first"))
fitlmlist3<- lmList(network_matching_all_species~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="specialist first"))
fitlmlist4<- lmList(network_matching_all_species~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="generalist first"))


# Extract model coefficients
slopedat <- matrix(NA,nrow = 24,ncol = 6)
colnames(slopedat) <- c("mean","low","high","response","anta_ratio","sequence")
slopedat<- as_data_frame(slopedat)
slopedat[,"mean"]<- c(coef(fitlmlist)[,2],coef(fitlmlist2)[,2],coef(fitlmlist3)[,2],coef(fitlmlist4)[,2])
slopedat[,"low"]<- c(confint(fitlmlist)[,1,2],confint(fitlmlist2)[,1,2],confint(fitlmlist3)[,1,2],confint(fitlmlist4)[,1,2])
slopedat[,"high"]<- c(confint(fitlmlist)[,2,2],confint(fitlmlist2)[,2,2],confint(fitlmlist3)[,2,2],confint(fitlmlist4)[,2,2])
slopedat[,"response"]<- rep(c("indirect effect","trait matching"),each=12)
slopedat[,"anta_ratio"]<- rep(seq(0,1,by=.2),4)
slopedat[,"sequence"]<- rep(c("specialist first","generalist first"),each=6,times=2)

# Plot indirect effects
figure_S3a <- ggplot(subset(slopedat, response == "indirect effect"), aes(x = anta_ratio, y = mean, color = as.factor(sequence))) +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.1), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_minimal() +
  theme(legend.title=element_blank(), text = element_text(size=16)) +
  scale_x_continuous(breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("\nFraction of antagonistic interactions") + ylab("Effect of network structure on indirect effects\n")  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1)

# Plot trait matching
figure_S3b <- ggplot(subset(slopedat,response=="trait matching"), aes(x = anta_ratio, y = mean, color = as.factor(sequence))) +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.1), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_minimal() + 
  theme(legend.title = element_blank(), text = element_text(size=16)) + 
  scale_x_continuous(breaks = c(0, 0.2, 0.4,0.6, 0.8, 1)) +
  xlab("\nFraction of antagonistic interactions") + ylab("Effect of network structure on trait matching\n") + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1)


figure_S3 <- ggarrange(figure_S3a, figure_S3b, # list of plots
                      labels = "AUTO", # labels
                      common.legend = T, # COMMON LEGEND
                      legend = "bottom", # legend position
                      ncol = 1,
                      hjust = -2)

###################
# figure S4

network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

figure_S4a <- ggplot(data = network_scale_results_for_plots, 
                     aes(x = pca_observed_values,
                         y = network_indirect_effects,
                         col = sequence)) + 
  geom_point() + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(text = element_text(size=16), legend.position = 'bottom', panel.margin.x=unit(1, "lines"), plot.subtitle = element_text(hjust = 0.5)) +
  xlab('\nNetwork structure') + ylab('Contribution of indirect effects\n to coevolution\n') + geom_smooth(method = "lm", se = FALSE) + scale_color_manual(values=color_pal) + 
  theme(legend.title=element_blank()) + scale_x_continuous(breaks=c(-5,0,4)) + labs(subtitle = 'Fraction of antagonistic interactions')

figure_S4b <- ggplot(data = network_scale_results_for_plots, 
                     aes(x = pca_observed_values,
                         y = network_matching_all_species,
                         col = sequence)) + 
  geom_point() + facet_grid(~anta_ratio) + theme_minimal() + 
  theme(text = element_text(size=16), legend.position = 'bottom', panel.margin.x=unit(1, "lines"), plot.subtitle = element_text(hjust = 0.5)) +
  xlab('\nNetwork structure') + ylab('Network trait matching\n') + geom_smooth(method = "lm", se = FALSE) + scale_color_manual(values=color_pal) + 
  theme(legend.title=element_blank()) + scale_x_continuous(breaks=c(-5,0,4)) + labs(subtitle = 'Fraction of antagonistic interactions')


figure_S4 <- ggarrange(figure_S4a, figure_S4b, # list of plots
                       labels = "AUTO", # labels
                       common.legend = T, # COMMON LEGEND
                       legend = "bottom", # legend position
                       nrow = 2)

###################
# figure S5

# Read data and wrangle
network_scale_results <- read.csv("~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv",header = T,row.names = 1)
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

# Fit models
fitlmlist <- lmList(network_indirect_effects~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="specialist first"))
fitlmlist2 <- lmList(network_indirect_effects~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="generalist first"))
fitlmlist3<- lmList(network_matching_all_species~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="specialist first"))
fitlmlist4<- lmList(network_matching_all_species~pca_observed_values_updated|anta_ratio,subset(network_scale_results_for_plots,sequence=="generalist first"))

# Extract model coefficients
slopedat <- matrix(NA,nrow = 24,ncol = 6)
colnames(slopedat) <- c("mean","low","high","response","anta_ratio","sequence")
slopedat<- as_data_frame(slopedat)
slopedat[,"mean"]<- c(coef(fitlmlist)[,2],coef(fitlmlist2)[,2],coef(fitlmlist3)[,2],coef(fitlmlist4)[,2])
slopedat[,"low"]<- c(confint(fitlmlist)[,1,2],confint(fitlmlist2)[,1,2],confint(fitlmlist3)[,1,2],confint(fitlmlist4)[,1,2])
slopedat[,"high"]<- c(confint(fitlmlist)[,2,2],confint(fitlmlist2)[,2,2],confint(fitlmlist3)[,2,2],confint(fitlmlist4)[,2,2])
slopedat[,"response"]<- rep(c("indirect effect","trait matching"),each=12)
slopedat[,"anta_ratio"]<- rep(seq(0,1,by=.2),4)
slopedat[,"sequence"]<- rep(c("specialist first","generalist first"),each=6,times=2)

# Plot indirect effects
figure_S5a <- ggplot(subset(slopedat, response == "indirect effect"), aes(x = anta_ratio, y = mean, color = as.factor(sequence))) +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.1), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_minimal() +
  theme(legend.title=element_blank(), text = element_text(size=16)) +
  scale_x_continuous(breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("\nFraction of antagonistic interactions") + ylab("Effect of network structure on indirect effects\n")  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1)

# Plot trait matching
figure_S5b <- ggplot(subset(slopedat,response=="trait matching"), aes(x = anta_ratio, y = mean, color = as.factor(sequence))) +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.1),  alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_minimal() + 
  theme(legend.title = element_blank(), text = element_text(size=16)) + 
  scale_x_continuous(breaks = c(0, 0.2, 0.4,0.6, 0.8, 1)) +
  xlab("\nFraction of antagonistic interactions") + ylab("Effect of network structure on trait matching\n") + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1)


figure_S5 <- ggarrange(figure_S5a, figure_S5b, # list of plots
                       labels = "AUTO", # labels
                       common.legend = T, # COMMON LEGEND
                       legend = "bottom", # legend position
                       ncol = 1,
                       hjust = -2)

###################
# figure S6
# Load data
full_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_mutualistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'antagonistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"),
         guild =  if_else(sp_names == 'C', 'consumer', 'resource'))

full_results_mutualistic_summary <- full_results_mutualistic %>% group_by(anta_ratio, sequence, guild) %>%
  summarise(iqr = IQR(z_final), sd = sd(z_final))

# figure S6
figure_S6 <- ggplot(data = full_results_mutualistic_summary, aes(x = as.factor(anta_ratio), y = sd, col = as.factor(sequence), group = interaction(sequence,guild))) + 
  geom_point() + geom_line(aes(linetype = guild)) +
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=16), legend.position = 'bottom')  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1) + ylab('Standard deviation of trait values\n after coevoluiton\n') + xlab('\nFraction of antagonistic interactions')

###
# Figure S7 & S8

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
         network_type == 'antagonistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first")) %>%
  rename(`indirect effects` = species_indirect_effects_col_sum) 
species_scale_results_for_plots <- species_scale_results_for_plots %>% mutate(guild = if_else(grepl('C',species_name), 'consumer', 'resource'))


## figure_S7
figure_S7 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = species_matching_all_species_no_diagonal, col = sequence))  + 
  scale_x_continuous(trans='log10')  +  geom_point(alpha = 0.3, shape = 21) + geom_smooth(method = 'lm') + theme_minimal() + theme(plot.subtitle = element_text(hjust = 0.5), legend.position = 'bottom', text = element_text(size=15), panel.margin.x=unit(1, "lines"), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species trait matching'), subtitle = 'Fraction of antagonistic interactions') + scale_color_manual(values=color_pal) 

## figure_S8
figure_S8 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = `indirect effects`, col = sequence))  + 
  scale_x_continuous(trans='log10') +  geom_point(alpha = 0.3, shape = 21)  + geom_smooth(method = 'lm') + theme_minimal() + theme(plot.subtitle = element_text(hjust = 0.5), legend.position = 'bottom', text = element_text(size=15), panel.margin.x=unit(1, "lines"), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species contribution to indirect effects'), subtitle = 'Fraction of antagonistic interactions') + scale_color_manual(values=color_pal) 






