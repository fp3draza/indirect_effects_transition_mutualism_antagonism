# supp manuscript figures
#clean up
rm(list=ls())
dev.off()

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

# source functions to fit linear models
source('~/mutualistic_antagonistic_indirect_effects/code/fit_linear_models.R')

# load network data
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'mutualistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# fit models with indirect effects as response variable
modelfit_TM <- rbind(model_fit(network_scale_results_for_plots, "specialist first", "network_matching_all_species"), 
                     model_fit(network_scale_results_for_plots, "generalist first", "network_matching_all_species")) %>%
  mutate(sig=ifelse(p_value<0.05, "p < 0.05", "p ≥ 0.05"))

# set parameters for plots 
pd <- position_dodge(0.1)
color_pal <- c("#0072B2", "#D55E00")

# plot indirect effects
figure_S1 <- ggplot(data=modelfit_TM %>% mutate(sequence = if_else(anta_ratio %in% c(0,1), "a", sequence)), aes(x=anta_ratio)) + 
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_high, col=sequence), width=0, position=pd) +
  geom_point(aes(y=slope_estimate, col=sequence, shape=sig, size=R2), position=pd) +
  facet_wrap(~explanatory_var, scales="free") +
  scale_shape_manual(values=c(19, 1)) +
  xlim(c(0,1)) +
  xlab("\nFraction of antagonistic interactions") +
  ylab("Effect of structural descriptor on trait matching\n")  + theme_minimal() +
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                      values = c("#0072B2", "#D55E00", "red")) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic')) + guides(shape="none")

figure_S1
ggsave("figure_S1.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)


###################
# figure S2

full_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_mutualistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'mutualistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))
gradient_colors <- colorRampPalette(c("#D4AF37", "#19543E"))
summary_results_mutualistic <- full_results_mutualistic %>% group_by(network_name, sequence, anta_ratio, replica) %>% 
  summarise(network_trait_matching = mean(network_matching_all_species_no_diagonal),
            network_indirect_effects = mean(network_indirect_effects))

figure_S2 <- ggplot(data = summary_results_mutualistic, aes(x = network_trait_matching, y = network_indirect_effects, 
                                    col = network_name)) + theme_minimal() + facet_grid(~sequence) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'none',
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic')) + 
  xlab('\nNetwork trait matching') + ylab('Contribution of indirect effects\n to coevolution\n')   + 
  geom_smooth(se = TRUE) 

figure_S2
ggsave("figure_S2.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=200, units="mm", dpi=600)

###################
# figure S3

# Load data
full_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/coevolution_network_structure_results.csv')
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

figure_S3 <- ggplot(data = full_results_mutualistic_summary, aes(x = as.factor(anta_ratio), y = sd, col = as.factor(sequence), group = interaction(sequence,guild))) + 
  geom_point() + geom_line(aes(linetype = guild)) +
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12), legend.position = 'bottom')  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1) + ylab('Standard deviation of trait values\n after coevoluiton\n') + xlab('\nFraction of antagonistic interactions')

figure_S3
ggsave("figure_S3.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)


###################
# figure S4

# load network data
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'antagonistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))
network_scale_results_for_plots <- network_scale_results_for_plots %>% 
  mutate(
    sequence = ifelse(anta_ratio %in% c(0,1), 'a', sequence)
  )


# summarise data
full_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/coevolution_network_structure_results.csv')
full_results <- full_results[,2:ncol(full_results)]
full_results_antagonistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'antagonistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

full_results_antagonistic <- full_results_antagonistic %>% 
  mutate(
    sequence = ifelse(anta_ratio %in% c(0,1), 'a', sequence)
  )

# figure S4a
figure_S4a <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_indirect_effects, 
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Contribution of indirect effects\n to coevolution\n')  + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

# figure S4b
figure_S4b <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = rate_of_adaptive_change,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Rate of adaptive change\n')  + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

# figure S4c
figure_S4c <- ggplot(data = full_results_antagonistic, aes(x = as.factor(anta_ratio), 
                                                          y = z_final, fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7)  + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Trait values after coevolution\n') + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)


# figure S4d
figure_S4d <- ggplot(data = network_scale_results_for_plots,
                     aes(x = as.factor(anta_ratio),
                         y = network_matching_all_species,
                         fill = as.factor(sequence))) + geom_boxplot(alpha = 0.6) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Network trait matching\n') + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)


# merge figures
figure_S4 <- ggarrange(figure_S4a, figure_S4d, figure_S4b,figure_S4c,  # list of plots
                      labels = "AUTO", # labels
                      common.legend = T, # COMMON LEGEND
                      legend = "bottom", # legend position
                      nrow = 2,
                      ncol = 2,
                      hjust = -2)
annotate_figure(figure_S4, top = text_grob("simulations with plant-herbivore and host-parasite networks", 
                                      face = "bold", size = 14))
ggsave("figure_S4.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)

###################
# figure S5

# source functions to fit linear models
source('~/mutualistic_antagonistic_indirect_effects/code/fit_linear_models.R')

# load network data
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'antagonistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# fit models with indirect effects as response variable
modelfit_IE <- rbind(model_fit(network_scale_results_for_plots, "specialist first", "network_indirect_effects"), 
                     model_fit(network_scale_results_for_plots, "generalist first", "network_indirect_effects")) %>%
  mutate(sig=ifelse(p_value<0.05, "p < 0.05", "p ≥ 0.05"))

# set parameters for plots 
pd <- position_dodge(0.1)
color_pal <- c("#0072B2", "#D55E00")

# plot indirect effects
figure_S5 <- ggplot(data=modelfit_IE %>% mutate(sequence = if_else(anta_ratio %in% c(0,1), "a", sequence)), aes(x=anta_ratio)) + 
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_high, col=sequence), width=0, position=pd) +
  geom_point(aes(y=slope_estimate, col=sequence, shape=sig, size=R2), position=pd) +
  facet_wrap(~explanatory_var, scales="free") +
  scale_shape_manual(values=c(19, 1)) +
  xlim(c(0,1)) +
  xlab("\nFraction of antagonistic interactions") + ggtitle("simulations with plant-herbivore and host-parasite networks") +
  ylab("Effect of structural descriptor on indirect effects\n")  + theme_minimal() +
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                      values = c("#0072B2", "#D55E00", "red")) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic'),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + guides(shape="none")
figure_S5
ggsave("figure_S5.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)

###################
# figure S6

# Read data and wrangle
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'antagonistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

# fit models with indirect effects as response variable
modelfit_TM <- rbind(model_fit(network_scale_results_for_plots, "specialist first", "network_matching_all_species"), 
                     model_fit(network_scale_results_for_plots, "generalist first", "network_matching_all_species")) %>%
  mutate(sig=ifelse(p_value<0.05, "p < 0.05", "p ≥ 0.05"))

# set parameters for plots 
pd <- position_dodge(0.1)
color_pal <- c("#0072B2", "#D55E00")

# plot indirect effects
figure_S6 <- ggplot(data=modelfit_TM %>% mutate(sequence = if_else(anta_ratio %in% c(0,1), "a", sequence)), aes(x=anta_ratio)) + 
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_high, col=sequence), width=0, position=pd) +
  geom_point(aes(y=slope_estimate, col=sequence, shape=sig, size=R2), position=pd) +
  facet_wrap(~explanatory_var, scales="free") +
  scale_shape_manual(values=c(19, 1)) +
  xlim(c(0,1)) +
  xlab("\nFraction of antagonistic interactions") + ggtitle("simulations with plant-herbivore and host-parasite networks") +
  ylab("Effect of structural descriptor on trait matching\n")  + theme_minimal() +
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                      values = c("#0072B2", "#D55E00", "red")) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic'),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + guides(shape="none")

figure_S6

ggsave("figure_S6.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)


###################
# figure S7

figure_S7 <- ggplot(data = full_results_antagonistic_summary, aes(x = as.factor(anta_ratio), y = sd, col = as.factor(sequence), group = interaction(sequence,guild))) + 
  geom_point() + geom_line(aes(linetype = guild)) +
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12), legend.position = 'bottom')  + 
  scale_color_manual(values=color_pal) + theme(aspect.ratio = 1) + 
  ylab('Standard deviation of trait values\n after coevoluiton\n') + 
  xlab('\nFraction of antagonistic interactions') + ggtitle("simulations with plant-herbivore and host-parasite networks") + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic'),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) 

figure_S7
ggsave("figure_S7.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)

###################
# figure S8

# Load data
species_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/species_scale_results_summary.csv')
species_scale_results <- species_scale_results[,2:ncol(species_scale_results)]
network_scale_results <- fread('~/mutualistic_antagonistic_indirect_effects/results/network_scale_results_summary.csv')
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

species_scale_results_for_plots <- species_scale_results_for_plots %>% 
  mutate(
    sequence = ifelse(anta_ratio %in% c(0,1), 'a', sequence)
  )
  


## figure_S8
figure_S8 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = species_matching_all_species_no_diagonal, col = sequence))  + 
  scale_x_continuous(trans='log10')  +  geom_point(alpha = 0.3, shape = 21) + geom_smooth(method = 'lm') + theme_minimal() + theme(plot.subtitle = element_text(hjust = 0.5), legend.position = 'bottom', text = element_text(size=12), panel.margin.x=unit(1, "lines"), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species trait matching'), subtitle = 'Fraction of antagonistic interactions')  + 
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) 
figure_S8
ggsave("figure_S8.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)

## figure_S9
figure_S9 <- species_scale_results_for_plots %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1)) %>% ggplot(data = ., aes(x = degree, y = `indirect effects`, col = sequence))  + 
  scale_x_continuous(trans='log10') +  geom_point(alpha = 0.3, shape = 21)  + geom_smooth(method = 'lm') + theme_minimal() + theme(plot.subtitle = element_text(hjust = 0.5), legend.position = 'bottom', text = element_text(size=12), panel.margin.x=unit(1, "lines"), legend.title=element_blank()) +
  facet_grid(guild~anta_ratio) + labs(x = '\nDegree (log)', y = expression('Species contribution to indirect effects'), subtitle = 'Fraction of antagonistic interactions')  + 
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                      values = c("#0072B2", "#D55E00", "red")) 
figure_S9
ggsave("figure_S9.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)


# figures sensitivity analyses

# load data
sensitivity_analyses <- read.csv("~/mutualistic_antagonistic_indirect_effects/results/sensitivity_analyses.csv", row.names = 1)
sensitivity_analyses <- sensitivity_analyses %>% 
  mutate(anta_ratio = if_else(
    anta_ratio ==  0, 0, if_else(
      anta_ratio == 2, 0.2, if_else(
        anta_ratio == 4, 0.4, if_else(
          anta_ratio == 6, 0.6, if_else(
            anta_ratio == 8, 0.8, 1
          )
        )
      )
    )
  ),
  strategy = ifelse(strategy == 'generalist', 'generalist first', 'specialist first')) %>% 
  mutate(
    strategy = ifelse(anta_ratio %in% c(0,1), 'a', strategy)
  )


## figure S10
ggplot(data = sensitivity_analyses %>% dplyr::filter(sigma == 10, alpha == 0.2), 
       aes(x = as.factor(anta_ratio), y = indirect_effects, fill = strategy)) +
  geom_boxplot(alpha = 0.6) + facet_grid(~m) + xlab('\nFraction of antagonistic interactions') + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12),
                          aspect.ratio = 1, legend.position = 'bottom',
                          plot.subtitle = element_text(size = 10, face = 'bold', hjust = 0.5)) +
  ylab("Contribution of indirect effects\n to coevolution") +   
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + labs(subtitle = expression(m~values))
ggsave("figure_S10.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=180, units="mm", dpi=600)

## figure S11
ggplot(data = sensitivity_analyses %>% dplyr::filter(sigma == 10, m == 0.7),
       aes(x = as.factor(anta_ratio), y = indirect_effects, fill = strategy)) +
  geom_boxplot(alpha = 0.6) + facet_grid(~alpha) + xlab('\nFraction of antagonistic interactions') + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12),
                          aspect.ratio = 1, legend.position = 'bottom',
                          plot.subtitle = element_text(size = 10, face = 'bold', hjust = 0.5)) +
  ylab("Contribution of indirect effects\n to coevolution") +   
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + labs(subtitle = expression(alpha~values))
ggsave("figure_S11.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=180, units="mm", dpi=600)

## figure S12
ggplot(data = sensitivity_analyses %>% dplyr::filter(alpha == 0.2, m == 0.7),
       aes(x = as.factor(anta_ratio), y = indirect_effects, fill = strategy)) +
  geom_boxplot(alpha = 0.6) + facet_grid(~sigma) + xlab('\nFraction of antagonistic interactions') + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12),
                          aspect.ratio = 1, legend.position = 'bottom',
                          plot.subtitle = element_text(size = 10, face = 'bold', hjust = 0.5)) +
  ylab("Contribution of indirect effects\n to coevolution") +   
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + labs(subtitle = expression(epsilon~values))
ggsave("figure_S12.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=180, units="mm", dpi=600)
