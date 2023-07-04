# manuscript figures

# requires
require(ggplot2)
require(data.table)
require(dplyr)
require(ggpubr)
require(lme4)
require(lmerTest)
require(reshape2)
require(paletteer)
require(forcats)

###################
# figure 3

# load network data
network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
                                    filter(sequence %in% c('generalist','specialist'),
                                    network_type == 'mutualistic_network') %>% 
                                    mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))

network_scale_results_for_plots <- network_scale_results_for_plots %>% 
  mutate(
    sequence = ifelse(anta_ratio %in% c(0,1), 'a', sequence)
  )

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

# modify sequence information to only show one treament when 
# fully mutualistic or antagonistic
full_results_mutualistic <- full_results_mutualistic %>% 
  mutate(
    sequence = ifelse(anta_ratio %in% c(0,1), 'a', sequence)
  )

# figure 3a
figure_3a <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_indirect_effects,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Contribution of indirect effects \nto coevolution\n') + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

# figure 3b
figure_3b <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = rate_of_adaptive_change,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Rate of adaptive change\n')  + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

  
# figure 3c
figure_3c <- ggplot(data = full_results_mutualistic, aes(x = as.factor(anta_ratio), 
                                            y = z_final, fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7)  + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Trait values after coevolution\n') +
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

# figure 3d
figure_3d <- ggplot(data = network_scale_results_for_plots,
                    aes(x = as.factor(anta_ratio),
                        y = network_matching_all_species,
                        fill = as.factor(sequence))) + geom_boxplot(alpha = 0.7) + 
  theme_minimal() + theme(legend.title=element_blank(), text = element_text(size=12)) +
  xlab('\nFraction of antagonistic interactions') + ylab('Network trait matching\n')  + 
  scale_fill_manual(breaks = c("generalist first", "specialist first"),
                    values = c("#0072B2", "#D55E00", "red")) + theme(aspect.ratio = 1)

# merge figures
figure_3 <- ggarrange(figure_3a, figure_3d, figure_3b, figure_3c, # list of plots
                      labels = "AUTO", # labels
                      common.legend = T, # COMMON LEGEND
                      legend = "bottom", # legend position
                      nrow = 2,
                      ncol = 2,
                      hjust = -2)
figure_3




###################

# figure 4

# source functions to fit linear models
source('~/indirect_effects_transition_mutualism_antagonism/code/fit_linear_models.R')

# load network data
network_scale_results <- fread('~/indirect_effects_transition_mutualism_antagonism/results/network_scale_results_summary.csv')
network_scale_results <- network_scale_results[,2:ncol(network_scale_results)]
network_scale_results_for_plots <- network_scale_results %>% 
  filter(sequence %in% c('generalist','specialist'),
         network_type == 'mutualistic_network') %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))


# fit models with indirect effects as response variable
modelfit_IE <- rbind(model_fit(network_scale_results_for_plots, "specialist first", "network_indirect_effects"), 
                     model_fit(network_scale_results_for_plots, "generalist first", "network_indirect_effects")) %>%
  mutate(sig=ifelse(p_value<0.05, "p < 0.05", "p ≥ 0.05"))

# set parameters for plots 
pd <- position_dodge(0.1)
color_pal <- c("#0072B2", "#D55E00")

# plot indirect effects
figure_4 <- ggplot(data=modelfit_IE %>% mutate(sequence = if_else(anta_ratio %in% c(0,1), "a", sequence)), aes(x=anta_ratio)) + 
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_high, col=sequence), width=0, position=pd) +
  geom_point(aes(y=slope_estimate, col=sequence, shape=sig, size=R2), position=pd) +
  facet_wrap(~explanatory_var, scales="free") +
  scale_shape_manual(values=c(19, 1)) +
  xlim(c(0,1)) +
  xlab("\nFraction of antagonistic interactions") +
  ylab("Effect of structural descriptor on indirect effects\n")  + theme_minimal() +
  scale_colour_manual(breaks = c("generalist first", "specialist first"),
                      values = c("#0072B2", "#D55E00", "red")) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic')) + guides(shape="none")

figure_4


###################

# figure 5

# define colours
gradient_colors <- colorRampPalette(c("#D4AF37", "#19543E"))
# reaload data
full_results_mutualistic <- full_results %>% na.omit() %>% 
  filter(network_type == 'mutualistic_network', sequence  %in% c('generalist','specialist')) %>% 
  mutate(sequence = if_else(sequence == 'generalist', "generalist first", "specialist first"))
# filter a particular network
data_figure_five <- full_results_mutualistic %>% filter(anta_ratio %in% c(0,0.2,0.4,0.6,0.8,1),
                                    network_name == 'M_SD_034') %>% 
  group_by(anta_ratio, sequence, replica) %>% 
  summarise(network_matching_all_species = mean(network_matching_all_species_no_diagonal),
            network_indirect_effects = mean(network_indirect_effects))
# build plot
figure_5 <- ggplot(data = data_figure_five, aes(x = as.numeric(network_matching_all_species), y = as.numeric(network_indirect_effects), 
    col = as.numeric(anta_ratio))) + geom_point(size = 2, alpha = 0.7) + theme_minimal() + facet_grid(~sequence) + 
  theme(aspect.ratio = 1, text = element_text(size=12),  legend.title=element_blank(), legend.position = 'top',
        strip.text.x = element_text(size = 12, face = 'italic'),
        strip.text.y = element_text(size = 12, face = 'italic'),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = 'italic')) + 
  xlab('\nNetwork trait matching') + ylab('Contribution of indirect effects\n to coevolution\n') +
  scale_color_gradientn(colours = gradient_colors(6)) + labs(subtitle = 'fraction of antagonistic interactions')
figure_5
ggsave("figure_5.pdf", path="~/mutualistic_antagonistic_indirect_effects/figures/amnat_revision/",
       width=180, height=220, units="mm", dpi=600)
                                                         

####################

# figure 6

# source
source("~/indirect_effects_transition_mutualism_antagonism/code/functions_for_trait_matching_figure.R")

# colours
pal_two <-  paletteer_d("rcartocolor::OrYel")

# run simulation
specialist_output <- run_simulation(0.2, 'specialist', 1, 1)
generalist_output <- run_simulation(0.2, 'generalist', 1, 1)

# split output
specialist_tmatrix <-specialist_output[[1]][[1]]
specialist_trait_matching <-specialist_output[[1]][[2]]
specialist_interaction_matrix <-specialist_output[[2]]
generalist_tmatrix <-generalist_output[[1]][[1]]
generalist_trait_matching <-generalist_output[[1]][[2]]
generalist_interaction_matrix <-generalist_output[[2]]

# rename specialist tmatrix
rownames(specialist_tmatrix)[rownames(specialist_tmatrix) == ""] <- paste0("R", 1:7)
rownames(specialist_tmatrix)[8:28] <- paste0("C", 1:21)
colnames(specialist_tmatrix) <- c(paste0("R", 1:7), paste0("C", 1:21))
# rename specialist trait matching matrix
colnames(specialist_trait_matching) <- c(paste0("R", 1:7), paste0("C", 1:21))
rownames(specialist_trait_matching) <- c(paste0("R", 1:7), paste0("C", 1:21))
# rename specialist interaction matrix
colnames(specialist_interaction_matrix) <- c(paste0("R", 1:7), paste0("C", 1:21))
rownames(specialist_interaction_matrix) <- c(paste0("R", 1:7), paste0("C", 1:21))

# convert specialist tmatrix to dataframe
specialist_tmatrix_dataframe <- setNames(melt(specialist_tmatrix), c('sp1', 'sp2', 'tmatrix_value'))
specialist_trait_matching_dataframe <- setNames(melt(specialist_trait_matching), c('sp1', 'sp2', 'trait_matching'))
specialist_interaction_dataframe <- setNames(melt(specialist_interaction_matrix), c('sp1', 'sp2', 'interaction'))

# join dataframes
specialist_dataframe <- specialist_tmatrix_dataframe %>% left_join(specialist_trait_matching_dataframe) %>% 
  left_join(specialist_interaction_dataframe) 

# rename generalist tmatrix
rownames(generalist_tmatrix)[rownames(generalist_tmatrix) == ""] <- paste0("R", 1:7)
rownames(generalist_tmatrix)[8:28] <- paste0("C", 1:21)
colnames(generalist_tmatrix) <- c(paste0("R", 1:7), paste0("C", 1:21))

# rename generalist trait matching matrix
colnames(generalist_trait_matching) <- c(paste0("R", 1:7), paste0("C", 1:21))
rownames(generalist_trait_matching) <- c(paste0("R", 1:7), paste0("C", 1:21))

# rename generalist interaction matrix
colnames(generalist_interaction_matrix) <- c(paste0("R", 1:7), paste0("C", 1:21))
rownames(generalist_interaction_matrix) <- c(paste0("R", 1:7), paste0("C", 1:21))

# convert specialist tmatrix to dataframe
generalist_tmatrix_dataframe <- setNames(melt(generalist_tmatrix), c('sp1', 'sp2', 'tmatrix_value'))
generalist_trait_matching_dataframe <- setNames(melt(generalist_trait_matching), c('sp1', 'sp2', 'trait_matching'))
generalist_interaction_dataframe <- setNames(melt(generalist_interaction_matrix), c('sp1', 'sp2', 'interaction'))

# join dataframes
generalist_dataframe <- generalist_tmatrix_dataframe %>% left_join(generalist_trait_matching_dataframe) %>%
  left_join(generalist_interaction_dataframe) 


# add information on strategy
specialist_dataframe$strategy <- 'specialist first'
generalist_dataframe$strategy <- 'generalist first'
df <- rbind(specialist_dataframe, generalist_dataframe)

# order for plots 
degree_sp1 <- df %>%  
  filter(strategy == 'specialist first') %>%
  group_by(sp1) %>% 
  summarise(
    degree_one = sum(abs(interaction)) 
  ) %>% arrange(degree_one)


# build order of species based on degree
order_degree_resources <- rev(as.character(degree_sp1$sp1[degree_sp1$sp1 %in% c('R1', 'R2' , 'R3',  'R4', 'R5', 'R6', 'R7')]))
order_degree_consumers <- as.character(degree_sp1$sp1[!degree_sp1$sp1 %in% c('R1', 'R2' , 'R3',  'R4', 'R5', 'R6', 'R7')])
order_degree <- as.factor(c(order_degree_resources, order_degree_consumers))

# determine which interactions are in the upper triangle
generalist_trait_matching_lower <- upper.tri(generalist_trait_matching[order_degree, order_degree], diag = TRUE)
# add names to upper triangle matrix
row.names(generalist_trait_matching_lower) <- as.character(order_degree)
colnames(generalist_trait_matching_lower) <- as.character(order_degree)
# convert to dataframe
generalist_trait_matching_lower <- setNames(melt(generalist_trait_matching_lower), c('sp1', 'sp2', 'true_is_lower'))

# order data based on degrees
df <- df %>% 
  arrange(factor(sp1, levels = order_degree), 
          factor(sp2, levels = order_degree))

# add information on whether interaction is in upper triangle
df <- df %>% left_join(generalist_trait_matching_lower)

# keep only data in upper triangle
df <- df %>%
  mutate(
    trait_matching = ifelse(true_is_lower == 1, trait_matching, NA),
    interaction = ifelse(true_is_lower == 1, interaction, NA))

# plots
figure_6 <- ggplot(data = df, aes(x = fct_inorder(sp1), y = fct_inorder(sp2), fill = trait_matching))  +
  geom_tile(alpha = 0.8) + facet_grid(~strategy) +  
  scale_fill_gradientn(colours = pal_two, na.value="#93B9ACFF") +
  xlab('') + ylab('')  + theme_minimal() +
  guides(fill=guide_legend(title="trait matching"))  + 
  theme(aspect.ratio = 1, legend.position = 'top',
        legend.title=element_blank(),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 10, face = 'italic'),
        strip.text.y = element_text(size = 10, face = 'italic'),
        plot.title = element_text(size = 12, face = 'bold'),
        plot.subtitle = element_text(size = 12, face = 'bold', hjust = 0.5),
        axis.text = element_text(size = 7, face = 'bold.italic'),
        axis.text.x = element_text(angle = 90),
        axis.ticks.x=element_blank(), 
        axis.ticks.y = element_blank())  +                  
  geom_point(aes(size=as.factor((abs(interaction))), col = as.factor(interaction))) +
  geom_point(aes(size=as.factor((abs(interaction)))), shape = 1, col = "black") +
  scale_size_manual(values=c(0,1.5),guide="none") + 
  scale_colour_manual(values=c("white","#A99CD9FF","black"),guide="none") + labs(subtitle = 'trait matching')
figure_6

