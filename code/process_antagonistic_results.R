# this script takes the raw output of all the simulations performed using the coevolutionary model
# on the antagonistic networks and compute a series of metrics for further analyses. It will
# store the post-processed results into a new file.


# clean up
rm(list=ls())
# load packages
require(doSNOW)
require(parallel)

# define function that will take the raw output of the antagonistic coevolution model 
# and calculate the relevant metrics discussed in the paper
process_antagonistic_results <- function(anta_sequence, anta_ratio){
  
  # Parameters
  alpha <- 0.2
  n_sim <- 100
  type <- 'antagonistic_network'
  # Initialize empty object
  df <- NULL
  
  # Network directory
  network_directory <- '~/indirect_effects_transition_mutualism_antagonism/data/networks/antagonistic/'
  # List files in directory
  net_files <- list.files(network_directory)
  
  # Loop through networks
  for (i in 1:length(net_files)) {
    
    # Current network name
    network_name <- tools::file_path_sans_ext(net_files[i])
    # Read matrix
    mat <- as.matrix(read.csv(paste0(network_directory, net_files[i]), row.names = 1))
    # Turn binary
    mat[mat>0] <- 1
    # Get number of plants and animals
    num_plants <- nrow(mat)
    num_animals <- ncol(mat)
    # Get degree of plants and animals
    plant_degree <- rowSums(mat)
    animal_degree <- colSums(mat)
    degree_vector <- as.numeric(c(plant_degree, animal_degree))
    
    
    # Build adjacency matrix
    f <- rbind(cbind(matrix(0,num_plants,num_plants),mat),cbind(t(mat),matrix(0,num_animals,num_animals)))
    
    # Extract fraction of antagonisitic interactions
    a_char = gsub(".", "", as.character(anta_ratio), fixed = TRUE)
    
    # Replica loop
    for (replica in 1:n_sim) {
      # Read the simulation result for the current network with specific replica, strength of mutualistic selection and simulation replica
      result_data <- read.csv(paste0('~/indirect_effects_transition_mutualism_antagonism/output/antagonistic_to_mutualistic/all_networks_m_07_a_',a_char,
                                     '_',anta_sequence,'/',network_name,'/',network_name,'.csv_m_07_alpha02_phi05_theta0-10_sim',
                                     replica,'_sigma10.csv'),row.names = 1)
      
      # Read trait data
      z_theta <- as.numeric(result_data[1,])
      z_initial <- as.numeric(result_data[2,])
      z_final <- as.numeric(result_data[3,])
      
      # Check if network achieved equilibrium
      if (!any(is.na(z_final))) {
        
        # if so measure degree of trait matching metrics
        matching_diff <- as.matrix(dist(z_final,upper=TRUE,diag=TRUE))
        matching <- (exp(-alpha*((matching_diff)^2))) 
        matching_only_interactions <- matching * f
        species_matching_only_interactions <- rowMeans(replace(matching_only_interactions, matching_only_interactions == 0, NA), na.rm = TRUE)
        matching_without_diagonal <- replace(matching, matching == 1, NA)
        species_matching_all_species_no_diagonal <- rowMeans(matching_without_diagonal,na.rm = TRUE)
        network_matching_only_interactions <- mean(species_matching_only_interactions, na.rm = TRUE)
        network_matching_only_interactions <- rep(network_matching_only_interactions, length(species_matching_only_interactions))
        network_matching_all_species_no_diagonal <- mean(species_matching_all_species_no_diagonal, na.rm = TRUE)
        network_matching_all_species_no_diagonal <- rep(network_matching_all_species_no_diagonal, length(species_matching_only_interactions))
        species_trait_change <- abs(z_final - z_initial)
        species_theta_matching <- exp(-alpha*(z_final - z_theta)^2)
        network_trait_change <- mean(species_trait_change, na.rm = TRUE)
        network_theta_change <- mean(species_theta_matching, na.rm = TRUE)
        
      }
      # otherwise set as NAs, simulation will be discarded
      else{
        species_matching_all_species_no_diagonal <- rep(NA, num_plants + num_animals)
        species_matching_only_interactions <- rep(NA, num_plants + num_animals)
        network_matching_all_species_no_diagonal <- rep(NA, num_plants + num_animals)
        network_matching_only_interactions <- rep(NA, num_plants + num_animals)
        species_trait_change <- rep(NA, num_plants + num_animals)
        species_theta_matching <- rep(NA, num_plants + num_animals)
        network_trait_change <- rep(NA, num_plants + num_animals)
        network_theta_change <- rep(NA, num_plants + num_animals)
      }
      
      # Extract results on indirect effects, rate of change and time to equilibrium
      indirect_effects_network <- as.numeric(result_data[4,])
      species_indirect_effects_col_sum <- as.numeric(result_data[5,])
      time_to_equilirbium <- as.numeric(result_data[11,])
      rate <- as.numeric(result_data[12,])
      z_final_raw <- z_final
      z_initial_raw <- z_initial
      theta_raw <- z_theta
      row_sp_names <- rep('R',num_plants)
      col_sp_names <- rep('C', num_animals)
      sp_names <- c(row_sp_names, col_sp_names)
      
      # Put everything together
      out_df <- as.data.frame(cbind(rep(network_name,length(degree_vector)), 
                                    rep(type, length(degree_vector)),
                                    rep(anta_ratio,length(degree_vector)), 
                                    rep(anta_sequence, length(degree_vector)), 
                                    rep(replica, length(degree_vector)), 
                                    indirect_effects_network, 
                                    time_to_equilirbium,
                                    rate,
                                    degree_vector, 
                                    sp_names, 
                                    species_matching_all_species_no_diagonal, 
                                    species_matching_only_interactions, 
                                    species_trait_change, 
                                    species_theta_matching, 
                                    species_indirect_effects_col_sum,
                                    network_trait_change,
                                    network_theta_change,
                                    network_matching_all_species_no_diagonal, 
                                    network_matching_only_interactions,
                                    z_final_raw,
                                    z_initial_raw,
                                    theta_raw))
      # append data to data frame
      df <- rbind(df, out_df)
      
    }
    
  }
  
  # create directory to store output
  dir.create('~/indirect_effects_transition_mutualism_antagonism/results/', showWarnings = FALSE)
  # store the output in a csv file
  write.csv(df, file = paste0('~/indirect_effects_transition_mutualism_antagonism/results/antagonistic_results_', anta_sequence, '_', anta_ratio,'.csv'))
  return(0)
}

# create temporary data frame for all parameter combinations
strategies <- rep(c('generalist', 'specialist'), each = 6)
fractions <- rep(seq(0,1,0.2), 2)
parameter_df <- data.frame(strategies, fractions)


# Run post-processing of results in parallel 
# initialise cluster for parallel processing
cl =  makeCluster(nrow(parameter_df))
registerDoSNOW(cl)
# run the post-processing of raw data
result = foreach(i=1:nrow(parameter_df)) %dopar% 
  process_antagonistic_results(anta_sequence = parameter_df$strategies[i], 
                              anta_ratio = parameter_df$fractions[i])
# stop the cluster
stopCluster(cl)
