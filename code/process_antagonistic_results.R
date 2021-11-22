rm(list=ls())
require(doSNOW)
require(parallel)

process_antagonistic_results <- function(anta_sequence){
  
  # Parameters
  alpha <- 0.2
  n_sim <- 100
  type <- 'antagonistic_network'
  # Initialize empty df
  df <- NULL
  
  # Network directory
  network_directory <- '~/indirect_effects_transition_mutualism_antagonism/data/networks/antagonistic/'
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
    
    # Create empty vector to identify specialist plants or plants which interact with specialist animals
    plant_is_specialist_or_interacts_with <- rep(0, num_plants)
    plant_is_specialist <- rep(0, num_plants)
    plant_interacts_with_specialist <- rep(0, num_plants)
    
    # Identify specialist plants (lowest degree in network)
    plant_is_specialist_or_interacts_with[which(plant_degree == min(plant_degree))] <- 1
    plant_is_specialist[which(plant_degree == min(plant_degree))] <- 1
    
    # Identify plants that interact with specialist animals 
    if (is.null(nrow(mat[,which(animal_degree == min(animal_degree))]))) {
      plant_is_specialist_or_interacts_with[which(mat[,which(animal_degree == min(animal_degree))] != 0)] <- 1
      plant_interacts_with_specialist[which(mat[,which(animal_degree == min(animal_degree))] != 0)] <- 1
    }
    
    else{
      plant_is_specialist_or_interacts_with[which(rowSums(mat[,which(animal_degree == min(animal_degree))]) != 0)] <- 1
      plant_interacts_with_specialist[which(rowSums(mat[,which(animal_degree == min(animal_degree))]) != 0)] <- 1
    }
    
    # Create empty vector to identify specialist animals or animals which interact with specialist plants
    animal_is_specialist_or_interacts_with <- rep(0, num_animals)
    animal_is_specialist <- rep(0, num_animals)
    animal_interacts_with_specialist <- rep(0, num_animals)
    
    # Identify specialist animals (lowest degree in network)
    animal_is_specialist_or_interacts_with[which(animal_degree == min(animal_degree))] <- 1
    animal_is_specialist[which(animal_degree == min(animal_degree))] <- 1
    
    # Identify plants that interact with specialist animals 
    if (is.null(nrow(mat[which(plant_degree == min(plant_degree)),]))) {
      animal_is_specialist_or_interacts_with[which(mat[which(plant_degree == min(plant_degree)),] != 0)] <- 1
      animal_interacts_with_specialist[which(mat[which(plant_degree == min(plant_degree)),] != 0)] <- 1
    }
    
    else{
      animal_is_specialist_or_interacts_with[which(colSums(mat[which(plant_degree == min(plant_degree)),]) != 0)] <- 1
      animal_interacts_with_specialist[which(colSums(mat[which(plant_degree == min(plant_degree)),]) != 0)] <- 1
    }
    
    # Put vectors together
    specialist_or_interacts_with_one <- c(plant_is_specialist_or_interacts_with, animal_is_specialist_or_interacts_with)
    is_specialist <- c(plant_is_specialist, animal_is_specialist)
    interacts_with_specialist <- c(plant_interacts_with_specialist, animal_interacts_with_specialist)
    
    # Build adjacency matrix
    f <- rbind(cbind(matrix(0,num_plants,num_plants),mat),cbind(t(mat),matrix(0,num_animals,num_animals)))
    
    # Antagonistic ratio loop
    for (anta in seq(0,1,0.2)) {
      
      # Extract fraction of antagonisitic interactions
      a_char = gsub(".", "", as.character(anta), fixed = TRUE)
      
      # Replica loop
      for (replica in 1:n_sim) {
        
        # Read the simulation result for the current network with specific connectance, replica, strength of mutualistic selection and simulation replica
        result_data <- read.csv(paste0('~/indirect_effects_transition_mutualism_antagonism/output/antagonistic_to_mutualistic/all_networks_m_07_a_',a_char,
                                       '_',anta_sequence,'/',network_name,'/',network_name,'.csv_m_07_alpha02_phi05_theta0-10_sim',
                                       replica,'_sigma10.csv'),row.names = 1)
        
        # Read trait data
        z_theta <- as.numeric(result_data[1,])
        z_initial <- as.numeric(result_data[2,])
        z_final <- as.numeric(result_data[3,])
        
        # Check if network achieved equilirbium
        if (!any(is.na(z_final))) {
          
          matching_diff <- as.matrix(dist(z_final,upper=TRUE,diag=TRUE))
          matching <- (exp(-alpha*((matching_diff)^2))) 
          
          species_matching_all_species_include_diagonal <- rowMeans(matching,na.rm = TRUE)
          matching_only_interactions <- matching * f
          species_matching_only_interactions <- rowMeans(replace(matching_only_interactions, matching_only_interactions == 0, NA), na.rm = TRUE)
          matching_without_diagonal <- replace(matching, matching == 1, NA)
          species_matching_all_species_no_diagonal <- rowMeans(matching_without_diagonal,na.rm = TRUE)
          network_matching_all_species_include_diagonal <- mean(species_matching_all_species_include_diagonal, na.rm = TRUE)
          network_matching_all_species_include_diagonal <- rep(network_matching_all_species_include_diagonal, length(species_matching_only_interactions))
          network_matching_only_interactions <- mean(species_matching_only_interactions, na.rm = TRUE)
          network_matching_only_interactions <- rep(network_matching_only_interactions, length(species_matching_only_interactions))
          network_matching_all_species_no_diagonal <- mean(species_matching_all_species_no_diagonal, na.rm = TRUE)
          network_matching_all_species_no_diagonal <- rep(network_matching_all_species_no_diagonal, length(species_matching_only_interactions))
          species_trait_change <- abs(z_final - z_initial)
          species_theta_matching <- exp(-alpha*(z_final - z_theta)^2)
          network_trait_change <- mean(species_trait_change, na.rm = TRUE)
          network_theta_change <- mean(species_theta_matching, na.rm = TRUE)
        }
        
        else{
          species_matching_all_species_include_diagonal <- rep(NA, num_plants + num_animals)
          species_matching_all_species_no_diagonal <- rep(NA, num_plants + num_animals)
          species_matching_only_interactions <- rep(NA, num_plants + num_animals)
          network_matching_all_species_include_diagonal <- rep(NA, num_plants + num_animals)
          network_matching_all_species_no_diagonal <- rep(NA, num_plants + num_animals)
          network_matching_only_interactions <- rep(NA, num_plants + num_animals)
          species_trait_change <- rep(NA, num_plants + num_animals)
          species_theta_matching <- rep(NA, num_plants + num_animals)
          network_trait_change <- rep(NA, num_plants + num_animals)
          network_theta_change <- rep(NA, num_plants + num_animals)
        }
        
        # Define the final trait value of the species after coevolution
        indirect_effects_network <- as.numeric(result_data[4,])
        indirect_effects_row_sum <- as.numeric(result_data[6,])
        indirect_effects_col_sum <- as.numeric(result_data[5,])
        indirect_effects_col_row_sums <- as.numeric(result_data[7,])
        indirect_effects_col_sums_over_ind <- as.numeric(result_data[8,])
        indirect_effects_row_sums_over_ind <- as.numeric(result_data[9,])
        indirect_effects_col_row_sums_over_ind <- as.numeric(result_data[10,])
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
                                      rep(anta,length(degree_vector)), 
                                      rep(anta_sequence, length(degree_vector)), 
                                      rep(replica, length(degree_vector)), 
                                      specialist_or_interacts_with_one, 
                                      is_specialist, 
                                      interacts_with_specialist,
                                      indirect_effects_network, 
                                      indirect_effects_row_sum, 
                                      indirect_effects_row_sums_over_ind,
                                      indirect_effects_col_sum, 
                                      indirect_effects_col_sums_over_ind,
                                      time_to_equilirbium,
                                      rate,
                                      degree_vector, 
                                      sp_names, 
                                      species_matching_all_species_include_diagonal, 
                                      species_matching_all_species_no_diagonal, 
                                      species_matching_only_interactions, 
                                      species_trait_change, 
                                      species_theta_matching, 
                                      network_trait_change,
                                      network_theta_change,
                                      network_matching_all_species_include_diagonal, 
                                      network_matching_all_species_no_diagonal, 
                                      network_matching_only_interactions,
                                      z_final_raw,
                                      z_initial_raw,
                                      theta_raw))
        #colnames(out_df) <- NULL
        df <- rbind(df, out_df)
        
      }
      
    }
    
  }
  
  write.csv(df, file = paste0('~/indirect_effects_transition_mutualism_antagonism/results/antagonistic_results_', anta_sequence, '.csv'))
  return(0)
}


# Run function in parallel
sequences <- c('random','generalist','specialist')
cl =  makeCluster(3)
registerDoSNOW(cl)
result = foreach(i=seq(1,3,1)) %dopar% 
  process_antagonistic_results(sequences[i])
stopCluster(cl)
