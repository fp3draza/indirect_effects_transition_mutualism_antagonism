# Generate parameter values for simulations
network_directory <- '~/indirect_effects_transition_mutualism_antagonism/data/networks/mutualistic/'
net_files <- list.files(network_directory)
n_replica <- 100

generate_parameter_dataframe <- function(network_directory, net_files, n_replica){
  
    # Create empty dataframe 
    param_df <- NULL
    
    for (network in net_files) {
    print(network)

    
    # Read matrix
    mat <- as.matrix(read.csv(paste(network_directory,network,sep=""),row.names =  1))
    mat <- bipartite::empty(mat)
    
    #turn into binary
    mat[mat!=0]<- 1
    
    # Define the number of rows, columns and species
    n_row = nrow(mat)
    n_col = ncol(mat)
    n_sp = n_row + n_col
    
    # Define alpha
    alpha = 0.2
    # Define phi
    phi_mean = 0.5
    phi_sd = 0.01
    # Define m sd
    m_mean = 0.7
    m_sd = 0.01
    # Define environemtnal optimum range
    theta_min = 0
    theta_max = 10
    
    
    for (replica in 1:n_replica) {
      
      # sample phi values
      phi = rnorm(n_sp, phi_mean, phi_sd)
      while(any(phi < 0)){
        phi = rnorm(n_sp, phi_mean, phi_sd)
      }
      
      # sample mutualism selection values
      m = rnorm(n_sp, m_mean, m_sd)
      while(any(m < 0)){
        m = rnorm(n_sp, m_mean, m_sd)
      }
      
      # sample theta
      theta = runif(n_sp, min = theta_min, max = theta_max)
      
      # sample initial tra  it values
      init = runif(n_sp, min = theta_min, max = theta_max)
      
      # Create df to append
      df_append <- data.frame(cbind(phi, m, theta, init))
      
      # Add network and replica information
      df_append$network_name <- network
      df_append$replica <- replica
      
      # Append to full data frame
      param_df <- rbind(param_df, df_append)
      
    }
    }
  
  # Save file
  write.csv(param_df, '~/indirect_effects_transition_mutualism_antagonism/data/parameters/mutualistic_parameters.csv')
}

# Run function
generate_parameter_dataframe(network_directory, net_files, n_replica)
