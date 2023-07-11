# this script runs the coevolutionary model on a set of mutualistic networks under
# different ratios of mutualists-to-antagonists. It will store the output of each 
# simulation to a file. 
# clean up
rm(list=ls())
# library
library(doSNOW)
library(parallel)
library(bipartite)

run_simulations_on_HPC = function(m,a,seq,network_directory,realanta=F){
  
  # Load required functions
  require(dplyr)
  
  ##############################################################################################
  # DEFINE PARAMETERS
  # Set number of replicas
  n_sim = 100
  
  
  # Define alpha
  alpha = 0.2
  # Define phi
  phi_mean = 0.5
  phi_sd = 0.01
  # Define m sd
  m_sd = 0.01
  # Define environemtnal optimum range
  theta_min = 0
  theta_max = 10
  # define threshold to stop simulations
  epsilon = 0.00001
  # define max number of generations
  t_max = 10000
  #critical mismatch of predator-prey interaction
  sigma=10
  # Define folder to store results
  folder = '~/indirect_effects_transition_mutualism_antagonism/output/antagonistic_to_mutualistic/'
  
  ##############################################################################################
  # RUN LOOP FOR SIMULATIONS
  # Set current m_mean value
  m_mean = m
  
  # Define metadata
  alpha_char = gsub(".", "", as.character(alpha), fixed = TRUE)
  phi_car = gsub(".", "", as.character(phi_mean), fixed = TRUE)
  m_char = gsub(".", "", as.character(m_mean), fixed = TRUE)
  theta_min_char = gsub(".", "", as.character(theta_min), fixed = TRUE)
  theta_max_char = gsub(".", "", as.character(theta_max), fixed = TRUE)
  a_char = gsub(".", "", as.character(a), fixed = TRUE)
  
  ########################################################
  # Create folder to store results
  dir.create(path = paste(folder,"all_networks","_m_",m_char,"_a_",a_char,"_",seq,sep=""), recursive = TRUE, showWarnings = FALSE)
  
  # all network files names
  net_files = list.files(network_directory,recursive = TRUE)
  # remove file extension
  net_names = net_files
  
  #store indirect effect value
  indirecttab <- matrix(,nrow=n_sim,ncol=length(net_files))
  
  # Read parameter dataframe 
  parameters <- read.csv('~/indirect_effects_transition_mutualism_antagonism/data/parameters/antagonistic_parameters.csv', row.names = 1)
  
  ########################################################
  # Reading network
  
  for (i in 1:length(net_files)) {#
    
    # create directory to store results from current network
    dir.create(path = paste(folder, "all_networks", "_m_", m_char, "_a_",a_char,"_",seq,"/", tools::file_path_sans_ext(net_names[i]), sep = ""), showWarnings = FALSE)
    print(net_files[i])
    mat = as.matrix(read.csv(paste(network_directory,net_files[i],sep=""),row.names =  1))
    mat <- bipartite::empty(mat)
    #turn into binary
    mat[mat!=0]<- 1
    # Define the number of rows, columns and species
    n_row = nrow(mat)
    n_col = ncol(mat)
    n_sp = n_row + n_col
    
    #number of links that need to change
    if(a!=0){
      Inum<-sum(mat)
      anum<-round(a*Inum)
      anidegree <- colSums(mat)
      degreetab <- as.data.frame(cbind(id=1:n_col,degree=anidegree))
      
      # choose those links that would be changed to antagonsitic 
      if(seq=="specialist"){
        degreetab$value<- anidegree
        mat<-fillfun(mat,degreetab,anum,decreasing = F)
        
      }else if(seq=="generalist"){
        degreetab$value<- anidegree
        mat<-fillfun(mat,degreetab,anum,decreasing = T)
        
      }else if(seq=="lownested"){
        nestre<-bipartite::nestedcontribution(mat)
        degreetab$value<- nestre$`higher level`
        mat<-fillfun(mat,degreetab,anum,decreasing = F)
      }else if(seq=="highnested"){
        nestre<-bipartite::nestedcontribution(mat)
        degreetab$value<- nestre$`higher level`
        mat<-fillfun(mat,degreetab,anum,decreasing = T)
      }
    }
    # Build the square adjacency matrix f
    csv = T
    
    #######################################################
    # Simulations
    for (k in 1:n_sim) {
      
      if (a!=0&seq == 'random') {
        Iall <- c(rep(-1,anum),rep(1,Inum-anum))
        mat[mat!=0]<- sample(Iall)
      }
      
      # Filter parameter dataframe 
      filtered_parameters <- parameters %>% filter(network_name == net_files[i], replica == k)
      phi <- filtered_parameters %>% select(phi) %>% unlist() %>% as.vector()
      theta <- filtered_parameters %>% select(theta) %>% unlist() %>% as.vector()
      init <- filtered_parameters %>% select(init) %>% unlist() %>% as.vector()
      m <- filtered_parameters %>% select(m) %>% unlist() %>% as.vector()
      #if it is real anta, then the sign of interaction should change
      if(realanta){
        mat<- -1*mat
      }
      
      # run simulations
      z = coevolution_modelhl(bi=mat, phi = phi, alpha = alpha, theta = theta, init = init, m = m, epsilon = epsilon, t_max = t_max,sigma = sigma)
      # build data frames
      colnames(z) = c(paste("R",1:n_row,sep = ""), paste("C", 1:n_col, sep = ""))
      # create folder to store output
      dir.create(paste0(folder,"all_networks","_m_",m_char, "_a_",a_char,"_",seq,"/",tools::file_path_sans_ext(net_names[i])), showWarnings = FALSE)
      # save species inital and final trait values
      write.csv(z, file = paste(folder,"all_networks","_m_",m_char, "_a_",a_char,"_",seq,"/",tools::file_path_sans_ext(net_names[i]),"/",net_names[i],
                                "_m_",m_char,"_alpha",alpha_char,"_phi",phi_car,"_theta",theta_min_char
                                ,"-",theta_max_char,"_sim",k,"_sigma",sigma,".csv",sep=""))
      
      
    }
  }
}
# source necessary functions
source('~/indirect_effects_transition_mutualism_antagonism/code/fillfun.R')
source('~/indirect_effects_transition_mutualism_antagonism/code/coevolution_function.R')

# create temporary data frame for all parameter combinations
strategies <- rep(c('generalist', 'specialist'), each = 6)
fractions <- rep(seq(0,1,0.2), 2)
parameter_df <- data.frame(strategies, fractions)

# set directory where networks are stored
dir = "~/indirect_effects_transition_mutualism_antagonism/data/networks/antagonistic/"

# Run simulations for each strategy and each fraction in parallel 
# initialise cluster for parallel processing
cl =  makeCluster(nrow(parameter_df))
registerDoSNOW(cl)
# run the post-processing of raw data
result = foreach(i=1:nrow(parameter_df)) %dopar% 
  run_simulations_on_HPC(m = 0.7, a = parameter_df$fractions[i], seq = parameter_df$strategies[i],
                         network_directory = dir, realanta = F)
# stop the cluster
stopCluster(cl)


