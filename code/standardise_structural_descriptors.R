# standardise nestedness values

# requires
library(bipartite)

# define functions

# null model 

ddnullad<-function(MATRIX){
  
  r<-dim(MATRIX)[1]
  c<-dim(MATRIX)[2]
  coldegreesprop<-(colSums(MATRIX>0))/(r-1)
  rowdegreesprop<-(rowSums(MATRIX>0))/(c-1)
  
  flag=0
  while (flag == 0) {
    
    
    #Fill up each matrix element probabilistically depending on the matrix dimensions and
    #degree distribution
    pmatrix <- 0.5* (array(rep(coldegreesprop,rep(r,c)), dim=c(r,c)) + array(rep(rowdegreesprop,c),dim=c(r,c)))
    diag(pmatrix) <- 0
    TEST<- 1* ( array(runif(r*c), dim=c(r,c)) < pmatrix ) 
    
    flag=1
    if (length(dim(TEST)) < 2) {flag=0}
    else result<-TEST
  }
  
  return(result)
}

# standardise nestedness
standardise_nestedness <- function(matrix){
  
  # remove weights
  matrix[matrix > 0] <- 1
  
  # remove empty rows/columns
  matrix <- bipartite::empty(matrix)
  
  # compute observed nestedness
  obs_nestedness <- nestedness_s(matrix)
  
  # compute nestedness for all randomised matrices
  nestedness_randomised_matrices <- unlist(replicate(100, nestedness_s(ddnullad(matrix))))
  
  # standardise the observed nestedness value
  standardised_nestedness <- (obs_nestedness - mean(nestedness_randomised_matrices))/sd(nestedness_randomised_matrices)
  
  # return
  return(standardised_nestedness)
}

# standardise modularity
standardise_modularity <- function(matrix){
  
  # remove weights
  matrix[matrix > 0] <- 1
  
  # remove empty rows/columns
  matrix <- bipartite::empty(matrix)
  
  # compute observed modularity
  net <- igraph::graph_from_incidence_matrix(matrix)
  obs_modularity <- igraph::cluster_louvain(net)
  obs_modularity <- igraph::modularity(obs_modularity)
  
  # compute modularity for all randomised matrices
  modularity_randomised_matrices <- unlist(replicate(100, igraph::modularity(igraph::cluster_louvain(igraph::graph_from_incidence_matrix(ddnullad(matrix))))))
  
  # standardise the observed modularity value
  standardised_modularity <- (obs_modularity - mean(modularity_randomised_matrices))/sd(modularity_randomised_matrices)
  
  # return
  return(standardised_modularity)
}


# nestedness function
nestedness_s <- function(M){ 
  
  # this code computes the nestedness of a given incident matrix M
  # according to the definition given in
  # Fortuna, M.A., et al.: Coevolutionary dynamics shape the structure of 
  # bacteria‐phage infection networks. Evolution 1001-1011 (2019).
  # DOI 10.1111/evo.13731
  
  # Make sure we are working with a matrix
  M <- as.matrix(M)
  # Binarize the matrix
  B <- as.matrix((M>0))
  class(B) <- "numeric"
  
  # Get number of rows and columns
  nrows <- nrow(B)
  ncols <- ncol(B)
  
  # Compute nestedness of rows
  nestedness_rows <- 0
  for(i in 1:(nrows-1)){
    for(j in (i+1): nrows){
      
      c_ij <- sum(B[i,] * B[j,])      # Number of interactions shared by i and j
      k_i <- sum(B[i,])               # Degree of node i
      k_j <- sum(B[j,])               # Degree of node j
      
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_rows <- nestedness_rows + o_ij
    }
  }
  
  # Compute nestedness of columns
  nestedness_cols <- 0
  for(i in 1: (ncols-1)){
    for(j in (i+1): ncols){
      
      c_ij <- sum(B[,i] * B[,j])      # Number of interactions shared by i and j
      k_i <- sum(B[,i])               # Degree of node i
      k_j <- sum(B[,j])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected.
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_cols <- nestedness_cols + o_ij
    }
  }
  
  # Compute nestedness of the network
  nestedness_val <- (nestedness_rows + nestedness_cols) / ((nrows * (nrows - 1) / 2) + 
                                                             (ncols * (ncols - 1) / 2))
  
  return(nestedness_val)
  
}

# read all mutualistic networks 
# change directory
setwd("~/mutualistic_antagonistic_indirect_effects/data/networks/mutualistic/")
# list all files in directory
mutualsitic_webname <- list.files()
mutualsitic_weball<- list()
for (i in 1:length(mutualsitic_webname)){
  tmp <- read.csv(mutualsitic_webname[i],header = T,row.names = 1)
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  mutualsitic_weball[[i]]<-tmp
  
  txtname<-gsub("csv","txt",mutualsitic_webname[i])
}

# read all antagonistic networks 
# change directory
setwd("~/mutualistic_antagonistic_indirect_effects/data/networks/antagonistic/")
# list all files in directory
antagonistic_webname <- list.files()
antagonistic_weball<- list()
for (i in 1:length(antagonistic_webname)){
  tmp <- read.csv(antagonistic_webname[i],header = T,row.names = 1)
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  antagonistic_weball[[i]]<-tmp
  
  txtname<-gsub("csv","txt",antagonistic_webname[i])
}

# standardise nestedness of mutualistic networks 
standardised_nested_values_mutualistic <- unlist(lapply(mutualsitic_weball, standardise_nestedness))

# standardise nestedness of antagonistic networks
standardised_nested_values_antagonistic <- unlist(lapply(antagonistic_weball, standardise_nestedness))

# standardise modularity of mutualistic networks 
standardised_modularity_values_mutualistic <- unlist(lapply(mutualsitic_weball, standardise_modularity))

# standardise modularity of antagonistic networks
standardised_modularity_values_antagonistic <- unlist(lapply(antagonistic_weball, standardise_modularity))

# load data with structure
load('~/mutualistic_antagonistic_indirect_effects/output/network_structure/structuretable.mut.RData')
load('~/mutualistic_antagonistic_indirect_effects/output/network_structure/structuretable.anta.RData')

# add standardised values to data
strtab$z_score_nestedness <- standardised_nested_values_mutualistic
strtab.anta$z_score_nestedness <- standardised_nested_values_antagonistic
strtab$z_score_modularity <- standardised_modularity_values_mutualistic
strtab.anta$z_score_modularity <- standardised_modularity_values_antagonistic

# store data
save(strtab, file = '~/mutualistic_antagonistic_indirect_effects/output/network_structure/structuretable.mut.RData')
save(strtab.anta, file = '~/mutualistic_antagonistic_indirect_effects/output/network_structure/structuretable.anta.RData')

