# this script contains the function used to code the coevolution model. 
# in short, the model will simulate the evolution of a single trait of each 
# species in a community over time as a result of selection imposed by 
# interactions and the environment. the model will return a dataframe 
# with trait values at steady state and other metrics

coevolution_modelhl = function(bi , phi, alpha, theta, init, m, epsilon, t_max,sigma){
  # Simulate coevolutionary dynamics of a mutualistic network of species interactions
  
  # Arguments:
  # n_sp: total number of species in the network
  # bi : bipartite matrix representing the mutualistic-antagonistic network (interaction typ: 1--mutualism,-1--antagonism,0.5--neutral)
  # f : square adjacency matrix representing the mutualistic-antagonistic network (presence-absence)
  # h : vector of heritability values
  # alpha : parameter alpha, sensitivity of selection to trait matching
  # theta : vector of environmental optimum values
  # init : vector of intial trait values
  # m : vector of proportion of selection due to mutualism
  # epsilon : value to determine when equilibrium is reached
  # t_max : maximum number of timesteps allowed
  # sigma : critical trait mismatch of antagonistic interaction!
  
  # Obs:
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1],..., row_sp[nrow], col_sp[1],...,col_sp[ncol]))#first plant and animal
  #
  # Returns:
  # A matrix containing in reach row ti, the trait values (z) of all species at time t.
  #turn f0 into binary matrix f with all interactions are set as one
  n_row = nrow(bi)#number of preys
  n_col = ncol(bi)
  n_sp = n_row + n_col
  f = rbind(cbind(matrix(0,n_row,n_row),bi),cbind(t(bi),matrix(0,n_col,n_col)))
  f[f!=0]<-1
  # setup matrix
  z_mat = matrix(NA,nrow=t_max,ncol=n_sp) # matrix to store z values
  # initial trait values
  z_mat[1,] = init # initial trait values
  #to restore indirect effect,if indirect remains NA, means communties did no reach equilibrium
  z_final <- rep(NA, n_sp)
  indirect <- rep(NA, n_sp)
  indirect_col_sums <- rep(NA, n_sp)
  indirect_row_sums <- rep(NA, n_sp)
  indirect_col_sums_over_ind <- rep(NA, n_sp)
  indirect_row_sums_over_ind <- rep(NA, n_sp)
  indirect_col_row_sums <- rep(NA, n_sp)
  indirect_col_row_sums_over_ind <- rep(NA, n_sp)
  current_time <- NA
  rate <- NA
  # simulation for a give maximum time
  for (t in 1:(t_max - 1)){ # simulation runs for a maximum of t_max timesteps
    # current trait
    z = z_mat[t,] # current z values
    # calculate differences in traits
    z_dif = t(f*z) - f*z # matrix with all trait differences #animal's trait -plant's trait
    # calculate q matrix
    q = f*(exp(-alpha*(z_dif^2))) # calculating matrix q (evolutionary effect of species j on i)###key step
    # standardize q matrix
    q_n = q / apply(q,1,sum) # normalizing matrix q
    # calculate selection from mutualism
    q_m = q_n * m # multiplying each row i of matrix q by m[i]
    # obtain selection differentials
    #####separate antagonistic and neutral interaction here####
    z.change <- z_dif
    if(any(bi==-1)){#when there are antagonistic
      #first identify prey
      z_difsub <-z.change[1:n_row,(n_row+1):n_sp][bi==-1]
      z_difsub.1 <- z_difsub
      #1. judge whether prey's trait is larger than predator or vice versa
      z_difsub.1[z_difsub>0]<- z_difsub.1[z_difsub>0]-sigma#when predator is larger
      z_difsub.1[z_difsub<0]<- z_difsub.1[z_difsub<0]+sigma#when predator is smaller #what about equals? maybe matters
      #2. judge whether the absolute difference of predator and prey's trait is larger than sigma=0.5
      u<- 1*(abs(z_difsub)<sigma)
      #put them back with mutualisitic interaction
      z.change[1:n_row,(n_row+1):n_sp][bi==-1] <- u*z_difsub.1
    }
    if(any(bi==0.5)){#when there are neutral
      z.change[1:n_row,(n_row+1):n_sp][bi==0.5] <- 0
    }
    
    sel_dif = q_m * z.change # calculating selection differentials
    # calculate selection by mutualism
    r_mut = phi * apply(sel_dif, 1, sum) # response to selection related to mutualism
    # calculate selection by environment
    r_env = phi * (1 - m)* (theta - z) # response to selection related to the environment
    # update trait values
    z_mat[t+1,] = z + r_mut + r_env # updating z values
    # calculate trait differences
    dif = mean(abs(z - z_mat[t+1,])) # computing the mean difference between old and new z values
    # evaluate if differences are at threshold
    if(dif < epsilon) {# check to see if equilibrium has been reached
      #quantify indirect effects
      current_time <- t
      Psi <- diag(1-m)# diag 1-m
      I <- diag(1,nrow=n_sp)
      Tm = solve(I - q_m) %*% Psi
      dir_ind = (1 - I) * Tm
      ind = (1 - f) * dir_ind# relative contribution of indirect and direct effects
      indirect = sum(ind) / sum(dir_ind)
      indirect_col_sums = colSums(ind) / colSums(dir_ind)
      indirect_row_sums = rowSums(ind) / rowSums(dir_ind)
      indirect_col_row_sums = rowMeans(cbind(indirect_col_sums, indirect_row_sums))
      indirect_col_sums_over_ind = colSums(ind) / sum(ind)
      indirect_row_sums_over_ind = rowSums(ind) / sum(ind)
      indirect_col_row_sums_over_ind = rowMeans(cbind(indirect_col_sums_over_ind, indirect_row_sums_over_ind))
      rate = mean(abs((z_mat[t+1,]-z_mat[1,])/current_time))
      z_final = z_mat[t+1,]
      break
    }
  }
  
  # return a matrix with trait values over time
  z.indi <- rbind(theta = theta, init = init, z_final = z_final, indirect=indirect,  indirect_col_sums=indirect_col_sums,
                                       indirect_row_sums=indirect_row_sums,  indirect_col_row_sums=indirect_col_row_sums,
                                       indirect_col_sums_over_ind = indirect_col_sums_over_ind, indirect_row_sums_over_ind = indirect_row_sums_over_ind,
                                       indirect_col_row_sums_over_ind = indirect_col_row_sums_over_ind, time_to_equi = rep(current_time, length(indirect_col_sums)),
                                       rate = rep(rate, length(indirect_col_sums)))
  return(z.indi)
}


