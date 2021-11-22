###strucutres####
library(bipartite)
#read sd networks and ouput network in txt for MODULAR
#MODULAR parameters: Newman index, Simulated annealing, initial temp:2, cooling temp:1.05,iteration factor:1,nullnum=100
setwd("~/indirect_effects_transition_mutualism_antagonism/data/networks/mutualistic/")
webname<-list.files()
weball<- list()
for (i in 1:length(webname)){
  tmp <- read.csv(webname[i],header = T,row.names = 1)
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  weball[[i]]<-tmp
  
  txtname<-gsub("csv","txt",webname[i])
  #write.table(tmp,txtname,row.names = F,col.names = F)
}

#save(weball,file="~/anta-mutual/structures/weball.RData")
####run nestedness in the terminal
require(doSNOW)
require(parallel)
library(bipartite)
load(file="~/anta-mutual/structures/weball.RData")

ddnullad<-function(MATRIX,returnp=F){
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
  if(returnp){
    return(list(result,pmatrix))
  }else{
    return(result)
  }
}

runnested<- function(web,nullnum){
  web[web>0]<-1
  web <- bipartite::empty(web)
  #obs <- vegan::nestednodf(web)$statistic["NODF"]
  obs <- nestedness_s(web)[[3]]
  nulls<- c()
  for (j in 1:nullnum){
    nulltmp <- ddnullad(web)
    #nulls[j] <- vegan::nestednodf(nulltmp)$statistic["NODF"]
    nulls[j] <- nestedness_s(web)[[3]]
  }
  return(c(obs=obs,mean=mean(nulls),sd=sd(nulls)))
}

nestedness_s <- function(x){ # x is the adjacency matrix (binary)
  # nestedness of rows 
  nested_rows <- 0
  for(i in 1:nrow(x)){
    for(j in i:nrow(x)){
      if(j>i){
        shared <- sum(x[i,]*x[j,]) # sum of common interactions
        min_shared <- min(sum(x[i,]),sum(x[j,])) # min of the degrees
        nested_rows <- nested_rows+(shared/min_shared)
      }
    }
  }
  
  nestedness_rows <- nested_rows/(nrow(x)*(nrow(x)-1)/2)
  
  # nestedness of columns
  nested_columns <- 0
  for(i in 1:ncol(x)){
    for(j in i:ncol(x)){
      if(j>i){
        shared <- sum(x[,i]*x[,j]) # sum of common interactions
        min_shared <- min(sum(x[,i]),sum(x[,j])) # min of the degrees
        nested_columns <- nested_columns+(shared/min_shared)
      }
    }
  }
  
  
  nestedness_columns <- nested_columns/(ncol(x)*(ncol(x)-1)/2)
  
  # nestedness of the network
  
  nestedness_network <- (nested_rows+nested_columns)/((nrow(x)*(nrow(x)-1)/2)+(ncol(x)*(ncol(x)-1)/2))
  
  return(list(nestedness_rows, nestedness_columns, nestedness_network))
  
}


cores = detectCores()
cl =  makeCluster(10)
registerDoSNOW(cl)
webnodf = foreach(i=1:length(weball),.inorder = TRUE,.combine = rbind) %dopar% 
  runnested(weball[[i]],nullnum = 100)#100 null models for each web
stopCluster(cl)
save(webnodf,file="~/indirect_effects_transition_mutualism_antagonism/data/network/structures/nodf.RData")


##### nodf return
hist((webnodf[,1]-webnodf[,2])/webnodf[,2])

###deal with MODULAR output###
setwd("~/indirect_effects_transition_mutualism_antagonism/structures/modularity")
mofile <- list.files()
mofile<- mofile[-40]
#things to do when the order is incorrect
#split <- strsplit(mofile, "OUT_2M_SD_")
#split <- as.numeric(sapply(split, function(x) x <- sub(".txt", "", x[2])))
#mofile <- mofile[order(split)]
molist<- list()
for (i in 1:39){
  molist[[i]]<- read.table(paste0(mofile[i]),header = T)
}
molist[[1]]
motab <- read.table("OUT_MOD.txt",header = T,stringsAsFactors = F)
#as above
#split <- strsplit(motab$File, "M_SD_")
#split <- as.numeric(sapply(split, function(x) x <- sub(".txt", "", x[2])))
#motab <- motab[order(split),]

motab$mean<-sapply(molist,function(x)mean(x[,"Modularity"])) 
motab$sd<-sapply(molist,function(x)sd(x[,"Modularity"])) 

test <- cbind(webnodf,motab)
colnames(test)[c(2:3,9:10)] <- c("nmean","nsd","mmean","msd")
ggplot(test,aes(x=obs.NODF,y=Modularity))+geom_point()#negative correlated
ggplot(strtab,aes(x=connectance,y=nsp))+geom_point()#negative correlated
summary(lm(obs.NODF~Modularity,test))

##adding species num and connectance
test$nsp<-NA
test$nrow<- NA
test$ncol<- NA
test$connectance <- NA
for (i in 1:39){
  tmp<-weball[[i]]
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  test$ncol[i] <-ncol(tmp)
  test$nrow[i]<- nrow(tmp)
  test$nsp[i]<- nrow(tmp)+ncol(tmp)
  test$connectance[i]<- networklevel(tmp,index = "connectance")
}

strtab <- test
rm(test)
save(strtab,file="structuretable.RData")


################# run the simulation for host-parasite and plant-herbivore ################
####### cleaning HP web before running #######
antafile <- list.files()
antaweb<- list()
for (i in 1:length(antafile)){
 tmp <- read.csv(antafile[i],header = T,row.names = 1,stringsAsFactors = F)
 if("Num..of.hosts.sampled"%in%colnames(tmp)){
  denum  <- which(colnames(tmp)=="Num..of.hosts.sampled")
  tmp <- tmp[,-denum]
 }
 write.csv(tmp,file = antafile[i])
 antaweb[[i]]  <- tmp
}

####for anta-web####
#modular calculation
#MODULAR parameters: Newman, Simulated annealing, initial temp:2, cooling temp:1.05,iteration factor:1,nullnum=100
for (i in 1:length(antaweb)){
  tmp<-antaweb[[i]]
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  write.table(tmp,gsub("csv","txt",antafile)[i],row.names = F,col.names = F)
}

##run in the terminal
require(doSNOW)
require(parallel)
library(bipartite)
load(file="~/indirect_effects_transition_mutualism_antagonism/structures/antaweb.RData")

ddnullad<-function(MATRIX,returnp=F){
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
  if(returnp){
    return(list(result,pmatrix))
  }else{
    return(result)
  }
}

runnested<- function(web,nullnum){
  web[web>0]<-1
  web <- bipartite::empty(web)
  obs <- nestedness_s(web)[[3]]
  nulls<- c()
  for (j in 1:nullnum){
    nulltmp <- ddnullad(web)
    nulls[j] <- nestedness_s(web)[[3]]
  }
  return(c(obs=obs,mean=mean(nulls),sd=sd(nulls)))
}

cores = detectCores()
cl =  makeCluster(5)
registerDoSNOW(cl)
antanodf = foreach(i=1:length(antaweb),.inorder = TRUE,.combine = rbind) %dopar% 
  runnested(antaweb[[i]],nullnum = 100)
stopCluster(cl)
save(antanodf,file="~/indirect_effects_transition_mutualism_antagonism/structures/antanodf.RData")

### combine all them in the table
###deal with MODULAR output###
setwd("~/indirect_effects_transition_mutualism_antagonism/structures/modularity-anta")
mofile <- list.files()
mofile<- mofile[-56]
molistanta<- list()
for (i in 1:55){
  molistanta[[i]]<- read.table(paste0(mofile[i]),header = T)
}
molistanta[[1]]
motabanta <- read.table("OUT_MOD.txt",header = T,stringsAsFactors = F)

motabanta$mean<-sapply(molistanta,function(x)mean(x[,"Modularity"])) 
motabanta$sd<-sapply(molistanta,function(x)sd(x[,"Modularity"])) 

test <- cbind(antanodf,motabanta)
colnames(test)[c(2:3,9:10)] <- c("nmean","nsd","mmean","msd")
ggplot(test,aes(x=obs.NODF,y=Modularity))+geom_point()#negative correlated
ggplot(strtab.anta,aes(x=connectance,y=nsp))+geom_point()#negative correlated
summary(lm(obs.NODF~Modularity,test))

##adding species num and connectance
test$nsp<-NA
test$nrow<- NA
test$ncol<- NA
test$connectance <- NA
for (i in 1:55){
  tmp<-antaweb[[i]]
  tmp[tmp>0]<-1
  tmp <- empty(tmp)
  test$ncol[i] <-ncol(tmp)
  test$nrow[i]<- nrow(tmp)
  test$nsp[i]<- nrow(tmp)+ncol(tmp)
  test$connectance[i]<- networklevel(tmp,index = "connectance")
}

strtab.anta <- test
rm(test)
save(strtab.anta,file="structuretable.anta.RData")
