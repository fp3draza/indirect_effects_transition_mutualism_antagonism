# this script converts the mutualistic interactions of a network 
# into antagonistic interactions based on the degree of the species. 

###### function used to FILL -1 link in species with specified sequence accordance#####
fillfun <- function(mat,degreetab,anum,decreasing){
  degreeord <-order(degreetab$value,decreasing = decreasing)
  degreetab <- degreetab[degreeord,]
  degreetab$cum <- cumsum(degreetab$degree)
  if(degreetab$cum[1]<=anum){
    fillsp <- which(degreetab$cum<=anum)
    maxsp <- max(fillsp)
    fillid<-degreetab$id[fillsp]
    maxid1<-degreetab$id[maxsp+1]
    mat[,fillid][mat[,fillid]!=0]<- -1
    left<- anum-degreetab$cum[maxsp]
  }else {
    left <- anum
    maxsp<-0
    maxid1<-degreetab$id[maxsp+1]
  }
  if(left>0){
    mat[,maxid1][mat[,maxid1]!=0]<- sample(c(rep(-1,left),rep(1,degreetab[maxsp+1,"degree"]-left)))
  }
  return(mat)
}
