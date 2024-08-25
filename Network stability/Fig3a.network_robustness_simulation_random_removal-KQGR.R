# raw contributor: Linwei Wu 
# Yuan, M. M. et al. Climate warming enhances microbial network complexity and stability. Nature Climate Change 11, 343-348, doi:10.1038/s41558-021-00989-9 (2021).

###### directely read in correlation matrix downloaded from SparCC. 

###### end of the two choices of correlation matrix ########
otutab<-read.table("net-otuc.txt",header = T,check.names = F,row.names=1,sep="\t")

#otutab<-read.table("net-raw.txt",header = T,check.names = F,row.names=1,sep="\t")

comm<-as.data.frame(t(otutab))

comm<-t(otutab)
sp.ra<-colMeans(comm)/40858  #relative abundance of each species; the 40858 was rarefied depth  


### run  three times for different glacial ecosystem
cormatrix <-read.table("network.cor-PRO.txt",header = T,check.names = F,row.names=1,sep="\t")
cormatrix <-read.table("network.cor-SUP.txt",header = T,check.names = F,row.names=1,sep="\t")
cormatrix <-read.table("network.cor-SUB.txt",header = T,check.names = F,row.names=1,sep="\t")

row.names(cormatrix)<-colnames(cormatrix)

diag(cormatrix)<-0    #no links for self-self    
sum(abs(cormatrix)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix[colSums(abs(cormatrix))>0,colSums(abs(cormatrix))>0]

sp.ra2 <- sp.ra[row.names(network.raw)]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

## robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
 
  sp.meanInteration<-colMeans(net.stength)

  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  #write.csv( net.Raw,"network pruned.csv")
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    #result<-c(remain.mean,remain.sd,remain.se)
    #names(result)<-c("remain.mean","remain.sd","remain.se")
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result <- cbind(result,remains)
    result
  })) 
}


Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)

Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))


currentdat<-dat1
### Integrate the results of 3 runs
rpro <- currentdat
rpro$type <- "PRO"
rsup <- currentdat
rsup$type <- "SUP"
rsub <- currentdat
rsub$type <- "SUB"
rdf <- rbind(rpro,rsub,rsup)


df <- rdf %>% dplyr::filter(Proportion.removed==0.5) %>% dplyr::select(2:4,208,ncol(rdf)) 
names(df)<-c("remain.mean","remain.sd","remain.se","weighted","type")

write.csv(df,"random_removal_result-0.5-n.csv")



dff$type <-  factor(dff$type, levels=c("SUP","SUB","PRO"))

