# raw contributor: Linwei Wu

otutab<-read.table("net-otuc.txt",header = T,check.names = F,row.names=1,sep="\t")

comm<-t(otutab)
sp.ra<-colMeans(comm)/40858  #relative abundance of each species


###### end of the two choices of correlation matrix ########
cormatrix <-read.table("network.cor-PRO.txt",header = T,check.names = F,row.names=1,sep="\t")
cormatrix <-read.table("network.cor-SUP.txt",header = T,check.names = F,row.names=1,sep="\t")
cormatrix <-read.table("network.cor-SUB.txt",header = T,check.names = F,row.names=1,sep="\t")

cormatrix2<-cormatrix*(abs(cormatrix)>=0.50)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sum(colSums(abs(cormatrix2))>0)
tt <- sp.ra[row.names(cormatrix2)]
sp.ra2<-tt[rowSums(abs(cormatrix2))>0]
#sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sp.ra3 <- sp.ra[row.names(network.raw)]
summary(sp.ra2==sp.ra3)

sum(row.names(network.raw)==names(sp.ra2))  #check if matched

sp.ra2 <- sp.ra[row.names(network.raw)]

sum(row.names(network.raw)==names(sp.ra2))  #check if matched



## robustness simulation 
#input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
#return the proportion of species remained

#get the keystone species list
#node.attri<-read.csv("input_file/NodeAttribute_Y14_W.txt",header = T, sep="\t")
#module.hub<-as.character(node.attri$Name[node.attri$Zi > 2.5 & node.attri$Pi <= 0.62])


#consider cascade effects: removed species will further influence the remaining nodes

rand.remov2.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

### random removel #####
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
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
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
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    #result<-c(remains,remain.mean,remain.sd,remain.se)
    #names(result)<-c("remain","remain.mean","remain.sd","remain.se")
    result
    #result <- cbind(result,remains)
  }))
}



### target removel #####

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=length(module.hub))
                 )

currentdat = dat1
write.csv(currentdat,"targeted_deletion_results.csv")

rsub <- currentdat
rsub$type <- 'sub'
rsup <- currentdat
rsup$type <- 'sup'
rdf <- rbind(rsub,rsup)



