### filtering
## calculate 98% central coverage per individual & coverage distribution

require(ggplot2)
require(gridExtra)
library("dplyr")
options("scipen"=100)

basedir="summary_stats/"
jval=as.numeric(unlist(commandArgs(TRUE)))[1]

indlist<-list()
spec<-c("Gorilla_beringei_beringei", "Gorilla_beringei_graueri","Gorilla_gorilla_diehli","Gorilla_gorilla_gorilla","Pan_paniscus","Pan_troglodytes_ellioti","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pongo_abelii","Pongo_pygmaeus","Pongo_tapanuliensis","captive_Pan","captive_Gorilla","captive_Pongo")
for (j in 1:length(spec)) {   indlist[[j]]<-list.dirs(paste("Result/Sequences/",spec[j],"",sep=""),recursive=F) }

namlist<-list()
for (j in (1:length(spec))) {   namlist[[j]]<-do.call(rbind,strsplit(indlist[[j]],split="/"))[,9] }

#distat<-list();distatX<-list();covdist<-list()
#for (j in (1:length(spec))) {  distat[[j]]<-list();distatX[[j]]<-list(); covdist[[j]]<-list(); for (k in (1:length(namlist[[j]]))) { distat[[j]][[k]]<-list();distatX[[j]][[k]]<-list() } }

#for (j in (1:length(spec))) {
  j=jval
  print(spec[j])
  for (k in (1:length(namlist[[j]]))) {
    print(namlist[[j]][k])
    adp=matrix(1:501,ncol=1)
    for (chr in 1:23) {
      chr<-ifelse(chr==23,"X",chr)
      dp<-try(strsplit(system(paste('bcftools stats --threads 4 -s - ',indlist[[j]][k],"/snpcalling/",namlist[[j]][k],"_",chr,".g.vcf.gz |",' grep "^DP"',sep=""),intern=T),split="\t"),silent=T)
      dp<-do.call(rbind,dp)
      dp[,3]<-gsub(">500","501",dp[,3]);dd<-cbind(as.numeric(dp[,3]),as.numeric(dp[,6]))
      adp<-merge(adp,dd,by.x=1,by.y=1, all.x=T,all.y=F)
    }
    adp[is.na(adp)]<-0
    kadp<-rowSums(adp[,-c(1,ncol(adp))])  
    dst<-cumsum(kadp/sum(kadp))
    dstX<-cumsum(adp[,ncol(adp)]/sum(adp[,ncol(adp)],na.rm=T))

    load(file=paste(basedir,"covdists.Robject",sep=""))
    covdist[[j]][[k]]<-kadp
    distat[[j]][[k]]<-adp[c(max(which(dst<=0.01)),min(which(dst>=0.99))),1]
    distatX[[j]][[k]]<-adp[c(max(which(dstX<=0.01)),min(which(dstX>=0.99))),1]
    save(distat,distatX,covdist,file=paste(basedir,"covdists.Robject",sep=""))
  }
#}

Sys.time()

q()


## write as table
load(file=paste(basedir,"covdists.Robject",sep=""))

allcov<-list()
for (j in (1:length(spec))) {
  print(spec[j])
  dk<-do.call(rbind,distat[[j]])
  dk[,1]<-ifelse(dk[,1]<5,5,dk[,1])
  dy<-do.call(rbind,distatX[[j]])
  dy[,1]<-ifelse(dy[,1]<5,5,dy[,1])
  dy[,1]<-ifelse(is.na(dy[,1])==T,5,dy[,1])
  allcov[[j]]<-cbind(spec[j],namlist[[j]],dk,dy)
  }  

acov<-do.call(rbind,allcov)
colnames(acov)<-c("Spec","Name","minAut","maxAut","minX","maxX")
write.table(acov,file=paste(basedir,"coverage_cutoffs.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)


