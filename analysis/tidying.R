##### tidying up!!

################################################################################
################ LIST OF INTERMEDIATE FILES TO REMOVE
indlist<-list()
spec<-c("Gorilla_beringei_beringei", "Gorilla_beringei_graueri","Gorilla_gorilla_diehli","Gorilla_gorilla_gorilla","Pan_paniscus","Pan_troglodytes_ellioti","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pongo_abelii","Pongo_pygmaeus","Pongo_tapanuliensis","captive_Pan","captive_Gorilla","captive_Pongo")
for (j in 1:length(spec)) {   indlist[[j]]<-list.dirs(paste("Result/Sequences/",spec[j],"",sep=""),recursive=F) }

cl<-list();dl<-list();i=0
for (j in (1:length(indlist))) {
  for ( k in (1:length(indlist[[j]]))) {
    i=i+1
    cl[[i]]<-c(system(paste("ls ",indlist[[j]][k],"/*.fastq.gz",sep=""),intern=T ),system(paste("ls ",indlist[[j]][k],"/*.fastq",sep=""),intern=T ));dl[[i]]<-c(system(paste("ls -d ",indlist[[j]][k],"/filtered*",sep="") ,intern=T),system(paste("ls -d ",indlist[[j]][k],"/mapping/bamread*",sep="") ,intern=T),system(paste("ls -d ",indlist[[j]][k],"/mapping/mark*",sep="") ,intern=T), system(paste("ls -d ",indlist[[j]][k],"/mapping/remove*",sep="") ,intern=T),system(paste("ls -d ",indlist[[j]][k],"/ERR*",sep="") ,intern=T),system(paste("ls -d ",indlist[[j]][k],"/SRR*",sep="") ,intern=T),system(paste("ls -d ",indlist[[j]][k],"/DRR*",sep="") ,intern=T)) 
  }
}

cl<-unlist(cl)
dl<-unlist(dl)
cl<-paste("rm ",cl,sep="")   
dl<-paste("rm -r ",dl,sep="")   

write.table(c("#!/bin/bash",cl,dl),file="Result/Sequences/tidy.txt",sep="\t",row.names=F,col.names=F,quote=F)
system(paste(dl))

################################################################################
################ SIZES OF DATASET & MD5SUM

basedir="summary_stats/"
indlist<-list()
spec<-c("Gorilla_beringei_beringei", "Gorilla_beringei_graueri","Gorilla_gorilla_diehli","Gorilla_gorilla_gorilla","Pan_paniscus","Pan_troglodytes_ellioti","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pongo_abelii","Pongo_pygmaeus","Pongo_tapanuliensis","captive_Pan","captive_Gorilla","captive_Pongo")
for (j in 1:length(spec)) {   indlist[[j]]<-list.dirs(paste("Result/Sequences/",spec[j],"",sep=""),recursive=F) }
options("scipen"=100)
load(file=paste(basedir,"covstats.Robject",sep=""))
indlist<-unlist(indlist)

gtd<-list();crd<-list();md5cram<-list();md5vcf<-list()
load(file=paste(basedir,"md5sums.Robject",sep=""))
for (i in (1:length(indlist))) {
    print(indlist[i])
    md5vcf[[i]]<-list();gtd[[i]]<-list()
    crd[[i]]<-system(paste("ls -l ",indlist[i],"/mapping/merge/*.merge.cram",sep=""),intern=T)
    md5cram[[i]]<-c(indlist[i],unlist(system(paste("md5sum ",indlist[i],"/mapping/merge/*.merge.cram",sep=""),intern=T)))
    for (chr in c(1:22,"X")) { 
      md5vcf[[i]][[chr]]<-system(paste("md5sum ",indlist[i],"/snpcalling/*_",chr,".g.vcf.gz",sep=""),intern=T) 
      gtd[[i]][[chr]]<-system(paste("ls -l ",indlist[i],"/snpcalling/*_",chr,".g.vcf.gz",sep=""),intern=T)
      }
    if(xv[i]=="M") { chr="Y"
      md5vcf[[i]][[chr]]<-system(paste("md5sum ",indlist[i],"/snpcalling/*_",chr,".g.vcf.gz",sep=""),intern=T)
      gtd[[i]][[chr]]<-system(paste("ls -l ",indlist[i],"/snpcalling/*_",chr,".g.vcf.gz",sep=""),intern=T)
      }
    md5vcf[[i]]<-unlist(md5vcf[[i]]);gtd[[i]]<-unlist(gtd[[i]])
    save(md5vcf,gtd,md5cram,crd,file=paste(basedir,"md5sums.Robject",sep=""))
  }


### get SRA IDs
basedir="summary_stats/"
load(file=paste(basedir,"covstats.Robject",sep=""))
indlist<-unlist(indlist)
sratab<-list()
for (i in 1:length(indlist)) {
  print(indlist[i])
  t1<-try(unlist(read.table(paste(indlist[i], "/SraAccList.txt",sep=""),sep="\t",as.is=T,fill=T)))
  if (class(t1)=="try-error") {
    l=0
    file=paste(indlist[i],"/mapping/merge/*.merge.cram",sep="")
    t1<-system(paste("samtools view -H ",file," | grep 'samtools merge' | tail -n 1",sep=""),intern=T)
    if(length(t1)==0) { t1<-system(paste("samtools view -H ",file," | grep 'samtools view' | tail -n 2 | head -n 1",sep=""),intern=T);l=1 }
    t1<-unlist(strsplit(t1,split=" "))
    t1<-t1[grep(".bam",t1)]
    t1<-t1[grep("merge",t1,invert=T)]
    t1<-unlist(do.call(rbind,strsplit(gsub("mapping/","",t1),split="/|\\."))[,ifelse(l==1,3,4)])
  }
  sratab[[i]]<-t1
}  

inds<-do.call(rbind,strsplit(indlist,split="/"))[,9]
sratabtab<-list()
for (j in (1:length(sratab))) { 
  sratabtab[[j]]<-c(cvl[j,1],inds[j],paste(sratab[[j]],collapse=", "))
}
sratabtab<-do.call(rbind,sratabtab)
write.table(sratabtab,file=paste(basedir,"sratab.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

## create per-individual list of files & md5sums (stratified by subspecies)
load(file=paste(basedir,"md5sums.Robject",sep=""))
i=0
ofil<-list()  
for (spe in 1:length(spec)) {
  for (j in (1:length(namlist[[spe]])))  {
    i=i+1
    ofil[[i]]<-cbind(do.call(rbind,strsplit(c(md5vcf[[i]],md5cram[[i]][2]),split="  /gpfs/|snpcalling/|/merge/"))[,c(3,1)], spec[[spe]], namlist[[spe]][j],do.call(rbind,strsplit(c(md5vcf[[i]],md5cram[[i]][2]),split="  "))[,c(2)]) }
}

otru<-c()
for (i in (1:length(ofil))) {
  otru[i]<-grepl(ofil[[i]][1,4],ofil[[i]][1,1])
  write.table(ofil[[i]],file=paste(basedir,"upload/",i,".txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
}

## Table for SI
fofil<-do.call(rbind,ofil)[,c(1:4)]
colnames(fofil)<-c("File","md5sum","Population","Name")
write.table(fofil,file=paste(basedir,"md5tab.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)


