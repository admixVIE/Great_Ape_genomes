#!/usr/bin/r

require(ggplot2)
require(gridExtra)
library("dplyr")
'%ni%' <- Negate('%in%')
options("scipen"=100)

# mosdepth?
mo=1

# bcfstats?
bc=1

# where to save the files?
basedir="summary_stats/"

## names
indlist<-list()
spec<-c("Gorilla_beringei_beringei", "Gorilla_beringei_graueri","Gorilla_gorilla_diehli","Gorilla_gorilla_gorilla","Pan_paniscus","Pan_troglodytes_ellioti","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pongo_abelii","Pongo_pygmaeus","Pongo_tapanuliensis","captive_Pan","captive_Gorilla","captive_Pongo")
for (j in 1:length(spec)) {   indlist[[j]]<-list.dirs(paste("Result/Sequences/",spec[j],"",sep=""),recursive=F) }
ptl<-unlist(indlist)
write.table(ptl,file=paste(basedir,"pathlist.txt",sep=""),row.names=F,col.names=F,quote=F)
capt<-read.table(paste("~/Great_Ape/genome_pipeline/files/captive_panel.txt",sep=""),as.is=T,header=F,sep="\t")

namlist<-list()
for (j in (1:length(spec))) {   namlist[[j]]<-do.call(rbind,strsplit(indlist[[j]],split="/"))[,9] }

typ=list()
for (j in 1:length(spec)) { typ[[j]]<-rep(spec[j],length(namlist[[j]])) }
typCo<-unlist(typ)
#typ[[13]]<- capt[match(namlist[[j]],capt[,1]),2]
typ=unlist(typ)
typC=unlist(typCo)
typ3<-typ
typ3<-gsub("_"," ",typ3);typC<-gsub("_"," ",typCo);specs<-gsub("_"," ",spec)
typ2=unlist(namlist)

mpty=c(" ","  ","   ","    ","     ","      ","       ","        ","         ","          ","           ","            ","             ","              ","               ","                ","                 ","                  ","                   ")

################################################################################
### MOSDEPTH coverage
if (mo==1) {

  ### mosdepth average coverage across genome
  print("mosdepth")
  covlist<-list()
  for (j in (1:length(spec))) {
    cvl<-rep(NA,length(indlist[[j]]))
    for (k in (1:length(indlist[[j]]))) {
      ccc<-try(read.table(paste(indlist[[j]][k],"/mapping/merge/mosdepth.mosdepth.summary.txt",sep=""),as.is=T,header=T,sep="\t"),silent=T)
      if (class(ccc) != "try-error") { cvl[k]<-ccc[grep("total",ccc[,1]),4] }
    }
    covlist[[j]]<-cbind(namlist[[j]],cvl)
  }
  
  write.table(do.call(rbind,covlist),paste(basedir,"mosdepths_coverage.tsv",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  
  ### mosdepth coverage per chromosome
  chrcov<-list()
  for (chr in c(1:22, "X","Y","M")) { chrcov[[chr]]<-list() }
  
  for (j in (1:length(spec))) {
    for (chr in 1:25) { chr<-ifelse(chr==23,"X",ifelse(chr==24,"Y",ifelse(chr==25,"M",chr)));  chrcov[[chr]][[j]]<-rep(NA,length(indlist[[j]])) }
    
    for (k in (1:length(indlist[[j]]))) {
      ccc<-try(read.table(paste(indlist[[j]][k],"/mapping/merge/mosdepth.mosdepth.summary.txt",sep=""),as.is=T,header=T,sep="\t"),silent=T)
      if (class(ccc) != "try-error") {  for (chr in 1:25) { chr<-ifelse(chr==23,"X",ifelse(chr==24,"Y",ifelse(chr==25,"M",chr)));  chrcov[[chr]][[j]][k]<-ccc[which(ccc[,1]==paste("chr",chr,sep="")),4]  } }
    }
  }

  ################
  ## Statistics ##
  
  ## genetic sexing
  cvl<-do.call(rbind,covlist)
  xv1=as.numeric(as.character(unlist(chrcov[["Y"]])))/as.numeric(as.character(unlist(chrcov[["X"]])))
  xv2=as.numeric(as.character(unlist(chrcov[["X"]])))/as.numeric(as.character(unlist(chrcov[[1]])))
  xv<-ifelse(xv2<0.75 & xv1>0.1,"M","F")
  tx<-table(xv);  tx/(sum(tx))

  ## metadata table
  typ4<-gsub("Gorilla_beringei_beringei","GBB",typCo)
  typ4<-gsub("Gorilla_beringei_graueri","GBG",typ4)
  typ4<-gsub("Gorilla_gorilla_diehli","GGD",typ4)
  typ4<-gsub("Gorilla_gorilla_gorilla","GGG",typ4)
  typ4<-gsub("Pan_paniscus","PPA",typ4)
  typ4<-gsub("Pan_troglodytes_ellioti","PTE",typ4)
  typ4<-gsub("Pan_troglodytes_schweinfurthii","PTS",typ4)
  typ4<-gsub("Pan_troglodytes_troglodytes","PTT",typ4)
  typ4<-gsub("Pan_troglodytes_verus","PTV",typ4)
  typ4<-gsub("Pongo_abelii","PA",typ4)
  typ4<-gsub("Pongo_pygmaeus","PP",typ4)
  typ4<-gsub("Pongo_tapanuliensis","PT",typ4)
  typ4<-gsub("captive_","CP_",typ4)
    
  cvl<-cbind(typCo,typ4,typ2,cvl[,2],xv)
  colnames(cvl)<-c("Species","Species_label","Name","Average_coverage","Sex")
  write.table(cvl,file=paste("~/Great_Ape/genome_pipeline/files/metadata.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  
  ## per-chromosome coverage metadata table
  filist<-list()
  for (chr in (1:length(chrcov))) { filist[[chr]]<-unlist(chrcov[[chr]]) }
  filist<-cbind(typ2,typ3,do.call(cbind,filist))  
  colnames(filist)<-c("Name","Species",paste("Coverage_chr",c(1:22,"X","Y","M"),sep=""))
  filist[,26]<-ifelse(xv=="F",filist[,26],NA)
  write.table(filist,file=paste(basedir,"coverage_per_chrom.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

  ## list of males
  mls<-ptl[which(xv=="M")]
  write.table(mls,file=paste(basedir,"pathlist_M.txt",sep=""),row.names=F,col.names=F,quote=F)
  
  ###########
  ## PLOTS ##
  
  ### plotting of average coverage
  cr=data.frame(typ=typC,tes=as.numeric(do.call(rbind,covlist)[,2]))
  cr<-cr[!is.na(cr[,2]),]
  pG<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=13),axis.text.x= element_text(size=14,angle=45,hjust=1,vjust=1.3),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15))  +xlab("")+ylab("") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=max(tes)*1.1))    + ggtitle(label="coverage")
  
  pdf(paste(basedir,"average_depth.pdf",sep=""),width=16,height=8)
  pG
  dev.off()
  
  ## plotting per chromosome
  ap<-list()
  sxp<-list();ct=0
  for (chr in c(1:22, "X","Y","M")) { 
    cr=data.frame(typ=typ3,tes=as.numeric(unlist(chrcov[[chr]])),sex=xv)
    cr<-cr[!is.na(cr[,2]),]
    ## subset by sex for chrX, chrY
    if (chr %in%c("X","Y")) {
      for (sx in c("F","M")) {
      ct=ct+1
      sxp[[ct]]<-cr %>% filter(cr$sex==sx) %>% ggplot(aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1) , axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15),plot.margin=margin(0.1,0.1,0.1,20) ) +xlab("")+ylab("") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=max(tes)*1.1)) + ggtitle(label=paste(ifelse(sx=="F","Females, ","Males, "),"coverage chr",chr,sep=""))
      }      
    }
    else 
    {
    p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1), axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,vjust=1.3,size=15), plot.margin=margin(0.1,0.1,0.1,20)) +xlab("")+ylab("") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=max(tes)*1.1))  + ggtitle(label=paste("coverage chr",chr,sep=""))
    ap[[chr]]<-p1
      }
    }

  # sex chromosomes
  pdf(paste(basedir,"average_depth_per_sex_chrom.pdf",sep=""),width=20,height=10)
  grid.arrange(sxp[[1]],sxp[[2]],ap[["M"]],sxp[[4]],nrow=2,ncol=2)
  dev.off()

  # autosomes
  pdf(paste(basedir,"average_depth_per_chrom.pdf",sep=""),width=24,height=50)
  grid.arrange(ap[[1]],ap[[2]],ap[[3]],ap[[4]],ap[[5]],ap[[6]],ap[[7]],ap[[8]],ap[[9]],ap[[10]],ap[[11]],ap[[12]],ap[[13]],ap[[14]],ap[[15]],ap[[16]],ap[[17]],ap[[18]],ap[[19]],ap[[20]],ap[[21]],ap[[22]],nrow=11,ncol=2)
  dev.off()
  
  
  ## plotting per individual
  chco<-list()
  for (chr in (1:22)) { chco[[chr]]<-unlist(chrcov[[chr]]) }
  chco<-do.call(rbind,chco)
  
  ai<-list();va=0
  for (spe in (1:length(specs))) { 
    cpa=ifelse(specs[[spe]]=="captive Pan",1,0)
    nave<-rep(unlist(namlist[[spe]]),22)
    chsub<-as.numeric(unlist(t(chco[,which(typC==specs[spe])])))
    pm<-ceiling(length(namlist[[spe]])/20)
    ain<-1:length(namlist[[spe]])
    if(cpa==1) { hc<-which(nave%in%c("Ai","Akira","Ayumu"));kain<-ain[which(namlist[[spe]]%ni%c("Ai","Akira","Ayumu"))];pm=3 }
    cel=ifelse(grepl("capt",specs[spe]),90,55)
    for (j in (1:pm)) {
      va=va+1
      if (cpa==1) { if (j<3) { suse<-namlist[[spe]][kain][which(kain>(j*20)-20 & kain<=(j*20))]  } else { suse<-namlist[[spe]][-kain];cel=200 } } else { suse<-namlist[[spe]][which(ain>(j*20)-20 & ain<=(j*20))]  }
      cr=data.frame(typ=rep(suse,22),tes=chsub[which(nave%in%suse)],cel=cel)
      cr<-cr[!is.na(cr[,2]),]
      if (length(suse)<20) { lvs<-unique(cr$typ); nr<-data.frame(typ=mpty[1:(20-length(suse))],tes=NA,cel=cel);cr<-rbind(cr,nr);cr$typ<-factor(cr$typ, levels=c(lvs,mpty[1:(20-length(suse))])) }
      p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none", axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45), axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) +xlab("")+ylab("")+  geom_segment(aes(x=0.0,y=0,xend=0.0,yend=cel*1.1)) + ggtitle(label=paste("coverage ",specs[[spe]],", part ",j,"/",pm,sep=""))
      ai[[va]]<-p1
      }
    }  
  
  ## plot per chrom per individual
  pdf(paste(basedir,"average_depth_per_individual.pdf",sep=""),width=24,height=57)
  grid.arrange(ai[[1]],ai[[2]],ai[[3]],ai[[4]],ai[[5]],ai[[6]],ai[[7]],ai[[8]],ai[[9]],ai[[10]],ai[[11]],ai[[12]],ai[[13]],ai[[14]],ai[[15]],ai[[16]],ai[[17]],ai[[18]],ncol=1)
  dev.off()

  ## save data
  save(cvl,typ3,typC,filist,indlist,namlist,covlist,chrcov,xv,xv1,xv2,file=paste("~/Great_Ape/genome_pipeline/files/covstats.Robject",sep=""))
  
  }

    

################################################################################
### bcftools stats values

#if (bc==1) {

  print("bcftools stats")
  sustat<-list()
  for (chr in c(1:22, "X","Y")) { sustat[[chr]]<-list() }
  
  for (j in (1:length(spec))) {
    print(spec[j])
    for (chr in 1:24) { chr<-ifelse(chr==23,"X",ifelse(chr==24,"Y",chr));  sustat[[chr]][[j]]<-matrix(NA,nrow=5,ncol=length(indlist[[j]])) }
    for (k in (1:length(namlist[[j]]))) {
      print(namlist[[j]][k])
      lapos<-try(read.table(paste(indlist[[j]][k],"/snpcalling/lastpos.txt",sep=""),sep="\t",as.is=T,header=F))
      for (chr in 1:24) {
        chrx<-chr
        chr<-ifelse(chr==23,"X",ifelse(chr==24,"Y",chr))
        tstv<-try(unlist(strsplit(system(paste('grep "^TSTV" ',indlist[[j]][k],"/snpcalling/bcftoolsstats_",namlist[[j]][k],"_",chr,".txt",sep=""),intern=T),split="\t"))[8],silent=T)
        sis<-try(unlist(strsplit(system(paste('grep "^SiS" ',indlist[[j]][k],"/snpcalling/bcftoolsstats_",namlist[[j]][k],"_",chr,".txt",sep=""),intern=T),split="\t"))[c(4,7)],silent=T)
        norec<-try(unlist(strsplit(system(paste('grep "number of records" ',indlist[[j]][k],"/snpcalling/bcftoolsstats_",namlist[[j]][k],"_",chr,".txt",sep=""),intern=T),split="\t"))[5],silent=T)
        sustat[[chr]][[j]][1,k]<-ifelse(inherits(lapos,"try-error"),NA,lapos[which(lapos[,1]==paste("chr",chr,sep="")),2])
        sustat[[chr]][[j]][2,k]<-ifelse(length(norec)==0,NA,norec)
        if (length(sis)>0) { sustat[[chr]][[j]][3:4,k]<-sis }
        sustat[[chr]][[j]][5,k]<-ifelse(length(tstv)==0,NA,tstv)
        }
      }
    }
  save(sustat,file=paste("~/Great_Ape/genome_pipeline/files/bcfstats.Robject",sep=""))
  print("done")
#  }

  ## plotting per individual across chroms: TSTV
  
  load(file=paste("~/Great_Ape/genome_pipeline/files/bcfstats.Robject",sep=""))
  load(file=paste("~/Great_Ape/genome_pipeline/files/covstats.Robject",sep=""))
  
  chco<-list()
  ai<-list();va=0
  for (spe in (1:length(spec))) { 
    nave<-rep(unlist(namlist[[spe]]),23)
    pm<-ceiling(length(namlist[[spe]])/20)
    ain<-1:length(namlist[[spe]])
    tsv<-list()
    for (chr in (1:23)) { tsv[[chr]]<-sustat[[chr]][[spe]][5,] }
    chsub<-as.numeric(unlist(do.call(cbind,tsv)))
    for (j in (1:pm)) {
      va=va+1
      suse<-namlist[[spe]][which(ain>(j*20)-20 & ain<=(j*20))] 
      cr=data.frame(typ=rep(suse,23),tes=chsub[which(nave%in%suse)])
      cr<-cr[!is.na(cr[,2]),]
      if (length(suse)<20) { lvs<-unique(cr$typ); nr<-data.frame(typ=mpty[1:(20-length(suse))],tes=rep(NA,(20-length(suse))));cr<-rbind(cr,nr);cr$typ<-factor(cr$typ, levels=c(lvs,mpty[1:(20-length(suse))])) }
      p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) + geom_segment(aes(x=0.0,y=min(tes)*0.9,xend=0.0,yend=max(tes)*1.1)) +xlab("")+ylab("")+ ggtitle(label=paste("ts/tv ratio ",spec[[spe]],", part ",j,"/",pm,sep=""))
      ai[[va]]<-p1
    }
  }  
  
  ## plot ts/tv
  pdf(paste(basedir,"tstv_per_individual.pdf",sep=""),width=24,height=57)
  grid.arrange(ai[[1]],ai[[2]],ai[[3]],ai[[4]],ai[[5]],ai[[6]],ai[[7]],ai[[8]],ai[[9]],ai[[10]],ai[[11]],ai[[12]],ai[[13]],ai[[14]],ai[[15]],ai[[16]],ai[[17]],ai[[18]],ncol=1) 
  dev.off()
  
  print("done")

  ## plotting across inds per chrom: last pos
  labl="last pos"
  ap<-list()
  for (chr in c(1:24)) { 
    crl<-list()
    for (spe in (1:length(spec))) { crl[[spe]]<-sustat[[chr]][[spe]][1,] }
    cr=data.frame(typ=typC,tes=as.numeric(unlist(crl)))
    if(chr==24) { cr<-cr[which(xv=="M"),] }
    cr<-cr[!is.na(cr[,2]),]
    p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1.3) ,plot.margin=margin(0.1,20,0.1,10),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) + geom_segment(aes(x=0.0,y=min(tes)*0.999,xend=0.0,yend=max(tes)*1.001))  +xlab("")+ylab("")    + ggtitle(label=paste(labl," chr",ifelse(chr==23,"X",ifelse(chr==24,"Y",chr)),sep=""))
    ap[[chr]]<-p1
  }  
  
  # plot last position for chromosome per individual
  pdf(paste(basedir,"last_pos_per_chrom.pdf",sep=""),width=24,height=52)
  grid.arrange(ap[[1]],ap[[2]],ap[[3]],ap[[4]],ap[[5]],ap[[6]],ap[[7]],ap[[8]],ap[[9]],ap[[10]],ap[[11]],ap[[12]],ap[[13]],ap[[14]],ap[[15]],ap[[16]],ap[[17]],ap[[18]],ap[[19]],ap[[20]],ap[[21]],ap[[22]],ap[[23]],ap[[24]],nrow=12,ncol=2)
  dev.off()
  

  ## plotting across inds per chrom: SNPs (SiS)
  labl="heterozygous sites"
  ap<-list()
  for (chr in c(1:23)) { 
    crl<-list()
    for (spe in (1:length(spec))) { crl[[spe]]<-sustat[[chr]][[spe]][4,] }
    cr=data.frame(typ=typC,tes=as.numeric(unlist(crl)))
    if(chr==23) { cr<-cr[which(xv=="F"),] }
    cr<-cr[!is.na(cr[,2]),]
    p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=0.95) ,plot.margin=margin(0.1,20,0.1,5),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) +xlab("")+ylab("")+ geom_segment(aes(x=0.0,y=min(tes)*0.9,xend=0.0,yend=max(tes)*1.1))  + ggtitle(label=paste(labl," chr",ifelse(chr==23,"X",ifelse(chr==24,"Y",chr)),sep=""))
    ap[[chr]]<-p1
  }  
  
  # plot number of hets
  pdf(paste(basedir,"hets_per_chrom.pdf",sep=""),width=24,height=52)
  grid.arrange(ap[[1]],ap[[2]],ap[[3]],ap[[4]],ap[[5]],ap[[6]],ap[[7]],ap[[8]],ap[[9]],ap[[10]],ap[[11]],ap[[12]],ap[[13]],ap[[14]],ap[[15]],ap[[16]],ap[[17]],ap[[18]],ap[[19]],ap[[20]],ap[[21]],ap[[22]],ap[[23]],nrow=12,ncol=2)
  dev.off()

  ## plotting across inds per chrom: Non-Reference positions (SiS)
  labl="Non-Reference sites"
  ap<-list()
  for (chr in c(1:24)) { 
    crl<-list()
    for (spe in (1:length(spec))) { crl[[spe]]<-sustat[[chr]][[spe]][3,] }
    cr=data.frame(typ=typC,tes=as.numeric(unlist(crl)))
    if(chr==24) { cr<-cr[which(xv=="M"),] }
    cr<-cr[!is.na(cr[,2]),]
    p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=0.95) ,plot.margin=margin(0.1,20,0.1,5),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) + geom_segment(aes(x=0.0,y=min(tes)*0.9,xend=0.0,yend=max(tes)*1.1))    +xlab("")+ylab("")+ ggtitle(label=paste(labl," chr",ifelse(chr==23,"X",ifelse(chr==24,"Y",chr)),sep=""))
    ap[[chr]]<-p1
  }  
  
  # plot number of SNPs
  pdf(paste(basedir,"norefs_per_chrom.pdf",sep=""),width=24,height=52)
  grid.arrange(ap[[1]],ap[[2]],ap[[3]],ap[[4]],ap[[5]],ap[[6]],ap[[7]],ap[[8]],ap[[9]],ap[[10]],ap[[11]],ap[[12]],ap[[13]],ap[[14]],ap[[15]],ap[[16]],ap[[17]],ap[[18]],ap[[19]],ap[[20]],ap[[21]],ap[[22]],ap[[23]],ap[[24]],nrow=12,ncol=2)
  dev.off()
  
  ## plotting across inds per chrom: records
  labl="Records"
  ap<-list()
  for (chr in c(1:24)) { 
    crl<-list()
    for (spe in (1:length(spec))) { crl[[spe]]<-sustat[[chr]][[spe]][2,] }
    cr=data.frame(typ=typC,tes=as.numeric(unlist(crl)))
    if(chr==24) { cr<-cr[which(xv=="M"),] }
    cr<-cr[!is.na(cr[,2]),]
    p1<-ggplot(cr,aes(factor(typ), tes)) + theme_bw() + geom_violin(scale="width",aes(fill=factor(typ)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=0.9) ,plot.margin=margin(0.1,10,0.1,5),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15)) + geom_segment(aes(x=0.0,y=min(tes)*0.9,xend=0.0,yend=max(tes)*1.1))  +xlab("")+ylab("")+ ggtitle(label=paste(labl," chr",ifelse(chr==23,"X",ifelse(chr==24,"Y",chr)),sep=""))
    ap[[chr]]<-p1
  }  
  
  pdf(paste(basedir,"records_per_chrom.pdf",sep=""),width=24,height=52)
  grid.arrange(ap[[1]],ap[[2]],ap[[3]],ap[[4]],ap[[5]],ap[[6]],ap[[7]],ap[[8]],ap[[9]],ap[[10]],ap[[11]],ap[[12]],ap[[13]],ap[[14]],ap[[15]],ap[[16]],ap[[17]],ap[[18]],ap[[19]],ap[[20]],ap[[21]],ap[[22]],ap[[23]],ap[[24]],nrow=12,ncol=2)
  dev.off()
  
#  q()
#  }


  ################################################################################
  ## coverage distribution across individuals
  
  load(file="~/Great_Ape/genome_pipeline/files/covdists.Robject")
  spcol=c("#0b9c2efe","#4aa761fe","#4c7b58fe","#0b5b1efe",
          "#9b001fff",
          "#cd022bff","#e54263ff","#cd6278ff","#d00905ff",
          "#b16800ff","#d59e4eff","#c6a06aff",
          "grey","grey","grey")
    
  ax<-1:501
  pdf(paste(basedir,"coverage_distributions.pdf",sep=""),width=24,height=52)
  par(mfcol=c(15,1))
  for (j in (1:length(spec))) {
    allm<-max(unlist(covdist[[j]]))*1.1
    k=1
    ay=covdist[[j]][[k]]
    par(new=F)
    cl=spcol[j]
    plot(y=ay,x=ax,type="l",xlim=c(0,ifelse(j==13,200,ifelse(j%in%c(14,15,7),90,60))),ylim=c(0,allm),xlab="Coverage",ylab="Count",bty="n",col=alpha(cl,.99),main=paste("Coverage ",spec[j],sep=""),cex.main=2,cex.lab=1.2,cex.axis=1.2)
    if(spec[j]=="Pongo_tapanuliensis") { next }
    for (k in (2:length(namlist[[j]]))) {
      ay=covdist[[j]][[k]]
      par(new=T)
      plot(y=ay,x=ax,type="l",xlim=c(0,ifelse(j==13,200,ifelse(j%in%c(14,15,7),90,60))),ylim=c(0,allm),xlab="",ylab="",bty="n",col=alpha(cl,.99),axes=F,main="")
    }
  }
  dev.off()
  
