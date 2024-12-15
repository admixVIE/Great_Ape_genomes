#!/usr/bin/r

## in preparation done once:
## wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToRheMac10.over.chain.gz; gunzip hg38ToRheMac10.over.chain.gz

method="getanc"
# getanc = create ancestral allele bed file 

options(scipen=100)
input=as.character(unlist((commandArgs(TRUE))))
print(input)

if (method=="getanc") {

  basedir=""
  library(rtracklayer)
  input<-unlist(strsplit(as.character(input),split=" "))
  chrom=as.numeric(as.character(input[1]))
  chrom<-ifelse(chrom==23,"X",chrom)
  spec=as.character(input[2])
  dir=as.character(input[3])
  print(c(spec,chrom))
  
  cm<-c("A","C","G","T","N",".","-")
  names(cm)<-c("T","G","C","A","N",".","-")
  cmfun<-function(ip) { as.character(cm[which(names(cm)==as.character(ip))]) }
  
  hgchain<-import.chain(paste(basedir,"hg",ifelse(spec=="human",19,38),"ToRheMac10.over.chain",sep=""))
  chainfunc<-function(intab) { liftOver(GRanges(seqnames=intab[,1],ranges=IRanges(start=as.numeric(intab[,2]),end=as.numeric(intab[,2]),width=rep(1,nrow(intab)),names=intab[,1]),strand=rep("+",nrow(intab)),mcols=intab[,2]),hgchain) }
  
  yy=0;yle=seq(0,250000000,1000000)
  repeat {
    yy=yy+1
    print(yy)
    if (yy==length(yle)) { break}
    intvl<-c(yle[yy],yle[yy+1])
    # 1) get segregating sites from hg38 coordinates into rhemac10 coordinates and write as bed
    vcf<-system(paste("ls ",dir,"/",spec,"/",ifelse(spec=="human","ALL.chr","chr"),chrom,ifelse(spec=="human",".*.vcf.gz",".vcf.gz"),sep=""),intern=T)
    vcftabb<-system(paste("Software/bcftools/bcftools-1.16/bcftools query -r ",ifelse(spec=="human","","chr"),chrom,":",intvl[1],"-",intvl[2]," -f '%POS\n' ",vcf,sep=""),intern=T)
    if (length(vcftabb)==0) { next }
    ttab<-cbind(paste("chr",chrom,sep=""),vcftabb)
    ttabch<-chainfunc(ttab)
    ttabch<-as.data.frame(ttabch)[,c(3:5,7,8),drop=F]
    if (nrow(ttabch)==0) { next }
    ttabch[,2]<-ttabch[,2]-1
    ttabch<-ttabch[order(ttabch[,2]),,drop=F]
    ttabch<-ttabch[order(ttabch[,3]),,drop=F]
    write.table(cbind(ttabch[,c(1:3),drop=F],NA,NA,ttabch[,4,drop=F]),paste(basedir,"tmp/tmp_",spec,chrom,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
    
    # 2) get rheMac genotypes
    ovrlift<-system(paste(basedir,"bedtools2/bin/./bedtools getfasta -tab -s -fi ",basedir,"rheMac10.fa -bed ",basedir,"/tmp/tmp_",spec,chrom,".bed -fo stdout",sep=""),intern=T)
    ovrlift<-do.call(rbind,strsplit(ovrlift,split="\t"))
    ovrlift<-cbind(do.call(rbind,strsplit(ovrlift[,1,drop=F],split=":")),ovrlift[,2])
    ovrlift<-cbind(ovrlift[,1,drop=F],do.call(rbind,strsplit(ovrlift[,2,drop=F],split="-"))[,1,drop=F],toupper(ovrlift[,3,drop=F]))
    system(paste("rm ",basedir,"/tmp/tmp_",spec,chrom,".bed",sep=""))
    
    # 3) turn rheMac genotypes into hg38/hg19 coordinates
    ttabhg<-cbind(paste("chr",chrom,sep=""),as.numeric(ttabch[,5])-1,ttabch[,5],ovrlift[,3])
    ttabhg<-ttabhg[order(as.numeric(ttabhg[,2])),,drop=F]
    ttabhg<-ttabhg[order(as.numeric(ttabhg[,3])),,drop=F]
    
    write.table(unique(ttabhg),paste(basedir,"ancgenome/anc_",spec,"_",chrom,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F,append=T)
  }
  
  print(c("done",spec,chrom))
  print(Sys.time())
  q()
  }  

