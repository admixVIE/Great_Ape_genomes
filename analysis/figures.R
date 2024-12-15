#!/usr/bin/r

'%ni%' <- Negate('%in%')
options("scipen"=100)

################################################################################
## Figure 2: plot of coverage & heterozygosity for species/subspecies (local)

load(file="genome_pipeline/files/bcfstats.Robject")
load(file="genome_pipeline/files/covstats.Robject")

require(ggplot2)
require(gridExtra)
library("dplyr")

cpt<-ifelse(grepl("captive",typC)==T,"CP","WB")
typ3[which(typ3=="Gorilla gorilla")]<- "Gorilla gorilla gorilla"
typ3[which(typ3=="Pan troglodytes")]<- "Pan troglodytes verus"

spec<-c("Gorilla_beringei_beringei", "Gorilla_beringei_graueri","Gorilla_gorilla_diehli","Gorilla_gorilla_gorilla","Pan_paniscus","Pan_troglodytes_ellioti","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pongo_abelii","Pongo_pygmaeus","Pongo_tapanuliensis","captive_Pan","captive_Gorilla","captive_Pongo")
aspecs<-c("Gorilla beringei beringei","Gorilla beringei graueri", "Gorilla gorilla diehli", "Gorilla gorilla gorilla",
          "Pan paniscus",
          "Pan troglodytes ellioti","Pan troglodytes schweinfurthii", "Pan troglodytes troglodytes"   ,"Pan troglodytes verus" ,
          "Pongo abelii" , "Pongo pygmaeus","Pongo tapanuliensis" ,
          "captive Pan","captive Gorilla","captive Pongo"
)
spcol=c("#0b9c2efe","#4aa761fe","#4c7b58fe","#0b5b1efe",
        "#9b001fff",
        "#cd022bff","#e54263ff","#cd6278ff","#d00905ff",
        "#b16800ff","#d59e4eff","#c6a06aff",
        "grey","grey","grey")
names(spcol)<-aspecs
scol <-  c("black","blue");names(scol)=c("WB","CP")
scls<-ifelse(cpt!="CP",scol[1],scol[2])


## coverage
cr=data.frame(typ=typ3,tes=as.numeric(do.call(rbind,covlist)[,2]),cpt=cpt)
sclsM<-scls; sclsM[which(cr$tes>70)]<-"dark blue"
cr$tes[which(cr$tes>70)]<-70
p1<-  ggplot(cr) + theme_bw() + geom_violin(mapping=aes(x=typ,y=tes,fill=typ), scale="width",adjust=1.0,draw_quantiles = c(0.5) , na.rm=T  ) + geom_jitter(mapping=aes(x=typ,y=tes,fill=cpt), height = 0, width = 0.2,na.rm=T,inherit.aes=F,show.legend=F,colour=sclsM) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=13),axis.text.x= element_text(size=14,angle=45,hjust=1,vjust=1.1),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15),plot.margin=margin(1,0.1,0.1,25) )  +xlab("")+ylab("") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=70)) +scale_fill_manual(values=spcol)    + ggtitle(label="coverage (fold)")

## (approximate) heterozygosity
cl<-list();for (spe in (1:length(spec))) { cl[[spe]]<-list() }
for (chr in c(1:23)) {  for (spe in (1:length(spec))) { cl[[spe]][[chr]]<-sustat[[chr]][[spe]][4,] } }
for (spe in (1:length(spec))) { cl[[spe]]<-colSums(matrix(as.numeric(do.call(rbind,cl[[spe]]))/3000000000*1000,nrow=length(cl[[spe]]))) }
cr2<-data.frame(typ=typ3,tes=unlist(cl),cpt=cpt)
scls2<-scls[!is.na(cr2[,2])]
cr2<-cr2[!is.na(cr2[,2]),]


p2<-ggplot(cr2) + theme_bw() + geom_violin(mapping=aes(x=typ,y=tes,fill=typ), scale="width",adjust=1.0,draw_quantiles = c(0.5) , na.rm=T  ) + geom_jitter(mapping=aes(x=typ,y=tes,fill=cpt), height = 0, width = 0.2,na.rm=T,inherit.aes=F,show.legend=F,colour=scls2) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=13),axis.text.x= element_text(size=14,angle=45,hjust=1,vjust=1.2),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15),plot.margin=margin(1,0.1,0.1,35) )  +xlab("")+ylab("") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=4)) +scale_fill_manual(values=spcol)    + ggtitle(label="heterozygosity (per 1,000bp)")

p1a <- p1 + labs( tag = "A") 
p2a <- p2 + labs( tag = "B") 

pdf(paste("genome_pipeline/files/Figure2.pdf",sep=""),7,9)
grid.arrange(p1a, p2a,ncol = 1)         
dev.off()


################################################################################
## SI figures on captive panel: f3-stats, rareCAGA for all

############### f3-stats
library(ggplot2)
library(tidyverse)
library(ggh4x)

sps<-c("pan_all","gorilla_all","pongo_all");names(sps)<-c("Pan","Gorilla","Pongo")
pops<-list(c("PTV","PTT","PTE","PTS","PPA"),c("GBB","GBG","GGD","GGG"),c("PT","PA","PP")); names(pops)<-sps
load(file="genome_pipeline/files/covstats.Robject")
cpt<-unlist(read.table(file=paste("genome_pipeline/files/captive.txt",sep=""),sep="\t",header=F,as.is=T))

## load all f3 results
fastats<-list();lablist<-c()
for (j in (1:length(cpt))) {
  samp=cpt[j]
  spec<-cvl[which(cvl[,3]==samp),2]
  fstats<-try(read.table(file=paste("analysis/f3_",samp,".tsv",sep=""),sep="\t",header=T,colClasses=c("character","character","character", "numeric","numeric","numeric","numeric")))
  spe=which(fstats$est==max(fstats$est))[1]
  colr<-rep("black",nrow(fstats));colr[spe]<-"darkred";fstats$colr<-colr
  fstats$labl=paste(samp,", ",gsub("CP_","",spec),", ", fstats[spe,]$pop2,sep="")
  fstats$spe=fstats[spe,]$pop2
  fastats[[j]]<-fstats
  lablist[j]<-fstats[spe,]$pop2
}

sustat<-do.call(rbind,fastats)
cps<-data.frame(cps=cpt,lablist=lablist) %>% arrange(across(lablist, ~factor(., levels=unlist(pops))))

sta<-tibble(Subspecies=sustat$pop2,f3=sustat$est,mn=sustat$est-sustat$se,spe=sustat$spe,mx=sustat$est+sustat$se,Sample=sustat$pop3,labl=sustat$labl,colr=sustat$colr) %>% mutate(across(Sample, ~factor(., levels=cps$cps)))
lb =  sta$labl; names(lb)=sta$Sample

## sort by subspecies
  
splots<-
    ggplot(sta, aes(x=Subspecies,y=f3,mn=mn,mx=mx))  +  geom_point(size=3,colour=sta$colr) + geom_errorbar(aes(ymin = mn,ymax = mx),linewidth=1,colour=sta$colr) +
    theme_minimal() + theme(panel.grid.major = element_blank(),strip.background= element_rect(fill = "transparent", colour = NA), panel.grid.minor = element_blank(),plot.title = element_text(face="bold",hjust=0.5,size=15),axis.text.x = element_text(angle = 45, vjust = 1.3, hjust=1),axis.line.y = element_line(color="black"),axis.ticks.y=element_line(color="black"),strip.text.x = element_text(size=12, color="darkblue")) +ylab("f3(Macaque;Ape,Test)") +xlab("Subspecies") +
    facet_wrap2(~Sample,scale="free", labeller = labeller(Sample = lb ), nrow=10,ncol=7, trim_blank = FALSE)
#+ ggtitle(label=labl) 
pdf(paste("genome_pipeline/files/f3_plots.pdf",sep=""),18,28)
print(splots) 
dev.off()

## table for SI
sra<-read.table("genome_pipeline/files/captive_sra.txt",sep="\t",header=T,as.is=T,fill=T)
f3<-read.table("genome_pipeline/files/f3_stats_captive.txt",sep="\t",header=T,as.is=T,fill=T)
load("genome_pipeline/files/hucontest.Robject")
huc<-cbind(names(huco),do.call(rbind,strsplit(as.character(huco),split=" "))[,1])
sfin<-merge(sra,f3,by.x=1,by.y=1,all=T)
sfin<-merge(sfin,huc,by.x=1,by.y=1,all=T)
sfin<-sfin[order(sfin[,3]),c(2,1,8,3,4,5)]
colnames(sfin)<-c("SRA_sample","Name","Human_contamination_(%)","SRA_species","SRA_short","f3_short")

write.table(sfin,file="genome_pipeline/files/captive_finalspec.txt",sep="\t",col.names=T,row.names=F,quote=F)

################################################################################
## Figure 5: plot of captive panel results

## A: f3-stats; B: HuConTest vioplot; C: rareCAGA
#module load R/4.2.3 bcftools conda; conda activate geo ; conda activate gdal;conda activate proj; LD_LIBRARY_PATH=$LD_LIBRARY_PATH: R --vanilla

## Libraries
library(gstat)
library(sp)
library(raster)
library(spData)
library(dplyr)
library(sf)
library(tmap)
library(tmaptools)
library(maptools)
library(automap)
library(spatialEco)
library(splancs)
library(ggplot2)
library(rgdal)
library(mapview)
library(TSCS)
library(tidyverse)
library(ggh4x)
library(gridExtra)
library(grid)

## preparation
sps<-c("pan_all","gorilla_all","pongo_all");names(sps)<-c("Pan","Gorilla","Pongo")
load(file="genome_pipeline/files/covstats.Robject")
cpt<-unlist(read.table(file=paste("genome_pipeline/files/captive.txt",sep=""),sep="\t",header=F,as.is=T))
'%ni%' <- Negate('%in%')

############
## HuConTest
load("genome_pipeline/files/hucontest.Robject")
## violin plots
aco<-data.frame(spec=cvl[grep("captive",cvl[,1]),1],cont=as.numeric(do.call(rbind,strsplit(as.character(huco),split=" "))[,1]))
aco$spec<-gsub("captive_","",aco$spec)
panel_a<-ggplot(aes(factor(spec), cont),data=aco) + theme_bw() + geom_violin(scale="width",aes(fill=factor(spec)),adjust=1.0,draw_quantiles = c(0.5)  ) + geom_jitter(height = 0, width = 0.1) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.title.y = element_text(size = 13,hjust=-5),axis.text.y= element_text(size=12),axis.text.x= element_text(size=12,angle=45,hjust=1,vjust=1) , axis.ticks.x = element_blank() ) +xlab("")+ylab("Human contamination (%)") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=1.5)) 

############
## f3-stats
samp<-cpt[which(cpt=="PD_0179")]
spec<-cvl[which(cvl[,3]==samp),2]
fstats<-try(read.table(file=paste("greatapes/analysis/f3_",samp,".tsv",sep=""),sep="\t",header=T,colClasses=c("character","character","character", "numeric","numeric","numeric","numeric")))
spe=which(fstats$est==max(fstats$est))[1]
colr<-rep("black",nrow(fstats));colr[spe]<-"darkred";fstats$colr<-colr
fstats$labl=paste(samp,", ",gsub("CP_","",spec),", ", fstats[spe,]$pop2,sep="")
fstats$spe=fstats[spe,]$pop2

sta<-tibble(Subspecies=fstats$pop2,f3=fstats$est,mn=fstats$est-fstats$se,spe=fstats$spe,mx=fstats$est+fstats$se,Sample=fstats$pop3,labl=fstats$labl,colr=fstats$colr)
panel_b<-  ggplot(sta, aes(x=Subspecies,y=f3,mn=mn,mx=mx))  +  geom_point(size=3,colour=sta$colr) + geom_errorbar(aes(ymin = mn,ymax = mx),linewidth=1,colour=sta$colr) +
  theme_minimal() + theme(panel.grid.major = element_blank(),strip.background= element_rect(fill = "transparent", colour = NA), panel.grid.minor = element_blank(),axis.text.y = element_text(size=12),axis.text.x = element_text(size=12,angle = 45, vjust = 1.3, hjust=1),axis.title.y = element_text(size = 13,hjust=5),axis.title.x = element_text(size = 13),axis.ticks.y=element_line(color="black")) +ylab("f3(Macaque;Test,PD_0179)") +xlab("Subspecies") + geom_segment(aes(x=0.0,y=0.44,xend=0.0,yend=0.52)) 


############
## rareCAGA
load(paste("rareCAGA/results/maptab_captive",sep=""))
load(paste("rareCAGA/results/R1_captive",sep=""))

ext<-t(matrix(c(-17,35,-12,17),2,2))
colnames(ext)<-c("min","max")
rownames(ext)<-c("coords.x1","coords.x2")
crs    <- CRS("+proj=longlat +ellps=WGS84") 

tn<-which(names(allTs)=="PD_0259")
panel_c<-allTs[[tn]]+tm_layout(title="PD_0259; Pan troglodytes") 

t2<-ggplot(x=0, y = 0) + annotate(geom="label", x = c(0.003,0.485), y = c(0.95,0.95),label = c("A","B"),colour = "black", size=7,fill="white",label.size=NA) + theme_void()  +scale_y_continuous(limits = c(0, 1)) +scale_x_continuous(limits = c(0, 1))
t3<-ggplot(x=0, y = 0) + annotate(geom="label", x = 0.01, y = 0.958,label = "C",colour = "black", size=7,fill="white",label.size=NA) + theme_void()  +scale_y_continuous(limits = c(0, 1)) +scale_x_continuous(limits = c(0, 1))


p_a <- panel_a 
p_b <- panel_b  
p_c <- panel_c 
lay <- grid.layout(nrow=2,ncol=1,heights=c(0.5,0.5))

pdf(paste("genome_pipeline/files/Figure5.pdf",sep=""),8,6)
pushViewport(viewport(layout=lay))
pab<-grid.arrange(arrangeGrob(p_a, p_b,ncol = 2,nrow=2,heights=c(0.5,0.5)),newpage=F)
print(pab,vp=viewport(layout.pos.col = 1,layout.pos.row = 1))
print(p_c, vp=viewport(layout.pos.col = 1,layout.pos.row = 2))
print(t2, vp=viewport(layout.pos.col = 1,layout.pos.row = 1))
print(t3, vp=viewport(layout.pos.col = 1,layout.pos.row = 2))

dev.off()


########################
## rareCAGA in one graph (but reordered)
nworder<-c("Ai","Akira","Ayumu","AG18354","Frits","Carolina","Carl","Simliki","A","B","C","D","E","F","G","H","I","Donald","CH114","CH391","PD_0257","PD_0258","PD_0259","PD_0427")
allT2<-list()
for (j in (1:length(nworder))) { allT2[[j]]<-allTs[[which(names(allTs)==nworder[j])]] }
allT2[[which(nworder=="C")]]<-allT2[[which(nworder=="C")]]+tm_layout(title="C") 

pdf(paste("museum_apes/rareCAGA/resultsgeoloc_captive_chimps.pdf",sep=""),12,24)
print(tmap_arrange(allT2,ncol=2,nrow=12)) 
dev.off()


################################################################################
## Relatedness & inbreeding results
#module load --auto bcftools R/4.2.3; LD_LIBRARY_PATH=$LD_LIBRARY_PATH: R --vanilla

library(tidyverse)
cvl<-read.table(file=paste("~/Great_Ape/genome_pipeline/files/metadata.txt",sep=""),header=T,as.is=T,sep="\t")

for (spec in c("pan","gorilla","pongo")) {
  tabl=read.table(paste("plink/ngsr/",spec,"_all.res",sep=""),sep="\t",header=T,as.is=T)
  nams<-system(paste("bcftools query -l ",spec,"_all/allchr.vcf.gz",sep=""),intern=T)
  names(nams)<-c(0:(length(nams)-1))
  specs<-cvl[match(nams,cvl$Name),c(2,3)]
  specs<-specs[order(specs$Species_label),]
  
  # make a heatmap
  tabl$a1<-nams[match(tabl$a,names(nams))]
  tabl$b1<-nams[match(tabl$b,names(nams))]
  tabl$KING<-ifelse(tabl$KING<0,0,tabl$KING)
  tbl2 <- data.frame(ida=tabl$a1, idb=tabl$b1, Relatedness=tabl$KING)
  tbl2<- tbl2 %>% mutate(across(ida, ~factor(.x, specs$Name))) %>%  mutate(across(idb, ~factor(.x, specs$Name))) 
  p1<-ggplot(tbl2,aes(ida, idb)) +
    geom_tile(aes(fill =  Relatedness), colour = "gray") +
    scale_fill_gradient("Relatedness", low = "grey90", high = "dark red") +
    scale_alpha_identity() + xlab("Individual A") + ylab("Individual B") + ggtitle(paste("Relatedness among ",spec,sep="")) +
    theme(plot.title = element_text(face="bold",hjust=0.5,size=13),axis.text.x= element_text(angle=90,hjust=1,vjust=.5)) +
    coord_equal(expand = FALSE)
  
  pdf(paste("analysis/relatedness_",spec,".pdf",sep=""),nrow(specs)*0.2,nrow(specs)*0.2)
  print(p1)
  dev.off()
  print(c(spec,nrow(specs),nrow(tbl2)))
  
}



########################################################
## plot of ROHs

sps<-c("pan","gorilla","pongo")
pops<-list(c("PTT","PTS","PTE","PTV","PPA","CP_Pan"),c("GGG","GGD","GBG","GBB","CP_Gorilla"),c("PA","PP","PT","CP_Pongo")); names(pops)<-sps
spcolx=c("#0b9c2efe","#4aa761fe","#4c7b58fe","#0b5b1efe",
        "#9b001fff","grey",
        "#cd022bff","#e54263ff","#cd6278ff","#d00905ff","grey",
        "#b16800ff","#d59e4eff","#c6a06aff",
        "grey")
names(spcolx)<-unlist(pops)

load(file="Great_Ape/genome_pipeline/files/covstats.Robject")
rtab<-list()
for (spe in sps) {
  rtab[[spe]]<-list()
  for (chr in c(1:22)) { 
    rtab[[spe]][[chr]]<-read.table(paste("analysis/",spe,".",chr,".rohs.txt",sep=""),sep="\t",header=F,as.is=T,fill=T)
    }
  }
  
allrohs<-list()
alens<-c()
alons<-c()
for (ind in (cvl[,3])) { 
  spec=cvl[which(cvl[,3]==ind),2]
  print(ind)
  spe=ifelse(spec%in%pops[[1]],"pan",ifelse(spec%in%pops[[2]],"gorilla",ifelse(spec%in%pops[[3]],"pongo",NA)))
  aroh=list()
  for (chr in c(1:22)) { 
    rohf<-rtab[[spe]][[chr]]
    aroh[[chr]]<-rohf[which(rohf[,2]==ind),6]
  }
  aroh<-unlist(aroh)
  alens[ind]<-sum(aroh)  
  alons[ind]<-sum(aroh[which(aroh>500000)])
  allrohs[[ind]]<-aroh
}

aplot<-data.frame(Species=cvl[,2],rohs=alens/1000000,rrohs=alons/1000000)
aplot$Species<-factor(aplot$Species,unlist(pops))

p1<-  
  ggplot(aplot) + theme_bw() + geom_violin(mapping=aes(x=Species,y=rrohs,fill=Species), scale="width",adjust=1.0,draw_quantiles = c(0.5) , na.rm=T  ) + geom_jitter(mapping=aes(x=Species,y=rrohs), height = 0, width = 0.2,na.rm=T,inherit.aes=F,show.legend=F) + theme(panel.border = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"),legend.position="none",axis.text.y= element_text(size=13),axis.text.x= element_text(size=14,angle=45,hjust=1,vjust=1.1),axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=15),plot.margin=margin(1,0.1,0.1,25),axis.title.y = element_text(size=14) )  +xlab("")+ylab("Mbp in RoHs > 500kbp") + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=400)) +scale_fill_manual(values=spcolx)    + ggtitle(label="Runs of homozygosity per population")

pdf(paste("genome_pipeline/files/Figure4.pdf",sep=""),7,9)
p1
dev.off()




################################################################################
## Tables SI



