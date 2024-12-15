##################################################################################
#### basic analyses for captive panel: f-statistics, hucontest, rarecaga

method=unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" "))[3]
print(method)

##### METHODS #######
# famcor = correct fam files for f-stats
# fstatc = precompute f-statistics for a given individual
# fstatbest = best fit f-statistics all individuals
# hucon = human contamination test

## method="hucon"

'%ni%' <- Negate('%in%')
options("scipen"=100)


################################################################################
### f3-statistics

if (method == "famcor") {
  ### first (and once), modify the fam files!
  sps<-c("pan_all","gorilla_all","pongo_all")
  load(file="Great_Ape/genome_pipeline/files/covstats.Robject")
  cvl[which(cvl[,3]=="C"),3]<-"C_a";  cvl[which(cvl[,3]=="Suma"),3]<-"Sum"
  for (spec in (sps)) { 
   snam=ifelse(spec=="pan_all","Pan",ifelse(spec=="gorilla_all","Gorilla",ifelse(spec=="pongo_all","Pongo",NA)))
   ofam<-read.table(paste("plink/dataset.",spec,".fam",sep=""),sep="\t",header=F,as.is=T)
   omet<-cvl[grep(snam,cvl[,1]),]
   ofo1<-match(ofam[,2],omet[,3])
   ofam2<-cbind(omet[ofo1,c(2,3)],ofam[,-c(1:2)])
   ofam2[which(ofam[,2]=="MAC"),]<-c("MAC",ofam[which(ofam[,2]=="MAC"),-1])
   ofam2[grep("CP",ofam2[,1]),1]<-ofam2[grep("CP",ofam2[,1]),2]
   system(paste("mv ","plink/dataset.",spec,".fam ", "plink/old_dataset.",spec,".fam",sep=""))
   write.table(ofam2,paste("plink/dataset.",spec,".fam",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
   cps<-cvl[grep("CP",cvl[,2]),3]
   write.table(cps,file=paste("Great_Ape/genome_pipeline/files/captive.txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
   }
  }

if (method == "fstatc") {
  library(admixtools)
  ip=unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" "))
  print(ip)
  samp=ip[1]
  dir=ip[2]
  print(samp)
  sps<-c("pan_all","gorilla_all","pongo_all");names(sps)<-c("CP_Pan","CP_Gorilla","CP_Pongo")
  load(file="Great_Ape/genome_pipeline/files/covstats.Robject")
  cvl[which(cvl[,3]=="C"),3]<-"C_a";  cvl[which(cvl[,3]=="Suma"),3]<-"Sum"

  # prepare subset metadata
  spec<-sps[cvl[which(cvl[,3]==samp),2]]
  pops<-list(c("PTT","PTV","PPA","PTE","PTS"),c("GBB","GBG","GGD","GGG"),c("PT","PA","PP")); names(pops)<-sps
  pops=pops[[spec]]
  ofam<-read.table(paste("plink/dataset.",spec,".fam",sep=""),sep="\t",header=F,as.is=T)
  
  # precalculate f2-stats
  loc<-paste(dir,"/",samp,"/",sep="")
  extract_f2(pref=paste("plink/dataset.",spec,sep=""),pops=c("MAC",pops,samp),outdir=loc,blgsize=500000)
  f2_blocks = f2_from_precomp(loc)
  ff<-f3(f2_blocks,pop1="MAC",pop2=pops,pop3=samp)
  ## save f3-stats per sample into plain text table
  write.table(ff,file=paste("analysis/f3_",samp,".tsv",sep=""),sep="\t",col.names=T,row.names=F,quote=F)

  print("done")
  print(Sys.time())

  q()
}


########################################################################
######## get the best fit in f3-stats

if (method == "fstatbest") {
  library(admixtools)
  sps<-c("pan_all","gorilla_all","pongo_all");names(sps)<-c("CP_Pan","CP_Gorilla","CP_Pongo")
  pops<-list(c("PTT","PTV","PPA","PTE","PTS"),c("GBB","GBG","GGD","GGG"),c("PT","PA","PP")); names(pops)<-sps
  load(file="Great_Ape/genome_pipeline/files/covstats.Robject")
  cps<-unlist(read.table(file=paste("Great_Ape/genome_pipeline/files/captive.txt",sep=""),sep="\t",header=F,as.is=T))

  bestfit<-list()
  for (j in (1:length(cps))) {
    samp=cps[j]
    spec<-sps[cvl[which(cvl[,3]==ifelse(samp=="C_a","C",samp)),2]]
    pop=pops[[spec]]
    fstats<-try(read.table(file=paste("analysis/f3_",samp,".tsv",sep=""),sep="\t",header=T,as.is=T))
    if(ncol(fstats)==1) {   bestfit[[j]]<-c(samp,NA, NA,NA); next }
    bestfit[[j]]<-c(ifelse(samp=="C_a","C",samp),unlist(fstats[which(fstats[,4]==max(fstats[,4])),c(2,4:5)][1,]))
  }
  
  bestfit<-do.call(rbind,bestfit)
  colnames(bestfit)<-c("Individual","Bestfit","f3","se")
  write.table(bestfit,file=paste("Great_Ape/genome_pipeline/files/f3_stats_captive.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  
  print("done")
  print(Sys.time())
  
  q()
  }


########################################################
######## human contamination test 

if (method == "hucon") {
  gnom=unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" "))[2]
  seqpath="Result/Sequences/"
  load(file="Great_Ape/genome_pipeline/files/covstats.Robject")
  ind=unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" "))[1]
#  huco<-list()
  bcft="Software/bcftools/bcftools-1.16/bcftools"
  setwd(paste("HuConTest/",sep=""))
  sps<-cbind(c("CP_Pan","CP_Gorilla","CP_Pongo"),c("captive_Pan","captive_Gorilla","captive_Pongo"),c("pan","gorilla","orang"))
  spec<-cvl[which(cvl[,3]==ind),1]
  tst<-sps[which(sps[,2]==spec),3]
  dfile=paste(seqpath,spec,"/",ind,"/mapping/merge/",ind,".merge.cram",sep="")
  hute<-system(paste("Rscript --vanilla contamination_test_ape.R ",gnom,",",tst,",",dfile," ", bcft," ",4,sep=""),intern=T)
  load(file=paste("summary_stats/hucontest.Robject",sep=""))
  huco[[ind]]<-hute
  save(huco,file=paste("summary_stats/hucontest.Robject",sep=""))

  print("done")
  print(Sys.time())
  
  q()
  
  }


