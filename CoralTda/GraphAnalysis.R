library(dplyr)
library(ggplot2)
library(plyr)
rm(list = ls())

projectDir<-"/home/jupiter/PeaceCorps/CoralTda/"
projectDir<-"C:/Users/ert/Dropbox/Code/PeaceCorps/CoralTda/"
projectDir<-"C:/Code/PeaceCorps/CoralTda/"

outputDir<-paste0(projectDir,"results/")
outputDir<-"C:/data/tdacoral/results/"
##############################################
#Graph and Betti analysis
##############################################
social <- read.delim(paste0(outputDir,"socialtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
enzyme <- read.delim(paste0(outputDir,"ENZYMEStimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
proteins <- read.delim(paste0(outputDir,"proteinstimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
redditbinary <- read.delim(paste0(outputDir,"REDDIT-BINARYtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
nci1 <- read.delim(paste0(outputDir,"NCI1timeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
msrc_21<-read.delim(paste0(outputDir,"MSRC_21timeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
dd<-read.delim(paste0(outputDir,"DDtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
firstmm_db<-read.delim(paste0(outputDir,"FIRSTMM_DBtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
ohsu<-read.delim(paste0(outputDir,"OHSUtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
synthetic<-read.delim(paste0(outputDir,"SYNTHETICtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
syntheticnew<-read.delim(paste0(outputDir,"SYNTHETICnewtimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)


datkernel<-dplyr::bind_rows(enzyme,redditbinary,proteins,dd,firstmm_db,ohsu,nci1,msrc_21,nci1,syntheticnew,synthetic)

colnames(datkernel)<-c("dataset","graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
datkernel$dataset <- gsub('proteins/proteins.', 'Proteins', datkernel$dataset)
datkernel$dataset <- gsub('DD/DD/DD', 'DD', datkernel$dataset)
datkernel$dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY.', 'Redditb', datkernel$dataset)
datkernel$dataset <- gsub('ENZYMES/ENZYMES.', 'Enzymes', datkernel$dataset)
datkernel$dataset <- gsub('MSRC_21/MSRC_21/MSRC_21', 'Msrc', datkernel$dataset)
datkernel$dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernel$dataset)
datkernel$dataset <- gsub('OHSU/OHSU/OHSU', 'OHSU', datkernel$dataset)
datkernel$dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM_DB', datkernel$dataset)
datkernel$dataset <- gsub('SYNTHETICnew/SYNTHETICnew/SYNTHETICnew', 'SynNew', datkernel$dataset)
datkernel$dataset <- gsub('SYNTHETIC/SYNTHETIC/SYNTHETIC', 'Syn', datkernel$dataset)
datkernel$dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernel$dataset)
datkernel$bettiNumber<-as.factor(datkernel$bettiNumber)

datsocial<-dplyr::bind_rows(social)
library(stringr)
tmp<-str_split_fixed(datsocial$V1, "//", 2)
datsocial$V0<-tmp[,1]
datsocial$V1<-tmp[,2]
colnames(datsocial)<-c("graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral","dataset")
datsocial$bettiNumber<-as.factor(datsocial$bettiNumber)
datsocial$dataset <- gsub('facebook', 'FB', datsocial$dataset)
datsocial$dataset <- gsub('twitter', 'Twitter', datsocial$dataset)

singledataset<-read.delim(paste0(outputDir,"singleDatasettimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingle<-(singledataset)
colnames(datsingle)<-c("dataset","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
datsingle$bettiNumber<-as.factor(datsingle$bettiNumber)
datsingle$dataset <- gsub('citeseer/citeseer.', 'citeseer', datsingle$dataset)
datsingle$dataset <- gsub('cora/cora.', 'cora', datsingle$dataset)


dat2kernel<-ddply(datkernel,.(dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2social<-ddply(datsocial,.(dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2single<-ddply(datsingle,.(dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))


##############################################
# Clustering coefficient analysis
##############################################
social <- read.delim(paste0(outputDir,"socialcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
enzyme <- read.delim(paste0(outputDir,"ENZYMEScluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
proteins <- read.delim(paste0(outputDir,"proteinscluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
redditbinary <- read.delim(paste0(outputDir,"REDDIT-BINARYcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
nci1 <- read.delim(paste0(outputDir,"NCI1cluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
msrc_21<-read.delim(paste0(outputDir,"MSRC_21cluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
dd<-read.delim(paste0(outputDir,"DDcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
firstmm_db<-read.delim(paste0(outputDir,"FIRSTMM_DBcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
ohsu<-read.delim(paste0(outputDir,"OHSUcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
synthetic<-read.delim(paste0(outputDir,"SYNTHETICcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
syntheticnew<-read.delim(paste0(outputDir,"SYNTHETICnewcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)


datkernelclus<-dplyr::bind_rows(enzyme,redditbinary,proteins,dd,firstmm_db,ohsu,nci1,msrc_21,nci1,syntheticnew,synthetic)

colnames(datkernelclus)<-c("dataset","graphId","cluscoeff","V","E")
datkernelclus$dataset <- gsub('proteins/proteins.', 'Proteins', datkernelclus$dataset)
datkernelclus$dataset <- gsub('DD/DD/DD', 'DD', datkernelclus$dataset)
datkernelclus$dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY.', 'Redditb', datkernelclus$dataset)
datkernelclus$dataset <- gsub('ENZYMES/ENZYMES.', 'Enzymes', datkernelclus$dataset)
datkernelclus$dataset <- gsub('MSRC_21/MSRC_21/MSRC_21', 'Msrc', datkernelclus$dataset)
datkernelclus$dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernelclus$dataset)
datkernelclus$dataset <- gsub('OHSU/OHSU/OHSU', 'OHSU', datkernelclus$dataset)
datkernelclus$dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM_DB', datkernelclus$dataset)
datkernelclus$dataset <- gsub('SYNTHETICnew/SYNTHETICnew/SYNTHETICnew', 'SynNew', datkernelclus$dataset)
datkernelclus$dataset <- gsub('SYNTHETIC/SYNTHETIC/SYNTHETIC', 'Syn', datkernelclus$dataset)
datkernelclus$dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernelclus$dataset)

datsocialclus<-dplyr::bind_rows(social)
library(stringr)
tmp<-str_split_fixed(datsocialclus$V1, "//", 2)
datsocialclus$V0<-tmp[,1]
datsocialclus$V1<-tmp[,2]
colnames(datsocialclus)<-c("graphId","cluscoeff","V","E","dataset")
datsocialclus$dataset <- gsub('facebook', 'FB', datsocialclus$dataset)
datsocialclus$dataset <- gsub('twitter', 'Twitter', datsocialclus$dataset)

singledataset<-read.delim(paste0(outputDir,"singleDatasetcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingleclus<-(singledataset)
colnames(datsingleclus)<-c("dataset","cluscoeff","V","E")
datsingleclus$dataset <- gsub('citeseer/citeseer.', 'citeseer', datsingleclus$dataset)
datsingleclus$dataset <- gsub('cora/cora.', 'cora', datsingleclus$dataset)


dat3kernel<-ddply(datkernelclus,.(dataset),summarize,E=mean(E),V=mean(V), C=mean(na.omit(cluscoeff)))
dat3social<-ddply(datsocialclus,.(dataset),summarize,E=mean(E),V=mean(V), C=mean(cluscoeff))
dat3single<-ddply(datsingleclus,.(dataset),summarize,E=mean(E),V=mean(V), C=mean(cluscoeff))

kernelchars<-merge(datkernel,datkernelclus, by=c("dataset","graphId"))
picDir<-paste0(projectDir,"figs")

for(betti in c(1,2,3,4,5)){
  pl=ggplot(data=kernelchars[kernelchars$bettiNumber==betti,], aes(x=cluscoeff,y=BStd))+geom_point()+labs(x ="clustering coeff", y = paste0("Betti",betti))
  p1=(ggplot(data=kernelchars[kernelchars$bettiNumber==betti,], aes(x=E/V,y=BStd))+geom_point()+labs(x ="E/V", y = paste0("Betti",betti)))
  
  ggsave(plot=pl,file=paste0("clusCoral",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  ggsave(plot=p1,file=paste0("avgCoral",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  
}
write.csv(kernelchars,file="allclus.csv",quote=F)

