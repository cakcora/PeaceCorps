
# this script requires results from coralTDA runs.
# We plot reductions in vertex, edge, clique size and time.
library(igraph)
library(dplyr)
library(ggplot2)
library(plyr)
library(stringr)
rm(list = ls())

projectDir<-"/home/jupiter/PeaceCorps/CoralTda/"
projectDir<-"C:/Code/PeaceCorps/CoralTda/"

outputDir<-paste0(projectDir,"results/")


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


datkernel1<-dplyr::bind_rows(dd,enzyme,firstmm_db,nci1,ohsu)
datkernel2<-dplyr::bind_rows(proteins,redditbinary,syntheticnew,synthetic)

colnames(datkernel1)<-c("Dataset","graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
colnames(datkernel2)<-c("Dataset","graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
datkernel1$Dataset <- gsub('DD/DD/DD', 'DD', datkernel1$Dataset)
datkernel1$Dataset <- gsub('ENZYMES/ENZYMES.', 'ENZYMES', datkernel1$Dataset)
datkernel1$Dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM', datkernel1$Dataset)
datkernel1$Dataset <- gsub('MSRC_21/MSRC_21/MSRC_21', 'MSRC', datkernel1$Dataset)
datkernel1$Dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernel1$Dataset)
datkernel1$Dataset <- gsub('OHSU/OHSU/OHSU', 'OHSU', datkernel1$Dataset)

datkernel2$Dataset <- gsub('proteins/proteins.', 'PROTEINS', datkernel2$Dataset)
datkernel2$Dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY.', 'REDDITB', datkernel2$Dataset)
datkernel2$Dataset <- gsub('SYNTHETICnew/SYNTHETICnew/SYNTHETICnew', 'SYNNEW', datkernel2$Dataset)
datkernel2$Dataset <- gsub('SYNTHETIC/SYNTHETIC/SYNTHETIC', 'SYN', datkernel2$Dataset)
#datkernel2$bettiNumber<-as.factor(datkernel2$bettiNumber)

datsocial<-dplyr::bind_rows(social)
tmp<-str_split_fixed(datsocial$V1, "//", 2)
datsocial$V0<-tmp[,1]
datsocial$V1<-tmp[,2]
colnames(datsocial)<-c("graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral","Dataset")
#datsocial$bettiNumber<-as.factor(datsocial$bettiNumber)
datsocial$Dataset <- gsub('facebook', 'FACEBOOK', datsocial$Dataset)
datsocial$Dataset <- gsub('twitter', 'TWITTER', datsocial$Dataset)

singledataset<-read.delim(paste0(outputDir,"singleDatasettimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingle<-(singledataset)
colnames(datsingle)<-c("Dataset","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
#datsingle$bettiNumber<-as.factor(datsingle$bettiNumber)
datsingle$Dataset <- gsub('citeseer/citeseer.', 'CITESEER', datsingle$Dataset)
datsingle$Dataset <- gsub('cora/cora.', 'CORA', datsingle$Dataset)
datsingle<-datsingle[datsingle$bettiNumber<=5,]

dat2kernel1<-ddply(datkernel1,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2kernel2<-ddply(datkernel2,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2social<-ddply(datsocial,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2single<-ddply(datsingle,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))

datr<-dplyr::bind_rows(dat2kernel1,dat2kernel2,dat2social,dat2single)
datr 

dat2s<-ddply(datr,.(bettiNumber),summarize,eReduction=1-mean(E), vReduction=1-mean(V),vMinReduction=1-max(V),vMaxReduction=1-min(V))

picDir<-paste0(projectDir,"figs")
pictureIndex=0
for(aDataset in list(dat2kernel1,dat2kernel2,dat2single,dat2social)){
  pictureIndex=pictureIndex+1
  #aDataset$bettiNumber=as.numeric(aDataset$bettiNumber)
  pl<-ggplot(data=aDataset,aes(x=bettiNumber,y=(1-(t_coral/t_std))*100,group=Dataset,color=Dataset))+ geom_line(size=2)+
    scale_x_continuous(name="Betti dim")+
    theme_minimal()+scale_y_continuous(name=("Time reduction (%)"))+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));
  print(pl)
  
  p2<-ggplot(data=aDataset,aes(x=bettiNumber,y=(1-V)*100,group=Dataset,color=Dataset))+ geom_line(size=2)+
    scale_x_continuous(name="Betti dim")+
    theme_minimal()+scale_y_continuous(name=("Vertex reduction (%)"))+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));p2;
  
  p3<-ggplot(data=aDataset,aes(x=bettiNumber,y=(1-E)*100,group=Dataset,color=Dataset))+ geom_line(size=2)+
    scale_x_continuous(name="Betti dim")+
    theme_minimal()+scale_y_continuous(name=("Edge reduction (%)"))+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));p3;
  
  p4<-ggplot(data=aDataset,aes(x=bettiNumber,y=(1-C)*100,group=Dataset,color=Dataset))+ geom_line(size=2)+
    scale_x_continuous(name="Betti dim")+
    theme_minimal()+scale_y_continuous(name=("Clique count reduction (%)"))+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));p4 
  
  
  
    ggsave(plot=pl,file=paste0("tCoral",pictureIndex,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p2,file=paste0("vCoral",pictureIndex,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p3,file=paste0("eCoral",pictureIndex,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p4,file=paste0("sCoral",pictureIndex,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  
}


