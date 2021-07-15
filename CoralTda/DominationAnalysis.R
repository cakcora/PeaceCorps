
# this script requires results from coralTDADominatingSet
# We plot reductions in vertex.
library(igraph)
library(dplyr)
library(ggplot2)
library(plyr)
library(stringr)
rm(list = ls())

projectDir<-"/home/jupiter/PeaceCorps/CoralTda/"
projectDir<-"C:/Code/PeaceCorps/CoralTda/"

outputDir<-paste0(projectDir,"results/")


datkernel1 <- read.delim(paste0(outputDir,"DominationKernelResults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)


colnames(datkernel1)<-c("Dataset","graphId","V","E","Reduction")
datkernel1$Dataset <- gsub('DD/DD', 'DD', datkernel1$Dataset)
datkernel1$Dataset <- gsub('ENZYMES/ENZYMES', 'ENZYMES', datkernel1$Dataset)
datkernel1$Dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM', datkernel1$Dataset)
datkernel1$Dataset <- gsub('MSRC_21/MSRC_21', 'MSRC', datkernel1$Dataset)
datkernel1$Dataset <- gsub('NCI1/NCI1', 'NCI1', datkernel1$Dataset)
datkernel1$Dataset <- gsub('OHSU/OHSU', 'OHSU', datkernel1$Dataset)

datkernel1$Dataset <- gsub('proteins/proteins', 'PROTEINS', datkernel1$Dataset)
datkernel1$Dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY/REDDIT-BINARY', 'REDDITB', datkernel1$Dataset)
datkernel1$Dataset <- gsub('SYNTHETICnew/SYNTHETICnew', 'SYNNEW', datkernel1$Dataset)
datkernel1$Dataset <- gsub('SYNTHETIC/SYNTHETIC', 'SYN', datkernel1$Dataset)

unique(datkernel1$Dataset)

datsocial<-read.delim(paste0(outputDir,"DominationSocialResults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
tmp<-str_split_fixed(datsocial$V1, "C://data/", 2)
datsocial$V1<-tmp[,2]
colnames(datsocial)<-c("Dataset","V","E","Reduction")
#datsocial$bettiNumber<-as.factor(datsocial$bettiNumber)
datsocial$Dataset <- gsub('facebook/', 'FB_', datsocial$Dataset)
datsocial$Dataset <- gsub('twitter/', 'TWITTER_', datsocial$Dataset)
datsocial$Dataset <- gsub('.edges', '', datsocial$Dataset)
tmp<-str_split_fixed(datsocial$Dataset, "_", 2)
datsocial$Dataset=tmp[,1]
datsocial$graphId=tmp[,2]
unique(datsocial$Dataset)

singledataset<-read.delim(paste0(outputDir,"DominationSingleResults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingle<-(singledataset)
colnames(datsingle)<-c("Dataset","V","E","Reduction")
datsingle$Dataset <- gsub('citeseer/citeseer.', 'CITESEER', datsingle$Dataset)
datsingle$Dataset <- gsub('cora/cora.', 'CORA', datsingle$Dataset)
unique(datsingle$Dataset)

dat2kernel1<-ddply(datkernel1,.(Dataset),summarize,meanReduction=mean(Reduction))
dat2social<-ddply(datsocial,.(Dataset),summarize,meanReduction=mean(Reduction))
dat2single<-ddply(datsingle,.(Dataset),summarize,meanReduction=mean(Reduction))

datall<-dplyr::bind_rows(dat2kernel1,dat2social,dat2single)
datall<-datall[order(datall$meanReduction,decreasing = T),]
rownames(datall)<-NULL
positions <- c(unique(datall$Dataset))
picDir<-paste0(projectDir,"figs")
 
  #aDataset$bettiNumber=as.numeric(aDataset$bettiNumber)
  pl<-ggplot(data=datall,aes(x=Dataset,y=(meanReduction)*100,order=meanReduction),color="red")+ geom_col(size=2,fill="lightblue")+
    theme_minimal()+scale_y_continuous(name=("Vertex reduction (%)"))+
    scale_x_discrete(limits=positions)+
    #geom_text(aes(label = 100*round(meanReduction,2)), vjust = -0.2,size=5)+
    theme(text = element_text(size=22,angle = 90, vjust = 0.5, hjust=0),
          axis.title.x=element_blank(),axis.text.x = element_text(color="black", size=15,margin=margin(-12,0,0,0)));pl
  
  
  ggsave(plot=pl,file=paste0("dominatingCoral.png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  


