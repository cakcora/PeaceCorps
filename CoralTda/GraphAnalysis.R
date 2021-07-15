##############################################
# Clustering coefficient analysis
##############################################


library(dplyr)
library(ggplot2)
library(plyr)
library(stringr)

rm(list = ls())

projectDir<-"/home/jupiter/PeaceCorps/CoralTda/"
projectDir<-"C:/Users/ert/Dropbox/Code/PeaceCorps/CoralTda/"
projectDir<-"C:/Code/PeaceCorps/CoralTda/"

outputDir<-paste0(projectDir,"results/")

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

colnames(datkernel)<-c("Dataset","graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
datkernel$Dataset <- gsub('proteins/proteins.', 'PROTEINS', datkernel$Dataset)
datkernel$Dataset <- gsub('DD/DD/DD', 'DD', datkernel$Dataset)
datkernel$Dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY.', 'REDDITB', datkernel$Dataset)
datkernel$Dataset <- gsub('ENZYMES/ENZYMES.', 'ENZYMES', datkernel$Dataset)
datkernel$Dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernel$Dataset)
datkernel$Dataset <- gsub('OHSU/OHSU/OHSU', 'OHSU', datkernel$Dataset)
datkernel$Dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM', datkernel$Dataset)
datkernel$Dataset <- gsub('SYNTHETICnew/SYNTHETICnew/SYNTHETICnew', 'SynNew', datkernel$Dataset)
datkernel$Dataset <- gsub('SYNTHETIC/SYNTHETIC/SYNTHETIC', 'SYN', datkernel$Dataset)
datkernel$bettiNumber<-as.factor(datkernel$bettiNumber)

datsocial<-dplyr::bind_rows(social)
tmp<-str_split_fixed(datsocial$V1, "//", 2)
datsocial$V0<-tmp[,1]
datsocial$V1<-tmp[,2]
colnames(datsocial)<-c("graphId","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral","Dataset")
datsocial$bettiNumber<-as.factor(datsocial$bettiNumber)
datsocial$Dataset <- gsub('facebook', 'FACEBOOK', datsocial$Dataset)
datsocial$Dataset <- gsub('twitter', 'TWITTER', datsocial$Dataset)

singleDataset<-read.delim(paste0(outputDir,"singleDatasettimeresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingle<-(singleDataset)
colnames(datsingle)<-c("Dataset","bettiNumber","VStd","VCoral","EStd","ECoral","CStd","CCoral","BStd","Bcoral","timeStd","timeCoral")
datsingle$bettiNumber<-as.factor(datsingle$bettiNumber)
datsingle$Dataset <- gsub('citeseer/citeseer.', 'CITESEER', datsingle$Dataset)
datsingle$Dataset <- gsub('cora/cora.', 'CORA', datsingle$Dataset)


dat2kernel<-ddply(datkernel,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2social<-ddply(datsocial,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))
dat2single<-ddply(datsingle,.(Dataset,bettiNumber),summarize,t_std=sum(timeStd), t_coral=sum(timeCoral),E=mean(ECoral/EStd),V=mean(VCoral/VStd), C=mean(CCoral/CStd))

################################################################################

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

datkernelclus<-dplyr::bind_rows(enzyme,redditbinary,proteins,dd,firstmm_db,ohsu,nci1,syntheticnew,synthetic)

colnames(datkernelclus)<-c("Dataset","graphId","cluscoeff","V","E")
datkernelclus$Dataset <- gsub('DD/DD/DD', 'DD', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('ENZYMES/ENZYMES.', 'ENZYMES', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB', 'FIRSTMM', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('MSRC_21/MSRC_21/MSRC_21', 'MSRC_21', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('NCI1/NCI1.', 'NCI1', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('OHSU/OHSU/OHSU', 'OHSU', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('proteins/proteins.', 'PROTEINS', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('REDDIT-BINARY/REDDIT-BINARY.', 'REDDITB', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('SYNTHETICnew/SYNTHETICnew/SYNTHETICnew', 'SYNNEW', datkernelclus$Dataset)
datkernelclus$Dataset <- gsub('SYNTHETIC/SYNTHETIC/SYNTHETIC', 'SYN', datkernelclus$Dataset)

datsocialclus<-dplyr::bind_rows(social)
tmp<-str_split_fixed(datsocialclus$V1, "//", 2)
datsocialclus$V0<-tmp[,1]
datsocialclus$V1<-tmp[,2]
colnames(datsocialclus)<-c("graphId","cluscoeff","V","E","Dataset")
datsocialclus$Dataset <- gsub('facebook', 'FACEBOOK', datsocialclus$Dataset)
datsocialclus$Dataset <- gsub('twitter', 'TWITTER', datsocialclus$Dataset)

singleDataset<-read.delim(paste0(outputDir,"singleDatasetcluscoeffresults.csv"),sep="\t",header=F,stringsAsFactors =FALSE)
datsingleclus<-(singleDataset)
colnames(datsingleclus)<-c("Dataset","cluscoeff","V","E")
datsingleclus$Dataset <- gsub('citeseer/citeseer.', 'CITESEER', datsingleclus$Dataset)
datsingleclus$Dataset <- gsub('cora/cora.', 'CORA', datsingleclus$Dataset)


dat3kernel<-ddply(datkernelclus,.(Dataset),summarize,E=mean(E),V=mean(V), C=mean(na.omit(cluscoeff)))
dat3social<-ddply(datsocialclus,.(Dataset),summarize,E=mean(E),V=mean(V), C=mean(cluscoeff))
dat3single<-ddply(datsingleclus,.(Dataset),summarize,E=mean(E),V=mean(V), C=mean(cluscoeff))

kernelchars<-merge(datkernel,datkernelclus, by=c("Dataset","graphId"))
picDir<-paste0(projectDir,"figs")

for(betti in c(1,2,3,4,5)){
  p1=ggplot(data=kernelchars[kernelchars$bettiNumber==betti,], aes(x=cluscoeff,y=BStd))+geom_point(color="blue")+
    labs(x ="Clustering coeff", y = paste0("Betti ",betti))+theme_minimal()+
    theme(text = element_text(size=22),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));
  
  p2=(ggplot(data=kernelchars[kernelchars$bettiNumber==betti,], aes(x=E/V,y=BStd))+geom_point(color="blue")+
        theme(text = element_text(size=22),
              legend.text = element_text(size=22))+theme_minimal()+
        labs(x ="Mean degree", y = paste0("Betti ",betti)))
  
  ggsave(plot=p1,file=paste0("clusCoralKernel",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  ggsave(plot=p2,file=paste0("avgCoralKernel",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  
}

socialchars<-merge(datsocial,datsocialclus, by=c("Dataset","graphId"))

for(betti in c(1,2,3,4,5)){
  p3=ggplot(data=socialchars[socialchars$bettiNumber==betti,], aes(x=cluscoeff,y=BStd))+geom_point(color="blue")+
    labs(x ="Clustering coeff", y = paste0("Betti ",betti))+theme_minimal()+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));p3;

    p4=(ggplot(data=socialchars[socialchars$bettiNumber==betti,], aes(x=E/V,y=BStd))+geom_point(color="blue")+
        labs(x ="mean degree", y = paste0("Betti ",betti)))+theme_minimal()+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));
    
  print(p3)
  print(p4)
  ggsave(plot=p3,file=paste0("clusCoralSocial",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  ggsave(plot=p4,file=paste0("avgCoralSocial",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  
}

singlechars<-merge(datsingle,datsingleclus, by=c("Dataset"))

for(betti in c(1,2,3,4,5)){
  p5=ggplot(data=singlechars[singlechars$bettiNumber==betti,], aes(x=cluscoeff,y=BStd))+geom_point(color="blue")+
    labs(x ="Clustering coeff", y = paste0("Betti ",betti))+theme_minimal()+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));
  
  p6=(ggplot(data=singlechars[singlechars$bettiNumber==betti,], aes(x=E/V,y=BStd))+geom_point(color="blue")+
        labs(x ="mean degree", y = paste0("Betti ",betti)))+theme_minimal()+
    theme(text = element_text(size=22),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size=5)));
  print(p5)
  print(p6)
  ggsave(plot=p5,file=paste0("clusCoralSingle",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  ggsave(plot=p6,file=paste0("avgCoralSingle",betti,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  
}



