library(igraph)
library(dplyr)
library(ggplot2)
library(plyr)
rm(list = ls())

projectDir<-"/home/jupiter/PeaceCorps/CoralTda/"
projectDir<-"C:/Users/ert/Dropbox/Code/PeaceCorps/CoralTda/"

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


picDir<-paste0(projectDir,"figs")
printPlots=T
i=0
for(dat2s in list(dat2kernel,dat2single,dat2social)){
  i=i+1
  
  pl<-ggplot(data=dat2s,aes(x=bettiNumber,y=1-(t_coral/t_std),group=dataset,color=dataset))+ geom_line(size=2)+
    scale_x_discrete(name="Betti")+
    theme_minimal()+scale_y_continuous(name=("Time gain"))+
    theme(text = element_text(size=12),#legend.position = c(0.8, 0.8),
          legend.text = element_text(size=12))+
    #scale_colour_manual(name='', values=c('Large tokens'='blue','Small tokens'='red'))+
    guides(color = guide_legend(override.aes = list(size=5)));
  print(pl)
  
  p1<-ggplot(data=dat2s,aes(x=bettiNumber,y=V,group=dataset,color=dataset))+ geom_line(size=2)+
    scale_x_discrete(name="Betti")+
    theme_minimal()+scale_y_continuous(name=("Vertices"))+
    theme(text = element_text(size=12),legend.position = c(0.8, 0.8),
          legend.text = element_text(size=12))+
    #scale_colour_manual(name='', values=c('Large tokens'='blue','Small tokens'='red'))+
    guides(color = guide_legend(override.aes = list(size=5)));p1;
  p2<-ggplot(data=dat2s,aes(x=bettiNumber,y=E,group=dataset,color=dataset))+ geom_line(size=2)+
    scale_x_discrete(name="Betti")+
    theme_minimal()+scale_y_continuous(name=("Edges"))+
    theme(text = element_text(size=12),legend.position = c(0.8, 0.8),
          legend.text = element_text(size=12))+
    #scale_colour_manual(name='', values=c('Large tokens'='blue','Small tokens'='red'))+
    guides(color = guide_legend(override.aes = list(size=5)));p2;
  p3<-ggplot(data=dat2s,aes(x=bettiNumber,y=C,group=dataset,color=dataset))+ geom_line(size=2)+
    scale_x_discrete(name="Betti")+
    theme_minimal()+scale_y_continuous(name=("Clique size"))+
    theme(text = element_text(size=12),legend.position = c(0.8, 0.8),
          legend.text = element_text(size=12))+
    #scale_colour_manual(name='', values=c('Large tokens'='blue','Small tokens'='red'))+
    guides(color = guide_legend(override.aes = list(size=5)));p3 
  
  
  if(printPlots){
    ggsave(plot=pl,file=paste0("tCoral",i,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p1,file=paste0("vCoral",i,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p2,file=paste0("eCoral",i,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
    ggsave(plot=p3,file=paste0("sCoral",i,".png"),device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
  }
}

# do low edge graphs have betti>2

datsingleHigh<-datkernel
datsingleHigh<-datsocial
datsingleHigh$bettiNumber<-as.numeric(datsingleHigh$bettiNumber)
plot(datsingleHigh$EStd/datsingleHigh$VStd,datsingleHigh$bettiNumber)
datsingleHigh<-datsingleHigh[datsingleHigh$bettiNumber>1,]
plot(datsingleHigh$EStd/datsingleHigh$VStd,datsingleHigh$bettiNumber)
datsingleHigh
plot(datsingleHigh$EStd/datsingleHigh$VStd,datsingleHigh$BStd)
datsingleHigh$avgDegree<-datsingleHigh$EStd/datsingleHigh$VStd
p3<-ggplot(data=datsingleHigh,aes(x=avgDegree,y=BStd,group=dataset,color=dataset))+ geom_point(size=2)+
  scale_x_continuous(name="Average node degree in a single graph")+
  theme_minimal()+scale_y_continuous(name=("Value of Betti 2 and higher"))+
  theme(text = element_text(size=12),legend.position = c(0.8, 0.8),
        legend.text = element_text(size=12))+
  #scale_colour_manual(name='', values=c('Large tokens'='blue','Small tokens'='red'))+
  guides(color = guide_legend(override.aes = list(size=5)));p3 
max(datsingleHigh[datsingleHigh$BStd>0&datsingleHigh$avgDegree<7,]$bettiNumber)




