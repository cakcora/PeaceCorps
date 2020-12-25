# requires results from GraphPersistenceDiagramExtraction.R
# This file contains code for exploratory analysis only
library(dplyr)
library(tsne)
library(ggplot2)
library(plyr)
library(igraph)
data<-read.table(file="graphMlResultsAnalysisSub.txt",header=F,sep="\t")
data<-na.omit(data)
colnames(data)<-c("dataset","graphId","betti","bettisignature")
hist(data$betti)

maxLength <- 0
vectorLength <- 300
discarded<-0
bettiHistogramData<-array()
sigVectors<-data.frame()
for(row in seq(seq_len(nrow(data)))){
  betti<-(data[row,]$betti)
  datasetName<-data[row,]$dataset
  graphId<-(data[row,]$graphId)
  bfunct<-(data[row,]$bettisignature)
  value <- strsplit(bfunct, split = " ")[[1]]
  sigArray <- value[1:vectorLength]
  sigArray[is.na(sigArray)]<-0
  sigArray<-as.integer(sigArray)

  lengthResult <- length(value)
  if(lengthResult>maxLength){
    maxLength <- lengthResult
  }
  if(lengthResult >vectorLength){
    discarded<-discarded+1
  }
  bettiHistogramData<-c(bettiHistogramData, lengthResult)
  df<- cbind.data.frame(name=datasetName,graphId=as.integer(graphId),betti=as.integer(betti))
  for(ind in seq(1:vectorLength)){
    df<-cbind(df,sigArray[ind])
  }
  sigVectors<-rbind(sigVectors, df)

}
message("We have at most ",maxLength," betti values in a signature")
message("if we use ",vectorLength," length, we ignore some signature values in ",discarded," graphs"," of total ", nrow(data))
hist(bettiHistogramData)


plot(sigVectors[sigVectors$betti==1,])
distw = dist(sigVectors[,4:(3+vectorLength)])
s<-list()
for(perp in c(5,10,15,20,25,30,50)){
  tsneResults = tsne(distw, epoch_callback = NULL, perplexity=perp, max_iter = 2500)
  s[[perp]]<-addresses
}
ts<-tsneResults
for(perp in c(5,10,15,20,25,30,50)){
  tsneResults<-s[[perp]]
  ggplot()+geom_point(aes(x=tsneResults[,1], y=tsneResults[,2],color=sigVectors$name, shape=sigVectors$name), size=3)
  p=p1+theme_classic()+labs(x="TSNE1",y="TSNE2")+
    theme(text = element_text(size=16),legend.title = element_blank(),legend.position = "none")
  show(p)
}
p=p1+theme_classic()+labs(x="TSNE1",y="TSNE2")+
  theme(text = element_text(size=16),legend.title = element_blank())
show(p)
ggsave(plot=p,file="tsnelegend.eps",width=6,height=4,units="in",path=imgDir)

