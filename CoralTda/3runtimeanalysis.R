rm(list = ls())
library(igraph)
library(TDA)
library(dplyr)

options(java.parameters = "-Xmx200g")

source("CoralTda/coralFunctions.R")

#  Graph ML datasets
dataPath <- "/home/jupiter/GraphML/"


loadGraph<-function(edgeFile){
  allData <- read.table(edgeFile, quote = "", comment.char = "", sep = " ")
  colnames(allData) <- c("source", "target")
  allData$source<-paste0("V",allData$source)
  allData$target<-paste0("V",allData$target)
  graph <- graph.data.frame(allData, directed = FALSE)
  return(graph)
}


compute<-function(dataPath, dataset,outputFile,rep=10,maxbetti=3){
  edgeFile <- paste0(dataPath, dataset, "edges")
  if (!file.exists(edgeFile)) {
    edgeFile=paste0(dataPath, dataset, "_A.txt")
  }
  graphFile <- paste0(dataPath, dataset, "graph_idx")
  if (!file.exists(graphFile)) {
    graphFile=paste0(dataPath, dataset, "_graph_indicator.txt")
  }
  message(graphFile)
  graphIdData <-
    read.table(graphFile,
               quote = "\"",
               comment.char = "",
               sep = ",")
  allData <-
    read.table(edgeFile,
               quote = "\"",
               comment.char = "",
               sep = ",")
  colnames(allData) <- c("source", "target")
  # I am converting graph vertex ids to strings that starts with V
  # This is done to avoid orig vertex ids where igraph given internal ids
  allData$source<-paste0("v",allData$source)
  allData$target<-paste0("v",allData$target)
  bettiFrame<-data.frame()
  graphs<-unique(graphIdData$V1)
  message("Processing ",dataset," has ",length(graphs)," graphs.")
  
  for(graphId in graphs){
    thisGraphNodes <- paste0("v",which(graphIdData$V1 == graphId))
    edgeData <-
      allData[allData$source %in% thisGraphNodes |
                allData$target %in% thisGraphNodes, ]
    graph <- graph.data.frame(edgeData, directed = FALSE)
    graph<-simplify(graph)
    #plot(graph)
    nodeCount <- vcount(graph)
    
    for(bettiTarget in seq(1,maxbetti,1)){
      #message(graphId,":",vcount(graph)," nodes, ",ecount(graph)," edges ",bettiTarget)
      graph2=graph
      results=computeTemporalTDA(graph2, bettiTarget,rep)
      #if(length(results)!=8)print(results)
      cliqStd = results[[1]]
      nodeStd = results[[2]]
      edgeStd = results[[3]]
      
      cliqCoral = results[[4]]
      nodeCoral = results[[5]]
      edgeCoral = results[[6]]
      
      bettiStd = results[[7]]
      bettiCoral = results[[8]]
      
      timeStd = results[[9]]
      timeCoral = results[[10]]
      str=paste0(dataset,"\t",graphId,"\t",bettiTarget, 
                 "\t",nodeStd,"\t", nodeCoral,"\t",edgeStd,"\t",edgeCoral,
                 "\t",cliqStd,
                 "\t",cliqCoral,
                 "\t",bettiStd,"\t",bettiCoral,"\t",timeStd,"\t",timeCoral)
      #cat(str,file=outputFile,append=TRUE)
      write(str, file=outputFile, append=T)
      if(bettiTarget>1&&bettiStd>0) {
        message(str)
      }
    }
  }
}
enzymes <-  "ENZYMES/ENZYMES."
proteins <-  "proteins/proteins."
redditbin <-  "REDDIT-BINARY/REDDIT-BINARY."
imdbm <-  "IMDB-MULTI/IMDB-MULTI."
nc1 <-  "NCI1/NCI1."
imdbbin <-  "IMDB-BINARY/IMDB-BINARY."
MSRC_21<-"MSRC_21/MSRC_21/MSRC_21"
DD<-"DD/DD/DD"
OHSU<-"OHSU/OHSU/OHSU"
FIRSTMM_DB<-"FIRSTMM_DB/FIRSTMM_DB/FIRSTMM_DB"
SYNTHETICnew<-"SYNTHETICnew/SYNTHETICnew/SYNTHETICnew"
SYNTHETIC<-"SYNTHETIC/SYNTHETIC/SYNTHETIC"




for(dataset in c(MSRC_21,DD,OHSU,FIRSTMM_DB,SYNTHETICnew,SYNTHETIC,redditbin)) {
  ind = gregexpr(pattern ='/',dataset)[[1]]-1
  outputFile = paste0("CoralTda/",substr(dataset,1,ind),"timeresults.csv")
  #Check if a previous result file  exists
  if (file.exists(outputFile)) {
    file.remove(outputFile)
  }
  compute(dataPath, dataset,outputFile,rep=30,maxbetti=5)
}


#---------------------------------

dataPathTwitter=paste0(dataPath,"twitter/")
filesTwitter=list.files(path = dataPathTwitter,recursive=T, pattern = ".edges", all.files = TRUE,
                        full.names = TRUE)

dataPathFB=paste0(dataPath,"facebook/")
filesFB=list.files(path = dataPathFB,recursive=T, pattern = ".edges", all.files = TRUE,
                   full.names = TRUE)

outputFile = "CoralTda/socialtimeresults.csv"#for twitter
#Check if a previous result file  exists
if (file.exists(outputFile)) {
  file.remove(outputFile)
}

maxbetti=5
for(fileName in c(filesFB,filesTwitter)){
  graph<-loadGraph(fileName)
  
  for(bettiTarget in seq(1,maxbetti,1)){
    graph2=graph
    
    if(ecount(graph2)<5000){
      
      results=computeTemporalTDA(graph2,bettiTarget,rep=1)
      cliqStd = results[[1]]
      nodeStd = results[[2]]
      edgeStd = results[[3]]
      
      cliqCoral = results[[4]]
      nodeCoral = results[[5]]
      edgeCoral = results[[6]]
      
      bettiStd = results[[7]]
      bettiCoral = results[[8]]
      
      timeStd = results[[9]]
      timeCoral = results[[10]]
      fname = substr(fileName,start=1+nchar(dataPath),stop=nchar(fileName)-5)
      str=paste0(fname,"\t",bettiTarget, 
                 "\t",nodeStd,"\t", nodeCoral,"\t",edgeStd,"\t",edgeCoral,
                 "\t",cliqStd,
                 "\t",cliqCoral,
                 "\t",bettiStd,"\t",bettiCoral,"\t",timeStd,"\t",timeCoral,"\r\n")
      cat(str,file=outputFile,append=TRUE)
      if(bettiStd>0) {
        message(str)
      }
      
    }
  }
}

## cora, citeseer, pubmed
maxbetti=10
outputSingleDatasetFile = paste0("CoralTda/singleDatasettimeresults.csv")
if (file.exists(outputSingleDatasetFile)) {
  file.remove(outputSingleDatasetFile)
}
for(dataset in c("citeseer/citeseer.","cora/cora.")){
  
  
edgeFile <- paste0(dataPath, dataset, "edges")

allData <-
  read.table(edgeFile,
             quote = "\"",
             comment.char = "",
             sep = ",")
colnames(allData) <- c("source", "target")
# I am converting graph vertex ids to strings that starts with V
# This is done to avoid orig vertex ids where igraph given internal ids
allData$source<-paste0("v",allData$source)
allData$target<-paste0("v",allData$target)
bettiFrame<-data.frame()

graph <- graph.data.frame(allData[,1:2], directed = FALSE)
graph<-simplify(graph)
nodeCount <- vcount(graph)
for(bettiTarget in seq(1,maxbetti,1)){
  graph2=graph
  
    
    results=computeTemporalTDA(graph2,bettiTarget,rep=10)
    cliqStd = results[[1]]
    nodeStd = results[[2]]
    edgeStd = results[[3]]
    
    cliqCoral = results[[4]]
    nodeCoral = results[[5]]
    edgeCoral = results[[6]]
    
    bettiStd = results[[7]]
    bettiCoral = results[[8]]
    
    timeStd = results[[9]]
    timeCoral = results[[10]]
    str=paste0(dataset,"\t",bettiTarget, 
               "\t",nodeStd,"\t", nodeCoral,"\t",edgeStd,"\t",edgeCoral,
               "\t",cliqStd,
               "\t",cliqCoral,
               "\t",bettiStd,"\t",bettiCoral,"\t",timeStd,"\t",timeCoral,"\r\n")
    cat(str,file=outputSingleDatasetFile,append=TRUE)
     
      message(str)
    
  }
}
