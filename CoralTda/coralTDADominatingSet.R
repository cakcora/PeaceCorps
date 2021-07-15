rm(list = ls())
library(igraph)
library(TDA)
library(dplyr)

options(java.parameters = "-Xmx200g")

#  Graph ML datasets
projectDir<-"C:/Code/PeaceCorps/CoralTda/"
source("CoralTda/coralFunctions.R")
dataPath <- "C://data/"



loadGraph<-function(edgeFile){
  allData <- read.table(edgeFile, quote = "", comment.char = "", sep = " ")
  colnames(allData) <- c("source", "target")
  allData$source<-paste0("V",allData$source)
  allData$target<-paste0("V",allData$target)
  graph <- graph.data.frame(allData, directed = FALSE)
  return(graph)
}


computeDominatingRatio<-function(dataPath, dataset,outputFile){
  edgeFile <- paste0(dataPath, dataset, "edges")
  if (!file.exists(edgeFile)) {
    edgeFile=paste0(dataPath, dataset, "_A.txt")
  }
  graphFile <- paste0(dataPath, dataset, "graph_idx")
  if (!file.exists(graphFile)) {
    graphFile=paste0(dataPath, dataset, "_graph_indicator.txt")
  }
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
    as_adjacency_matrix(graph)
    nodeCount <- vcount(graph)
    mat = as_adjacency_matrix(graph,sparse=F)
    ratio = dominationRatio(mat)
    str=paste0(dataset,"\t",graphId,"\t",vcount(graph),"\t",ecount(graph),"\t",ratio,"\r")
    cat(str,file=outputFile,append=TRUE)
  }
}
enzymes <-  "ENZYMES/ENZYMES"
proteins <-  "proteins/proteins"
redditbin <-  "REDDIT-BINARY/REDDIT-BINARY/REDDIT-BINARY"
imdbm <-  "IMDB-MULTI/IMDB-MULTI"
nc1 <-  "NCI1/NCI1"
imdbbin <-  "IMDB-BINARY/IMDB-BINARY"
MSRC_21<-"MSRC_21/MSRC_21"
DD<-"DD/DD"
OHSU<-"OHSU/OHSU"
FIRSTMM_DB<-"FIRSTMM_DB/FIRSTMM_DB"
SYNTHETICnew<-"SYNTHETICnew/SYNTHETICnew"
SYNTHETIC<-"SYNTHETIC/SYNTHETIC"



outputFile = paste0("CoralTda/results/DominationKernelResults.csv")
  #Check if a previous result file  exists
  if (file.exists(outputFile)) {
    file.remove(outputFile)
  }

# processed: proteins,enzymes,redditbin,nc1,MSRC_21,
for(dataset in c(DD,OHSU,FIRSTMM_DB,SYNTHETICnew,SYNTHETIC)) {
  computeDominatingRatio(dataPath, dataset,outputFile)
}


dataPathTwitter=paste0(dataPath,"twitter/")
filesTwitter=list.files(path = dataPathTwitter,recursive=T, pattern = ".edges", all.files = TRUE,
                        full.names = TRUE)

dataPathFB=paste0(dataPath,"facebook/")
filesFB=list.files(path = dataPathFB,recursive=T, pattern = ".edges", all.files = TRUE,
                   full.names = TRUE)

outputFile = "CoralTda/results/DominationSocialResults.csv"
#Check if a previous result file  exists
if (file.exists(outputFile)) {
  file.remove(outputFile)
}

for(fileName in c(filesTwitter,filesFB)){
  graph<-loadGraph(fileName)
  graph<-simplify(graph)
  mat = as_adjacency_matrix(graph,sparse=F)
  if(!isSymmetric(mat)){
    message("Error: The adj matrix is not symmetric.")
  }else{
  ratio = dominationRatio(mat)
  str=paste0(fileName,"\t",vcount(graph),"\t",ecount(graph),"\t",ratio,"\r")
  cat(str,file=outputFile,append=TRUE)
  }
}

## cora, citeseer, pubmed
outputSingleDatasetFile = paste0("CoralTda/results/DominationSingleResults.csv")
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
  
  graph <- graph.data.frame(allData[,1:2], directed = FALSE)
  graph<-simplify(graph)
  mat = as_adjacency_matrix(graph,sparse=F)
  ratio = dominationRatio(mat)
  str=paste0(dataset,"\t",vcount(graph),"\t",ecount(graph),"\t",ratio,"\r")
  cat(str,file=outputSingleDatasetFile,append=TRUE)
  
}

#small running example for domination:

for(gInd in seq(5,50)){
  g <- sample_gnp(gInd, 3/10)
  mat = as_adjacency_matrix(g,sparse=F)
  ratio = dominationRatio(mat)
  message(gInd," ",ratio)
}
