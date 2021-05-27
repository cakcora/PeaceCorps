rm(list = ls())
library(igraph)
source("CoralTda/coralFunctions.R")
compute<-function(dataPath, dataset,outputFile){
  edgeFile <- paste0(dataPath, dataset, "edges")
  graphFile <- paste0(dataPath, dataset, "graph_idx")
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
    if(ecount(graph)>2000){
      message(graphId," will not be computed:",ecount(graph))
    }
    else {
      for(bettiTarget in seq(1,10,1)){
        #message(graphId,":",vcount(graph)," nodes, ",ecount(graph)," edges ",bettiTarget)
        graph2=graph
        results=computeCoralTDA(graph2, bettiTarget)
        if(length(results)!=8)print(results)
        cliqStd = results[[1]]
        nodeStd = results[[2]]
        edgeStd = results[[3]]
        
        cliqCoral = results[[4]]
        nodeCoral = results[[5]]
        edgeCoral = results[[6]]
        
        bettiStd = results[[7]]
        #print(bettiStd)
        bettiCoral = results[[8]]
        str=paste0(dataset,"\t",graphId,"\t",bettiTarget, 
                   "\t",nodeStd,"\t", nodeCoral,"\t",edgeStd,"\t",edgeCoral,
                   "\t",cliqStd,
                   "\t",cliqCoral,
                   "\t",bettiStd,"\t",bettiCoral)
        #cat(str,file=outputFile,append=TRUE)
        write(str, file=outputFile, append=T)
        if(bettiTarget>1&&bettiStd>0) {
          message(str)
        }
      }
    }
  }
}
dataPath <- "/home/jupiter/GraphML/"
enzymes <-  "ENZYMES/ENZYMES."
proteins <-  "proteins/proteins."
redditbin <-  "REDDIT-BINARY/REDDIT-BINARY."
imdbm <-  "IMDB-MULTI/IMDB-MULTI."
nc1 <-  "NCI1/NCI1."
imdbbin <-  "IMDB-BINARY/IMDB-BINARY."


for(dataset in c(nc1,imdbbin)) {
  ind = gregexpr(pattern ='/',dataset)[[1]]-1
  outputFile = paste0("CoralTda/coralTDA",substr(dataset,1,ind),"results.csv")
  #Check if a previous result file  exists
  if (file.exists(outputFile)) {
    file.remove(outputFile)
  }
  compute(dataPath, dataset,outputFile)
}

