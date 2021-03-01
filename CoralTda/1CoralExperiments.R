library(igraph)
library(TDA)
library(dplyr)
rm(list = ls())
options(java.parameters = "-Xmx200g")


#  Graph ML datasets
dataPath <- "C:/data/"

loadGraph<-function(dataset, dataAlias){
  whichOutputFile<-paste0(outputPath,dataAlias,outputFile)
  # Remove earlier result files, if they exist
  if (file.exists(whichOutputFile)) {
    file.remove(whichOutputFile)
  }
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
  graphs <- unique(graphIdData$V1)
  message("Processing ",dataAlias," for feature ",feature,". Num of Graphs: ",length(graphs))
  results<-data.frame()
  for (graphId in graphs) {
    # locate nodes of the graph
    thisGraphNodes <- which(graphIdData$V1 == graphId)
    # load edges of the graph
    edgeData <-
      allData[allData$source %in% thisGraphNodes |
                allData$target %in% thisGraphNodes, ]
    graph <- graph.data.frame(edgeData, directed = FALSE)
    nodeCount <- vcount(graph)
     
      # below we use sub level filtrations.
      
      #sublevel filtration
      # compute node function values 
      tryCatch({
        F.values = computeNodeVals(feature, graph,dataset=dataset,dataPath=dataPath)
      }, error = function(error_condition) {
        message("Error in extracting PDs for ",dataAlias)
        return(-1);
      })
      # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
      cmplx <- cliques(as.undirected(graph), min = 0, max = maxClique)
      # use sublevel=T for sublevel, sublevel=F for superlevel filtration
      # F.values are node values. At these values, node appear in the complex,
      # and their edges to other active nodes also activate.
      FltRips <- funFiltration(FUNvalues = F.values,
                               cmplx = cmplx,
                               sublevel = T) # Construct filtration using F.values
      
      #extract the persistence diagram
      persistenceDiagram <-
        filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
      
      #extract betti number - unused
      # b0OrigGraph <- sum(persistenceDiagram[, 1] == 0)
      # b1OrigGraph <- sum(persistenceDiagram[, 1] == 1)
      # b2OrigGraph <- sum(persistenceDiagram[, 1] == 2)
      
      #extract birth and death times
      pd2 <-cbind(persistenceDiagram,GraphId = graphId,dataset = dataAlias)
      results<-rbind(results,pd2)
      
    
  }
  write.table(
    results,
    file = whichOutputFile,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    append = T,
    quote = FALSE
  )
}


loadGraph("NCI1/NCI1.","NCI1",feature=f)
loadGraph("ENZYMES/ENZYMES.","Enzyme",feature=f)
loadGraph("BZR/BZR.","BZR",feature=f)
loadGraph("COX2/COX2.","COX2",feature=f)
loadGraph("DHFR/DHFR.","DHFR",feature=f)
loadGraph("FRANKENSTEIN/FRANKENSTEIN.","FRANKENSTEIN",feature=f)
#extractPD("REDDIT-MULTI-5K/REDDIT-MULTI-5K.","REDDIT5K",feature=f)

loadGraph("proteins/proteins.","Protein",feature=f)
loadGraph("REDDIT-BINARY/REDDIT-BINARY.","RedditBinary",feature=f)
loadGraph("IMDB-MULTI/IMDB-MULTI.","IMDBMulti",feature=f)

loadGraph("IMDB-BINARY/IMDB-BINARY.","IMDBBinary",feature=f)