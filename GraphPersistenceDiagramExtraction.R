# Title     : GraphPersistenceDiagramExtraction
# Objective : Extract sub level or power filtration persistence diagrams from graphs
# requires dataset folders from http://networkrepository.com/labeled.php
# Created by: Cuneyt Akcora
# Created on: 2020-12-18
library(igraph)
library(TDA)
rm(list = ls())
options(java.parameters = "-Xmx200g")


# 3 Graph ML datasets
dataPath <- "C:/Users/akkar/Documents/GraphML/"

p1 <-  "ENZYMES/ENZYMES."
p2 <-  "proteins/proteins."
p3 <-  "REDDIT-BINARY/REDDIT-BINARY."
p4 <-  "IMDB-MULTI/IMDB-MULTI."
p5 <-  "NCI1/NCI1."
p6 <-  "IMDB-BINARY/IMDB-BINARY."

# maximum simplex dimension
maxDimension <- 3
# node count thresholds to ignore a graph
minNodeCount <- 4

nodeFeatures <- c("betweenness","closeness","authority","degree")

subresultsFile <- "PDSubFiltration.txt"
powresultsFile <- "PDPowFiltration.txt"

useSubLevel <- TRUE
# if subLevel=FALSE, consider reducing maxNodeCount because power filtration
# does not work with large graphs
maxNodeCount <- 6000

computeNodeVals <- function(feature, thisGraph) {
  if (feature == "degree") {
    nodeValues <- apply(as_adjacency_matrix(thisGraph), 1, sum)
    # normalize values. This remains unused because we had better results without normalization
    #nodeValues=nodeValues/max(nodeValues)
  }else if (feature == "authority") {
    nodeValues = authority_score(thisGraph)$vector
    nodeValues<-nodeValues/max(nnodeValues)
     
  }else if (feature == "closeness") {
    nodeValues = closeness(thisGraph,normalized=TRUE)
     
  }else if (feature == "betweenness") {
    nodeValues = betweenness(thisGraph,normalized=TRUE)
    
  }else {
    message(nodeFeature," has not ben implemented as a node function in computeNodeVals()")
    }
  return(nodeValues)
}


for(nodeFeature in nodeFeatures){
  whichOutputFile <- if (useSubLevel){subresultsFile}else{powresultsFile}
  whichOutputFile<-paste0(nodeFeature,whichOutputFile)
  #Check its existence
  if (file.exists(whichOutputFile)) {
    file.remove(whichOutputFile)
  }
  
  for (path in c(p1, p2, p3, p4, p5, p6)) {
    edgeFile <- paste0(dataPath, path, "edges")
    graphFile <- paste0(dataPath, path, "graph_idx")
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
    message("Processing ",path," for feature ",nodeFeature)
    for (graphId in unique(graphIdData$V1)) {
      thisGraphNodes <- which(graphIdData$V1 == graphId)
      edgeData <-
        allData[allData$source %in% thisGraphNodes |
                  allData$target %in% thisGraphNodes, ]
      graph <- graph.data.frame(edgeData, directed = FALSE)
      originalV <- vcount(graph)
      # if the graph is not too small nor too big, compute filtration
      
      if (originalV > minNodeCount && originalV < maxNodeCount) {
        # below we use power filtration or sub/super level filtrations.
        # You need to use one of them to created the filtration object
        if (useSubLevel) {
          #1- sublevel filtration
          F.values = computeNodeVals(nodeFeature, graph)
          
          # for maxDimension=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
          cmplx <-
            cliques(as.undirected(graph), min = 0, max = maxDimension)
          # use sublevel=T for sublevel, sublevel=F for superlevel filtration
          FltRips <-
            funFiltration(FUNvalues = F.values,
                          cmplx = cmplx,
                          sublevel = T) # Construct filtration using F.values
        } else{
          # 2 - power filtration
          distanceMatrix <-
            shortest.paths(graph, v = V(graph), to = V(graph))
          maxScale <- max(distanceMatrix)
          FltRips <-
            ripsFiltration(
              distanceMatrix,
              maxdimension = maxDimension,
              maxscale = maxScale,
              dist = 'arbitrary',
              printProgress = FALSE
            )
        }
        #extract the persistence diagram
        persistenceDiagram <-
          filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
        
        #extract betti number
        b0OrigGraph <- sum(persistenceDiagram[, 1] == 0)
        b1OrigGraph <- sum(persistenceDiagram[, 1] == 1)
        b2OrigGraph <- sum(persistenceDiagram[, 1] == 2)
        
        #extract birth and death times
        pd2 <-
          cbind(persistenceDiagram,
                GraphId = graphId,
                dataset = path)
        
        write.table(
          pd2,
          file = whichOutputFile,
          sep = "\t",
          row.names = FALSE,
          col.names = FALSE,
          append = T,
          quote = FALSE
        )
        
      } else {
        message("Ignoring ",
                path,
                " graph ",
                graphId,
                " Node count:",
                originalV)
      }
    }
  }
}