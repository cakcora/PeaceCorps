# Title     : PeaceCorps
# Objective : Extract sub level or power filtration persistence diagrams from graphs
# requires dataset folders from http://networkrepository.com/labeled.php
# Created by: Cuneyt Akcora
# Created on: 2020-12-18
library(igraph)
library(TDA)
rm(list = ls())
options(java.parameters = "-Xmx200g")

# 3 Graph ML datasets
p1 <- "GraphML/ENZYMES/ENZYMES."
p2 <- "GraphML/proteins/proteins."
p3 <- "GraphML/REDDIT-BINARY/REDDIT-BINARY."
p4 <- "GraphML/IMDB-MULTI/IMDB-MULTI."
p5 <- "GraphML/NCI1/NCI1."
p6 <- "GraphML/IMDB-BINARY/IMDB-BINARY."

# maximum simplex dimension
maxDimension <- 3
# node count thresholds to ignore a graph
maxNodeCount <- 800
minNodeCount <- 4

subresultsFile <- "graphMlResultsSubFiltration.txt"
powresultsFile <- "graphMlResultsPowFiltration.txt"

useSubLevel<-TRUE

#Check its existence
if(useSubLevel<-TRUE){
  if (file.exists(subresultsFile)) {
    file.remove(subresultsFile)
  }
}else if (file.exists(powresultsFile)) {
  file.remove(powresultsFile)
}


for (path in c(p1, p3, p2, p6, p4, p5)) {
  edgeFile <- paste0(path, "edges")
  graphFile <- paste0(path, "graph_idx")
  graphIdData <- read.table(graphFile, quote = "\"", comment.char = "", sep = ",")
  allData <- read.table(edgeFile, quote = "\"", comment.char = "", sep = ",")
  colnames(allData) <- c("source", "target")
  for (graphId in unique(graphIdData$V1)) {
    thisGraphNodes <- which(graphIdData$V1 == graphId)
    edgeData <- allData[allData$source %in% thisGraphNodes | allData$target %in% thisGraphNodes,]
    graph <- graph.data.frame(edgeData, directed = FALSE)
    originalV <- vcount(graph)
    # if the graph is not too small nor too big, compute filtration

    if (originalV > minNodeCount && originalV < maxNodeCount) {
      # below we use power filtration or sub/super level filtrations.
      # You need to use one of them to created the filtration object
      if(useSubLevel){
        #1- sublevel filtration
        F.values = apply(as_adjacency_matrix(graph), 1, sum) # F.value_i=amount sent by node i
        #F.values=F.values/max(F.values) # normalize F.values
        # for maxDimension=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
        cmplx <- cliques(as.undirected(graph), min = 0, max = maxDimension)
        FltRips <- funFiltration(FUNvalues = F.values, cmplx = cmplx, sublevel = T) # Construct filtration using F.values
      }else{
        # 2 - power filtration
        distanceMatrix <- shortest.paths(graph, v = V(graph), to = V(graph))
        maxScale <- max(distanceMatrix)
        FltRips<-ripsFiltration(distanceMatrix, maxdimension = maxDimension, maxscale = maxScale, dist = 'arbitrary',printProgress = FALSE)
      }
      #extract the persistence diagram
      persistenceDiagram <- filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram

      #extract betti number
      b0OrigGraph <- sum(persistenceDiagram[, 1] == 0)
      b1OrigGraph <- sum(persistenceDiagram[, 1] == 1)
      b2OrigGraph <- sum(persistenceDiagram[, 1] == 2)

      #extract birth and death times
      pd2 <- cbind(persistenceDiagram, GraphId = graphId, dataset = path)
      whichFile<-if(useSubLevel) subresultsFile else powresultsFile
      write.table(pd2, file = whichFile, sep = "\t", row.names = FALSE, col.names = FALSE, append = T, quote = FALSE)

      message(path, "\t", graphId, "\t", originalV, "\t",
              b0OrigGraph, "\t", b1OrigGraph, "\t")
    }else {
      message("Ignoring ", path, " graph ", graphId, " Node count:", originalV)
    }
  }
}
