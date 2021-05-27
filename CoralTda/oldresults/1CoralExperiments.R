library(igraph)
library(TDA)
library(dplyr)
rm(list = ls())
options(java.parameters = "-Xmx200g")


#  Graph ML datasets
dataPath <- "C:/data/"
dataPath <- "/home/jupiter/GraphML/"
feature = "degree"
# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
#maxDimension <- 2

# upper limit on the simplex dimension, size of the cliques to find
maxClique <-5

#betti 2 icin max simple size 3 veya ustunde olmasi bir sey farkettirmeyecek.
#betti 3 icin 4 ve uzeri olmasi farketmez.

testRange<-c(1,2,3,4)

dataset="NCI1/NCI1."
dataAlias="NCI1"
source("nodeFunctions.R")
loadGraph<-function(dataset="NCI1/NCI1.", dataAlias="NCI1"){
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
  times<-c()
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
    nGraph<-as.undirected(graph)
    timer<-system.time({
      for(maxDimension in testRange){
        
        #sublevel filtration
        # compute node function values 
        tryCatch({
          F.values = degree(nGraph, v = V(nGraph), mode = c("all"));
        }, error = function(error_condition) {
          message("Error in extracting PDs for ",dataAlias)
          return(-1);
        })
        # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
       
        cmplx <- cliques(nGraph, min = 0, max = maxClique)
        # use sublevel=T for sublevel, sublevel=F for superlevel filtration
        # F.values are node values. At these values, node appear in the complex,
        # and their edges to other active nodes also activate.
        FltRips <- funFiltration(FUNvalues = F.values,
                                 cmplx = cmplx,
                                 sublevel = T) # Construct filtration using F.values
        
        #extract the persistence diagram
        persistenceDiagram <-
          filtrationDiag(filtration = FltRips, library="Dionysus", maxdimension = maxDimension)$diagram
        b0OrigGraph <- sum(persistenceDiagram[, 1] == maxDimension)
        
        #extract birth and death times
        pd0 <-data.frame(count=b0OrigGraph,betti=maxDimension,GraphId = graphId,dataset = dataAlias,style="all")
        results<-rbind(results,pd0)
      }
    })
    
    timeAll<-(timer[[1]]);
    ### coral computations
    
    cGraph<-nGraph
     timer<-system.time({
    coreness <- graph.coreness(cGraph)
    for(maxDimension in testRange){
        cGraph <- induced.subgraph(graph=cGraph,vids=names(which(coreness>maxDimension)))
        #message(maxDimension, " ",vcount(cGraph))
        if(vcount(cGraph)>0){
          tryCatch({
            #F.values =  apply(as_adjacency_matrix(cGraph), 1, sum)
            F.values = degree(cGraph, v = V(cGraph), mode = c("all"))
          }, error = function(error_condition) {
            message("Error in extracting PDs for ",dataAlias)
            message(error_condition)
            return(-1);
          })
          
          cmplx <- cliques(cGraph, min = 0, max = maxClique)
          # use sublevel=T for sublevel, sublevel=F for superlevel filtration
          # F.values are node values. At these values, node appear in the complex,
          # and their edges to other active nodes also activate.
          FltRips <- funFiltration(FUNvalues = F.values,
                                   cmplx = cmplx,
                                   sublevel = T) # Construct filtration using F.values
          #extract the persistence diagram
          persistenceDiagram <-
            filtrationDiag(filtration = FltRips, library="Dionysus",maxdimension = maxDimension)$diagram
          bCore <- sum(persistenceDiagram[, 1] == maxDimension)
          #if(maxDimension>1 & bCore>0)
            #message(maxDimension," ---------------------",bCore)
          pdCore <-data.frame(count=bCore,betti=maxDimension,GraphId = graphId,dataset = dataAlias,style="core")
          results<-rbind(results,pdCore)
        }
      
    }
    })
    timeCore<-(timer[[1]]);
    message(graphId," ", timeAll," ", timeCore/timeAll, " ",vcount(cGraph)/vcount(nGraph)," ",
            ecount(cGraph)/ecount(nGraph))
    times<-c(times,timeCore/timeAll)
  }
  message(dataAlias," ",mean(times,rm.na=T))
}


loadGraph("NCI1/NCI1.","NCI1")
loadGraph("proteins/proteins.","Protein")
loadGraph("REDDIT-BINARY/REDDIT-BINARY.","RedditBinary")
loadGraph("IMDB-MULTI/IMDB-MULTI.","IMDBMulti")
loadGraph("IMDB-BINARY/IMDB-BINARY.","IMDBBinary")