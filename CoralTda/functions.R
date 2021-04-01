computeStdTda<-function(dat,verbose,bettiTarget){
  allData<-dat
  allData$source<-paste0("V",allData$source)
  allData$target<-paste0("V",allData$target)
  graph <- graph.data.frame(allData, directed = FALSE)
  #graph<-simplify(graph)
  nodeCount=vcount(graph)
  edgeCount = ecount(graph)
  if(verbose){
    plot(graph)
    message("Graph has ",vcount(graph)," nodes and ",ecount(graph)," edges")
  }
  #if(vcount(graph)==0) return(c())
  #timer<-system.time({})
  #1 ordinary sublevel filtration
 
  tryCatch({
    # F.values are node values. At these values, node appear in the complex,
    # and their edges to other active nodes also activate.
    F.values = degree(graph, v = V(graph), mode = c("all"));
  }, error = function(error_condition) {
    message("Error in finding node activation values for  ",dataAlias)
    return(-1);
  })
  if(verbose){
    print("computed node values")
    print(F.values)
    }
  
  nodeCountinClique<-bettiTarget+2
  cmplx <- cliques(graph, min = 1, max = nodeCountinClique)
  cliqueSize=length(cmplx)
  FltRips <- funFiltration(FUNvalues =F.values, cmplx = cmplx,sublevel = T) 
  if(verbose)print("extracted the filtration")
  #extract the persistence diagram
  tryCatch({
    persistenceDiagram <-
      filtrationDiag(filtration = FltRips, library="dionysus", maxdimension = bettiTarget)$diagram
    if(verbose){
      print("Standard TDA values:")
      print((persistenceDiagram))
    }
  }, error = function(error_condition) {
    #print(error_condition)
  })
  
  #extract birth and death times
  res=data.frame(persistenceDiagram[persistenceDiagram[, 1] == bettiTarget,])
  return(list(cliqueSize,res,edgeCount,nodeCount))
}

computeCoralTda<-function(allData,verbose=F,bettiTarget){
  allData$source<-paste("V",allData$source)
  allData$target<-paste("V",allData$target)
  graph <- graph.data.frame(allData, directed = FALSE)
  
  graph<-simplify(graph)
  coreness <- graph.coreness(graph)
  #message("\t","decomposition ",timer[[1]]," ",vcount(cGraph)," ",ecount(cGraph))
  #F.values =  apply(as_adjacency_matrix(cGraph), 1, sum)
  F2.values = degree(graph, v = V(graph), mode = c("all"))
  graph <- induced.subgraph(graph=graph,vids=names(which(coreness>bettiTarget)))
  cliqueSize=0
  nodeCount=0
  edgeCount=0
  if(vcount(graph)>0){
    if(verbose)plot(graph)
    F.values = (F2.values[intersect(names(F2.values),names(V(graph)))])
    nodesinClique=bettiTarget+2
    cmplx <- cliques(graph, min = 0, max = nodesinClique)
    cmplx2 <- cliques(graph, min = bettiTarget, max = nodesinClique)
    cliqueSize=length(cmplx2)
    nodeCount=vcount(graph)
    edgeCount = ecount(graph)
    if(verbose)
      message(bettiTarget," number of cliques are ",length(cmplx))
    
    FltRips <- funFiltration(FUNvalues = F.values,
                             cmplx = cmplx,  
                             sublevel = T) 
    #extract the persistence diagram
    persistenceDiagram <-
      filtrationDiag(filtration = FltRips, library="gudhi",maxdimension = bettiTarget)$diagram
    if(verbose){
      print("Coral TDA values:")
      print(persistenceDiagram)
    }
    bCore <- sum(persistenceDiagram[, 1] == bettiTarget)
    if(verbose&bettiTarget>1 & bCore>0)
      message(bettiTarget," ---------------------",bCore)
    res<-data.frame(persistenceDiagram[persistenceDiagram[, 1] == bettiTarget,])
     
  }else{
    if(verbose)message("graph does not have enough edges")
    res<-data.frame()
  }  
  return(list(cliqueSize,res,edgeCount,nodeCount))
}

