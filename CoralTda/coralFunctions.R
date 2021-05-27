library(TDA)
library(igraph)
computeCoralTDA<-function(graph,bettiTarget,runTDA=TRUE){
  nodeCountStd=vcount(graph)
  edgeCountStd = ecount(graph)
  cmplxStd <- cliques(graph, min = 1, max = bettiTarget+2)
  cliqueSizeStd=length(cmplxStd)
  bettiValStd=0
  if(runTDA){
    Flt <- funFiltration(FUNvalues = degree(graph), cmplx = cmplxStd, sublevel = T)
    PD <- filtrationDiag(filtration = Flt, maxdimension = bettiTarget,library = 'Dionysus', location=T)$diagram
    vpe=PD[PD[, 1] == bettiTarget,]
    if(!inherits(vpe, "data.frame")) 
      bettiValStd =length(vpe)
    else
      bettiValStd=nrow(PD[PD[, 1] == bettiTarget,])
  }
  
  coreness <- graph.coreness(graph)
  Fvalues2 = degree(graph)
  coralGraph<- induced.subgraph(graph=graph,vids=names(which(coreness>bettiTarget)))
  nodeCountCoral=vcount(coralGraph)
  edgeCountCoral = ecount(coralGraph)
  cliqueSizeCoral=0
  bettiValCoral=0
  if(vcount(coralGraph)>0){
    
    cmplxCoral <- cliques(coralGraph, 
                          min = 1, max = bettiTarget+2 )
    cliqueSizeCoral=length(cmplxCoral)
    bettiValCoral=0
    if(runTDA){
      Fvalues = (Fvalues2[intersect(names(Fvalues2),names(V(coralGraph)))])
      Flt2 <- funFiltration(FUNvalues = Fvalues, cmplx = cmplxCoral, sublevel = T)
      if(vcount(coralGraph)!=length(Fvalues))
        message("error: ",vcount(coralGraph)," is not equal to ",length(cmplxCoral))
      PD2 <- filtrationDiag(filtration = Flt2, maxdimension = bettiTarget,library = 'Dionysus', location=T)$diagram
      vpd=PD2[PD2[, 1] == bettiTarget,]
      if(!inherits(vpd, "data.frame")) {
        bettiValCoral =length(vpd)
        #print(PD2[PD2[, 1] == bettiTarget,])
      }
      else
        bettiValCoral=nrow(PD2[PD2[, 1] == bettiTarget,])
      #print(PD2[PD2[, 1] == bettiTarget,])
    }
  }else{
    #message("Betti ",bettiTarget," graph does not have enough edges")
  }  
  return(list(cliqueSizeStd,nodeCountStd,edgeCountStd,
              cliqueSizeCoral,nodeCountCoral,edgeCountCoral,bettiValStd,bettiValCoral))
}

computeTemporalTDA<-function(graph,bettiTarget,MCRuns=10){
  nodeCountStd=vcount(graph)
  edgeCountStd = ecount(graph)
  cmplxStd <- cliques(graph, min = 1, max = bettiTarget+2)
  cliqueSizeStd=length(cmplxStd)
  bettiValStd=0
  runVector<-rep(0,MCRuns)
  for(i in seq(1,MCRuns,1)){
    tstd<-system.time({
      Flt <- funFiltration(FUNvalues = degree(graph), cmplx = cmplxStd, sublevel = T)
      PD <- filtrationDiag(filtration = Flt, maxdimension = bettiTarget,library = 'Dionysus', location=T)$diagram
    })[[1]]
    runVector[[i]]=tstd
  }
  
  
  vpe=PD[PD[, 1] == bettiTarget,]
  if(!inherits(vpe, "data.frame")) 
    bettiValStd =length(vpe)
  else
    bettiValStd=nrow(PD[PD[, 1] == bettiTarget,])
  
  
  coreness <- graph.coreness(graph)
  Fvalues2 = degree(graph)
  coralGraph<- induced.subgraph(graph=graph,vids=names(which(coreness>bettiTarget)))
  nodeCountCoral=vcount(coralGraph)
  edgeCountCoral = ecount(coralGraph)
  cliqueSizeCoral=0
  bettiValCoral=0
  tcoral=0
  if(vcount(coralGraph)>0){
    
    cmplxCoral <- cliques(coralGraph, 
                          min = 1, max = bettiTarget+2 )
    cliqueSizeCoral=length(cmplxCoral)
    bettiValCoral=0
    
    Fvalues = (Fvalues2[intersect(names(Fvalues2),names(V(coralGraph)))])
    if(vcount(coralGraph)!=length(Fvalues))
      message("error: ",vcount(coralGraph)," is not equal to ",length(cmplxCoral))
      runCoralVector = rep(0,MCRuns)
      for(i in seq(1,MCRuns,1)){
        tcoral<-system.time({
        Flt2 <- funFiltration(FUNvalues = Fvalues, cmplx = cmplxCoral, sublevel = T)
        PD2 <- filtrationDiag(filtration = Flt2, maxdimension = bettiTarget,library = 'Dionysus', location=T)$diagram
      })[[1]]
        runCoralVector[[i]]=tcoral
        
    }
    #message(tstd," ",tcoral," seconds")
    vpd=PD2[PD2[, 1] == bettiTarget,]
    if(!inherits(vpd, "data.frame")) {
      bettiValCoral =length(vpd)
    }
    else
      bettiValCoral=nrow(PD2[PD2[, 1] == bettiTarget,])
    
  }   
  return(list(cliqueSizeStd,nodeCountStd,edgeCountStd,
              cliqueSizeCoral,nodeCountCoral,edgeCountCoral,bettiValStd,bettiValCoral,median(tstd),median(tcoral)))
}
