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

# Given an adjacency matrix, find how many vertices are dominated by each vertex
# we return a vector of domination counts
findDominatingVertex<-function(adjMatrix){
  len1=ncol(adjMatrix)
  colnames(adjMatrix)<-NULL
  rownames(adjMatrix)<-NULL
  for(k in 1:(len1-1)){
    vec1=adjMatrix[k,]
    
    for(j in (k+1):(len1)){
      if(adjMatrix[k,j]==1){#only run over neighbors
        vec2=adjMatrix[j,]
        if(sum(vec1>=vec2)==(len1-1)){
          return(k)
        }else if(sum(vec1<=vec2)==(len1-1)){
          return(j)
        }else if(sum(vec1==vec2)==(len1-1)){
          return(j)
        }
      }
    }
  }
  # no dominating vertex was found
  return(-1)
}

# given an adjacency matrix and a dominating vertex,
# delete vertices that are dominated.
deleteDominatedBy<-function(adjMatrix,dominator){
  vec1=adjMatrix[dominator,]
  len1=length(vec1)
  deleteThese=c()
  for(j in 1:len1){
    if(j!=dominator&&adjMatrix[j,dominator]==1){
      vec2=adjMatrix[j,]
      sub1_2=FALSE
      sub2_1=FALSE
      for(i in seq(1,len1)){
        if(i!=dominator&&i!=j){
          if(vec1[[i]]>vec2[[i]]){
            sub1_2=T
          }
          else if(vec1[[i]]<vec2[[i]]){
            sub2_1=T
          }
        }
      }
      
      if(sub1_2==T&&sub2_1==F)
        deleteThese<-c(deleteThese,j)
      
      else if(sub2_1==F&&sub1_2==F) {
        deleteThese<-c(deleteThese,j)
      }
      
    }
  }
  
  deleteThese<-deleteThese[deleteThese!=dominator]
  if(length(deleteThese)==0){
    # this should never happen
    message("error: dominating vertex leads to no vertex deletion.")
    return(NULL)
  }
  # delete both rows and columns of the selected vertices
  return(adjMatrix[-deleteThese,-deleteThese])
}


# the main algorithm that iterates over an adjacency matrix, and deletes 
# dominated vertices. 
dominationRatio<-function(adjMatrix){
  message("graph order:",nrow(adjMatrix))
  if(nrow(adjMatrix)<2) return(0)
  if(!isSymmetric(adjMatrix)){
    message("Error: The adj matrix is not symmetric.")
    return(-1)
  }
  if(length(which(adjMatrix>1, arr.ind=TRUE)>0)){
    message("Error: The adj matrix should contain 0s and 1s only.")
    return(-1)
  }
  tempMatrix=adjMatrix
  startvertex=1
  repeat{
    dominatingVertex=findDominatingVertex(tempMatrix)
    if(dominatingVertex==-1)break;
    tempMatrix = deleteDominatedBy(tempMatrix,dominatingVertex)
    # these break conditions are complex but required.
    if(sum(is.na(tempMatrix))>0)break;
    if(length(tempMatrix)<=1)break;
    if(dim(tempMatrix)<3)break;
  }
  size2=length(tempMatrix)
  if(size2>1) size2=nrow(tempMatrix)
  # reduction rate: 0.8 means 20% of the vertices have been removed.
  return((nrow(adjMatrix)-size2)/nrow(adjMatrix))
}

 
