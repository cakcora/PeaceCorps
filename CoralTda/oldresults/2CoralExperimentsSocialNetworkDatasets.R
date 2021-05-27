rm(list = ls())
library(igraph)
library(TDA)
library(dplyr)

options(java.parameters = "-Xmx200g")

#betti 2 icin clique min=3, max= 4 node olmali 
#betti 2yi k=3 core graph uzerinde hesapliyoruz
source("CoralTda/functions.R")
computeActivation<-function(cmplx,F.values){
  
  edgeFile <- paste0(dataPath, "twitter/5.", "edges")
  allData <- read.table(edgeFile, quote = "", comment.char = "", sep = " ")
  bettiTarget=2
  colnames(allData) <- c("source", "target")
  gr <- graph.data.frame(allData, directed = FALSE)
  plot(gr)
  F.values = degree( gr, v = V( gr), mode = c("all"));
  grcmplx <- cliques(gr, min = 1, max =bettiTarget+2)
  message("max clique size is ",length(grcmplx))
  sublvlfiltration<-funFiltration(FUNvalues =F.values , cmplx = grcmplx)
  filtrationDiag(filtration = sublvlfiltration, location=T,
                 library="Dionysus",maxdimension = bettiTarget)$diagram
  # 6 edgeli graphta bettitarget+1 ve bettitarget+2 ayni sonuclari vermeyecek
  
}

source("CoralTda/coralFunctions.R")

#  Graph ML datasets
dataPath <- "/home/jupiter/GraphML/"
dataAlias="Twitter"


loadGraph<-function(edgeFile){
  allData <- read.table(edgeFile, quote = "", comment.char = "", sep = " ")
  colnames(allData) <- c("source", "target")
  allData$source<-paste0("V",allData$source)
  allData$target<-paste0("V",allData$target)
  graph <- graph.data.frame(allData, directed = FALSE)
  return(graph)
}

dataPathTwitter=paste0(dataPath,"twitter/")
filesTwitter=list.files(path = dataPathTwitter,recursive=T, pattern = ".edges", all.files = TRUE,
                        full.names = TRUE)

dataPathFB=paste0(dataPath,"facebook/")
filesFB=list.files(path = dataPathFB,recursive=T, pattern = ".edges", all.files = TRUE,
                   full.names = TRUE)
outputtwitterFile = "CoralTda/coralTDAtwitterresults.csv"#for twitter
outputfbFile = "CoralTda/coralTDAfacebookresults.csv"#for twitter
outputFile=outputfbFile
#Check if a previous result file  exists
if (file.exists(outputFile)) {
  file.remove(outputFile)
}

toys<-c(1,2,3,4,5)
for(fileName in c(filesFB)){
  graph<-loadGraph(fileName)
  runTDA=TRUE
  for(bettiTarget in seq(1,10,1)){
    graph2=graph
    
    if(ecount(graph2)<7000){
      
      results=computeCoralTDA(graph2,bettiTarget,runTDA)
      cliqStd = results[[1]]
      nodeStd = results[[2]]
      edgeStd = results[[3]]
      
      cliqCoral = results[[4]]
      nodeCoral = results[[5]]
      edgeCoral = results[[6]]
      
      bettiStd = results[[7]]
      bettiCoral = results[[8]]
      fname = substr(fileName,start=1+nchar(dataPath),stop=nchar(fileName)-5)
      str=paste0(fname,"\t",bettiTarget, 
                       "\t",nodeStd,"\t", nodeCoral,"\t",edgeStd,"\t",edgeCoral,
                       "\t",cliqStd,
                       "\t",cliqCoral,
                       "\t",bettiStd,"\t",bettiCoral,"\r\n")
      cat(str,file=outputFile,append=TRUE)
      if(bettiStd>0) {
        message(str)
      }
      if(bettiStd==0&&bettiTarget>3){
        runTDA=FALSE
      }
    }
  }
}

