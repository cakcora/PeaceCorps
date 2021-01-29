# Title     : PDSignatureExtractionForDiscreteFeatures.R
# Objective : Extract graph betti signatures from persistence diagrams
# requires results from GraphPersistenceDiagramExtraction.R
# Created by: Cuneyt Akcora
# Created on: 2020-12-18
rm(list = ls())
library(dplyr)
library(plyr)
library(igraph)


useSubLevel <- TRUE
subresultsFile <- "PDSubFiltration.txt"
powresultsFile <- "PDPowFiltration.txt"
dataInputPath="C://Users/akkar/IdeaProjects/PeaceCorps/"
dataOutputPath="C://Users/akkar/IdeaProjects/PeaceCorps/"
nodeFeatures <- c("degree","eccentricity")
subSignatureFile <- "SignatureSub.txt"
powSignatureFile <- "SignaturePow.txt"

sLen=100
signatureLength<-3*sLen

for(nodeFeature in nodeFeatures){
  whichInputFile <- if (useSubLevel)  subresultsFile else powresultsFile
  whichInputFile<-paste0(dataInputPath
                        ,nodeFeature,whichInputFile)
  
  data <- read.table(file = whichInputFile, header = F, sep = "\t")
  # the first line contains NA due to creating a data frame
  data <- na.omit(data)
  colnames(data) <- c("betti", "birth", "death", "graphId", "dataset")
   
   
  # does every graph have noth betti 0 and 1?
  sanityCheck <- ddply(data, .(dataset, graphId), summarize, l = length(unique(betti)))
  table(sanityCheck$l)
  
  whichOutputFile <- if (useSubLevel) subSignatureFile else powSignatureFile
  whichOutputFile<-paste0(dataOutputPath,nodeFeature,whichOutputFile)
  if (file.exists(whichOutputFile)) {
    #Delete file if it exists
    file.remove(whichOutputFile)
  }
  
  # format persistent shapes (that end at infinity)
  data$death[!is.finite(data$death)] <- NA
  
  for (datasetName in unique(data$dataset)) {
    message("Processing ", datasetName, " for feature ",nodeFeature)
    dataOfADataset <- data[data$dataset == datasetName,]
    
    for (graphId in unique(dataOfADataset$graphId)) {
      specificGraphData_ <- dataOfADataset[dataOfADataset$graphId == graphId,]
      for (bettiNum in unique(specificGraphData_$betti)) {
        specificGraphData = specificGraphData_[specificGraphData_$betti == bettiNum,]
        # Normalize birth and death threshod values to an interval between 0 and sLen
        maxThres <- max(specificGraphData[,c("birth","death")],na.rm=T)
        specificGraphData$birth<-as.integer(sLen*specificGraphData$birth/maxThres)
        specificGraphData$death<-as.integer(sLen*specificGraphData$death/maxThres)
        if(max(specificGraphData$death,na.rm = TRUE)>sLen||min(specificGraphData$birth,na.rm = TRUE)<0){
          message("Warning: ", nodeFeature," values for ",datasetName," is not within [0,sLen]")
        }
        # normalization ends here
        
        starts <- unique(specificGraphData$birth)
        # count how many births occur at what threshold
        birthData <- ddply(specificGraphData, .(birth), summarise, bs = length(birth))
        # count how many deaths occur at what threshold
        deathData <- ddply(specificGraphData, .(death), summarise, bs = length(birth))
        deaths <- unique(deathData$death)
        nrowResult <- nrow(birthData)
        
        thresholdValues <- sort(unique(c(starts, deaths)))
        thresholdValues <- thresholdValues[!is.na(thresholdValues)]
        bettiSignatureArray = array(c(0))
        
        # create the betti step signature.
        # we use an index scheme of 3* threshold
        # index 3*threhold-1 holds bettis that exist before deaths at the threshold
        # index 3*threhold holds bettis that remain after deaths at the threshold
        # index 3*threshold+1 holds bettis that remain after deaths + those that are newly born at the threshold
        if (nrowResult > 0) {
          latestCountValue <- 0
          minThreshold=min(starts,na.rm = TRUE)
          for (threshold in thresholdValues) {
            birthCount <- birthData[birthData$birth == threshold,]$bs[1]
            deathCount <- deathData[deathData$death == threshold,]$bs[1]
            if(is.na(deathCount)) deathCount=0
            if(is.na(birthCount)) birthCount=0
            if(threshold==minThreshold){
              bettiSignatureArray[3*minThreshold+1] <- birthCount-deathCount
              bettiSignatureArray[3*minThreshold] <- birthCount
            }else{
              bettiSignatureArray[(threshold * 3) - 1] <- latestCountValue
              bettiSignatureArray[threshold * 3] <- latestCountValue - deathCount
              bettiSignatureArray[threshold * 3 + 1] <- bettiSignatureArray[threshold * 3] + birthCount
            }
            latestCountValue = bettiSignatureArray[threshold * 3 + 1]
          }
          # in our 3* thresholds scheme, some values will remain NA because
          # no new bettis are born or die around some thresholds.
          # We will replace these NA values with values stored in earlier cells.
          # Also,
          # 1 - If all bettis die, the signature will end in 0
          # 2 - If some bettis survive, the signature will end in a non-zero value
          # For example, a signature that is 0 NA 1 8 NA NA NA 8 7 11 NA NA
          # will become 0 0 1 8 8 8 8 8 7 11 11 11
          for (i in seq(1, (1+signatureLength))) {
            if (is.na(bettiSignatureArray[i])) {
              if(i==1){
                bettiSignatureArray[i]=0
              }else {
                bettiSignatureArray[i] <- bettiSignatureArray[i - 1]
              }
            } 
          }
          
          
          # write signatures to the file
          bettiSignatureArray = cbind(datasetName, graphId, bettiNum, paste(bettiSignatureArray, collapse = " "))
          write.table(bettiSignatureArray, file = whichOutputFile, sep = "\t", row.names = FALSE, col.names = FALSE, append = T, quote = FALSE)
        }
      }
    }
  }
  }
  