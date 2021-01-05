# Title     : GraphSignatureExtraction
# Objective : Extract graph betti signatures from persistence diagrams
# requires results from GraphPersistenceDiagramExtraction.R
# Created by: Cuneyt Akcora
# Created on: 2020-12-18
rm(list = ls())
library(dplyr)
library(plyr)
library(igraph)
useSubLevel <- TRUE
subresultsFile <- "graphMlResultsSubFiltration.txt"
powresultsFile <- "graphMlResultsPowFiltration.txt"
whichResulFile <- if (useSubLevel)  subresultsFile else powresultsFile
data <- read.table(file = whichResulFile, header = F, sep = "\t")
# the first line contains NA due to creating a data frame a
data <- na.omit(data)
colnames(data) <- c("betti", "birth", "death", "graphId", "dataset")
# does every graph have noth betti 0 and 1?
sanityCheck <- ddply(data, .(dataset, graphId), summarize, l = length(unique(betti)))
table(sanityCheck$l)

subSignatureFile <- "graphMlResultsAnalysisSub.txt"
powSignatureFile <- "graphMlResultsAnalysisPow.txt"

whichFile <- if (useSubLevel) subSignatureFile else powSignatureFile
if (file.exists(whichFile)) {
  #Delete file if it exists
  file.remove(whichFile)
}

# format persistent shapes (end at infinity)
data$death[!is.finite(data$death)] <- NA


for (datasetName in unique(data$dataset)) {
  message("Processing ", datasetName)
  dataOfADataset <- data[data$dataset == datasetName,]

  for (graphId in unique(dataOfADataset$graphId)) {
    specificGraphData_ <- dataOfADataset[dataOfADataset$graphId == graphId,]
    for (bettiNum in unique(specificGraphData_$betti)) {
      specificGraphData = specificGraphData_[specificGraphData_$betti == bettiNum,]
      starts <- unique(specificGraphData$birth)
      # count how many births occur at what threshold
      birthData <- ddply(specificGraphData, .(birth), summarise, bs = length(birth))
      # count how many deaths occur at what threshold
      deathData <- ddply(specificGraphData, .(death), summarise, bs = length(birth))
      deaths <- unique(deathData$death)
      nrowResult <- nrow(birthData)

      thresholdValues <- sort(unique(c(starts, deaths)))
      thresholdValues <- thresholdValues[!is.na(thresholdValues)]
      bettiSignatureArray = array(c(0, 0))

      # create the betti signature.
      # we use an index scheme of 3* threshold
      # index 3*threhold-1 holds bettis that exist before deaths at the threshold
      # index 3*threhold holds bettis that remain after deaths at the threshold
      # index 3*threshold+1 holds bettis that remain after deaths + those that are newly born at the threshold
      if (nrowResult > 0) {
        latestCountValue <- 0
        for (threshold in thresholdValues) {
          birthCount <- birthData[birthData$birth == threshold,]$bs[1]
          deathCount <- deathData[deathData$death == threshold,]$bs[1]
          bettiSignatureArray[threshold * 3] <- latestCountValue
          if (!is.null(deathCount) & !is.na(deathCount)) {
            bettiSignatureArray[threshold * 3] = bettiSignatureArray[threshold * 3] - deathCount
          }
          if (!is.null(birthCount) & !is.na(birthCount)) {
            bettiSignatureArray[threshold * 3 + 1] = bettiSignatureArray[threshold * 3] + birthCount
          }else bettiSignatureArray[threshold * 3 + 1] = bettiSignatureArray[threshold * 3]
          latestCountValue = bettiSignatureArray[threshold * 3 + 1]
        }
      }
      # in our 3* thresholds scheme, some values will remain NA because
      # no new bettis are born or die around some thresholds.
      for (i in seq(1, length(bettiSignatureArray))) {
        if (is.na(bettiSignatureArray[i])) {
          bettiSignatureArray[i] <- bettiSignatureArray[i - 1]
        }
      }
      # if all bettis die, the signature will end in 0
      # otherwise, if some bettis survive, the signature will end in a value>0

      # write signatures to the file
      bettiSignatureArray = cbind(datasetName, graphId, bettiNum, paste(bettiSignatureArray, collapse = " "))
      write.table(bettiSignatureArray, file = subSignatureFile, sep = "\t", row.names = FALSE, col.names = FALSE, append = T, quote = FALSE)
    }
  }
}


