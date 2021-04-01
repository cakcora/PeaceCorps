library(ggplot2)
resKernel <- read.delim("~/PeaceCorps/CoralTda/coralTDAkernelresults.csv", header=FALSE)
colnames(resKernel)<-c("dataset","graphId","betti","nodeStd","nodeCoral","edgeStd","edgeCoral","cliqueStd","cliqueCoral","bettiStd","bettiCoral")
ggplot(data=resKernel,aes(betti,nodeCoral,group=dataset))+geom_line()
