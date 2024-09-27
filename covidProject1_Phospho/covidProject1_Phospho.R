# Mfuzz

load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/newPhosphoAll.rda")

newPhosphoRows <- rownames(newPhosphoAll)

phosphoPVals <- newPhosphoAll[, 6:9]
phosphoPVals <- apply(X = phosphoPVals, MARGIN = 2, FUN = as.numeric)
rownames(phosphoPVals) <- newPhosphoRows

phosphoData <- newPhosphoAll[,10:24]
phosphoData <- apply(X = phosphoData, MARGIN = 2, FUN = as.numeric)
rownames(phosphoData) <- newPhosphoRows

load(file = "C:/Users/usuario/Downloads/fatPaper/fatima/rda/descriptionSARS.rda")

proteomeDescription <- description
rm(description)
addToDescription <- proteomeDescription[newPhosphoAll[,2],]
phosphoDescription <- cbind(newPhosphoAll[,1:2],
                            addToDescription[,-1],
                            newPhosphoAll[,3:5])
colnames(phosphoDescription)[4] <- "Gene.ID"
phosphoDescription[is.na(phosphoDescription)] <- ""
phosphoDescription[,5] <- as.logical(phosphoDescription[,5])

saveRowNames <- rownames(phosphoDescription)
tabl <- table(phosphoDescription[,2])
tabl <- as.numeric(tabl[phosphoDescription[,2]])

catchSites <- regexpr(pattern = "[0-9]*xPhospho",
                      text = phosphoDescription$Modifications)
nSites <- substr(x = phosphoDescription$Modifications,
                 start = catchSites,
                 stop = catchSites + attr(x = catchSites, which = "match.length") - 9)
nSites <- as.numeric(nSites)
aggDF <- data.frame("Protein" = phosphoDescription$Protein, "Phosphosites" = nSites)

aggResult <- aggregate(x = aggDF[,-1], by = list(aggDF$Protein), FUN = sum)
aggResult[which(aggResult[1,] == ""), 2] <- NA
namedPhosphosites <- aggResult[,2]
names(namedPhosphosites) <- aggResult[,1]
namedPhosphosites <- namedPhosphosites[phosphoDescription[,2]]

phosphoDescription <- cbind.data.frame(phosphoDescription[,1:2],
                                       "All.phosphopeptides" = NA,
                                       "Significant.phosphopeptides" = tabl,
                                       "All.phosphosites" = NA,
                                       "Significant.phosphosites" = namedPhosphosites,
                                       phosphoDescription[,3:ncol(phosphoDescription)])

phosphoDescription <- phosphoDescription[order(phosphoDescription[,"Protein"]),]

saveRows <- rownames(phosphoData)
aggSum <- matrix(data = c(rowSums(phosphoData[,1:3]),
                          rowSums(phosphoData[,4:6]),
                          rowSums(phosphoData[,7:9]),
                          rowSums(phosphoData[,10:12]),
                          rowSums(phosphoData[,13:15])),
                 ncol = 5)
rownames(aggSum) <- saveRows
colnames(aggSum) <- c("MOCK.0h",
                      "SARS.CoV2.3h",
                      "SARS.CoV2.6h",
                      "SARS.CoV2.9h",
                      "SARS.CoV2.16h")
aggSum <- log2(aggSum)

allPhosphoDataPath <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/allR1.xlsx"
allPhosphoDescPath <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/allR2.xlsx"

allRawPhosphoData <- read.xlsx2(file = allPhosphoDataPath, sheetIndex = 1)
descAllPhosphoData <- read.xlsx2(file = allPhosphoDescPath, sheetIndex = 1)

firstCondition <- grepl(x = descAllPhosphoData[,"Modifications"], pattern = "Phospho")
secondCondition <- nchar(descAllPhosphoData[,3]) > 1
thirdCondition <- !(grepl(x = descAllPhosphoData[,3], pattern = "; "))
allConditions <- firstCondition & secondCondition & thirdCondition

filteredDesc <- descAllPhosphoData[allConditions,]

allRawPhosphoData <- allRawPhosphoData[,-ncol(allRawPhosphoData)]
dataColNames <- colnames(allRawPhosphoData)
dataColNames <- gsub(x = dataColNames, fixed = "TRUE", pattern = "X", replacement = "")
dataColNames <- gsub(x = dataColNames, fixed = "TRUE", pattern = "..", replacement = ".")
colnames(allRawPhosphoData) <- dataColNames

filteredData <- allRawPhosphoData[allConditions,]

rownames(filteredData) <- sprintf("ID.%s", 1:nrow(filteredData))
rownames(filteredDesc) <- rownames(filteredData)

sigPhosphoPath <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoReduced.xlsx"

sigDesc <- read.xlsx2(file = sigPhosphoPath, sheetIndex = 1)
sigData <- read.xlsx2(file = sigPhosphoPath, sheetIndex = 2)
colnames(sigData) <- gsub(pattern = "X",
                          replacement = "",
                          fixed = TRUE,
                          x = colnames(sigData))
colnames(sigData) <- gsub(pattern = "..",
                          replacement = ".",
                          fixed = TRUE,
                          x = colnames(sigData))
sigData <- sigData[,-ncol(sigData)]

rownames(sigData) <- sprintf("sig.ID.%s", 1:nrow(sigData))
rownames(sigDesc) <- rownames(sigData)

splitList <- strsplit(x = sigDesc[,"Protein"], split = "; ", fixed = TRUE)
firstSigCondition <- lapply(X = splitList, FUN = length)
firstSigCondition <- as.numeric(firstSigCondition)
firstSigCondition <- firstSigCondition > 1
firstSigCondition <- !(firstSigCondition)
secondSigCondition <- nchar(sigDesc[,3]) > 1
allSigConditions <- firstSigCondition & secondSigCondition

sigData <- sigData[allSigConditions,]
sigDesc <- sigDesc[allSigConditions,]
savedRows <- rownames(sigData)
sigData <- as.matrix(sigData)
sigData <- apply(X = sigData, MARGIN = 2, FUN = as.numeric)
rownames(sigData) <- savedRows

aggSum <- matrix(data = c(rowSums(sigData[,1:3])/3,
                          rowSums(sigData[,4:6])/3,
                          rowSums(sigData[,7:9])/3,
                          rowSums(sigData[,10:12])/3,
                          rowSums(sigData[,13:15])/3),
                 nrow = nrow(sigData),
                 ncol = 5)
colnames(aggSum) <- c("MOCK.t0",
                      "SARS.CoV2.t3",
                      "SARS.CoV2.t6",
                      "SARS.CoV2.t9",
                      "SARS.CoV2.t16")
rownames(aggSum) <- rownames(sigData)
aggSum <- log2(aggSum)

phosphoFuzz <- mfuzzing(aggregatedResult = aggSum, chosenRange = 2:10)

resPhosphoFuzz <- clAndPlot(mfuzzRes = phosphoFuzz,
                            chosenC = 6,
                            chosenMF = c(3, 2),
                            chosenTimes = c(0, 3, 6, 9, 16),
                            xlab2 = "Time (h)")

sigPhosphoProteinsClusters <- list()
for(i in 1:max(resPhosphoFuzz$cluster)){
  inames <- names(resPhosphoFuzz$cluster)[resPhosphoFuzz$cluster == i]
  iprots <- sigDesc[inames, 3]
  clusterName <- paste("cluster", as.numeric(i), sep = "")
  sigPhosphoProteinsClusters[[clusterName]] <- iprots
}

################################################################################

# Phospho-processes

allProteinsInPhospho <- unique(filteredDesc[,3])

# "UniProtR" was a library upon which the user-defined function "bioProcesses2"
# was defined. It can be possible that this function can't be run due to the
# changes that happened in the library. The result was saved in a variable that
# can be loaded, even if the original function can't be executed.

refinedPhosphoMatrix <- bioProcesses2(uniProtAccessionVector = allProteinsInPhospho,
                                      pvalVec = rep(0.01, length(allProteinsInPhospho)))

# save(x = refinedPhosphoMatrix,
#      file = "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")

# load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")

allSigPhospho <- unique(unlist(sigPhosphoProteinsClusters))
phosphoAlteredProcesses <- newAlteredProcesses(refinedPhosphoMatrix, allSigPhospho)

c1.3 <- unique(c(sigPhosphoProteinsClusters$cluster1,
                 sigPhosphoProteinsClusters$cluster3))
c2.5 <- unique(c(sigPhosphoProteinsClusters$cluster2,
                 sigPhosphoProteinsClusters$cluster5))
c3.6 <- unique(c(sigPhosphoProteinsClusters$cluster3,
                 sigPhosphoProteinsClusters$cluster6))
c4.5 <- unique(c(sigPhosphoProteinsClusters$cluster4,
                 sigPhosphoProteinsClusters$cluster5))

phosphoAlteredC1 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster1)
phosphoAlteredC2 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster2)
phosphoAlteredC3 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster3)
phosphoAlteredC4 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster4)
phosphoAlteredC5 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster5)
phosphoAlteredC6 <- newAlteredProcesses(refinedPhosphoMatrix,
                                        sigPhosphoProteinsClusters$cluster6)

phosphoAlteredC13 <- newAlteredProcesses(refinedPhosphoMatrix,
                                         c1.3)
phosphoAlteredC25 <- newAlteredProcesses(refinedPhosphoMatrix,
                                         c2.5)
phosphoAlteredC45 <- newAlteredProcesses(refinedPhosphoMatrix,
                                         c4.5)
phosphoAlteredC36 <- newAlteredProcesses(refinedPhosphoMatrix,
                                         c3.6)

allPhosphoProcesses <- list("cluster1" = phosphoAlteredC1,
                            "cluster2" = phosphoAlteredC2,
                            "cluster3" = phosphoAlteredC3,
                            "cluster4" = phosphoAlteredC4,
                            "cluster5" = phosphoAlteredC5,
                            "cluster6" = phosphoAlteredC6,
                            "clusters13" = phosphoAlteredC13,
                            "clusters25" = phosphoAlteredC25,
                            "clusters36" = phosphoAlteredC36,
                            "clusters45" = phosphoAlteredC45
                            )

allPValsPhospho <- list("cluster1" = c(),
                        "cluster2" = c(),
                        "cluster3" = c(),
                        "cluster4" = c(),
                        "cluster5" = c(),
                        "cluster6" = c(),
                        "clusters13" = c(),
                        "clusters25" = c(),
                        "clusters36" = c(),
                        "clusters45" = c()
                        )

allLog2FCPhospho <- allPValsPhospho

relevantClusters <- list("cluster1" = sigPhosphoProteinsClusters$cluster1,
                         "cluster2" = sigPhosphoProteinsClusters$cluster2,
                         "cluster3" = sigPhosphoProteinsClusters$cluster3,
                         "cluster4" = sigPhosphoProteinsClusters$cluster4,
                         "cluster5" = sigPhosphoProteinsClusters$cluster5,
                         "cluster6" = sigPhosphoProteinsClusters$cluster6,
                         "clusters13" = c1.3,
                         "clusters25" = c2.5,
                         "clusters36" = c3.6,
                         "clusters45" = c4.5
                         )

for(i in 1:length(relevantClusters)){
  
  designPVal <- rep(1, length(allProteinsInPhospho))
  names(designPVal) <- allProteinsInPhospho
  designLog2FC <- rep(0, length(allProteinsInPhospho))
  names(designLog2FC) <- allProteinsInPhospho
  
  designPVal[relevantClusters[[i]]] <- 0.01
  designLog2FC[relevantClusters[[i]]] <- 1
  
  allPValsPhospho[[i]] <- designPVal
  allLog2FCPhospho[[i]] <- designLog2FC
  
}

names(allPValsPhospho) <- names(relevantClusters)
names(allLog2FCPhospho) <- names(relevantClusters)

load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")

# Individual clusters

phosphocluster1APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                            descript = completeDescription,
                                            chosenPValue = 0.05,
                                            pValuesVector = allPValsPhospho[["cluster1"]],
                                            logFCVector = allLog2FCPhospho[["cluster1"]])

plotProcessBars3(processesArr = allPhosphoProcesses$cluster1,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 1)",
                 processesExplainedList = phosphocluster1APE,
                 obviate.log2FC = TRUE)

phosphocluster2APE <- allProcessesExplained(
  refinedArr = refinedPhosphoMatrix,
  descript = completeDescription,
  chosenPValue = 0.05,
  pValuesVector = allPValsPhospho[["cluster2"]],
  logFCVector = allLog2FCPhospho[["cluster2"]]
  )

plotProcessBars3(processesArr = allPhosphoProcesses$cluster2,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 2)",
                 processesExplainedList = phosphocluster2APE,
                 obviate.log2FC = TRUE)

phosphocluster3APE <- allProcessesExplained(
  refinedArr = refinedPhosphoMatrix,
  descript = completeDescription,
  chosenPValue = 0.05,
  pValuesVector = allPValsPhospho[["cluster3"]],
  logFCVector = allLog2FCPhospho[["cluster3"]]
)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster3,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 3)",
                 processesExplainedList = phosphocluster3APE,
                 obviate.log2FC = TRUE)

phosphocluster4APE <- allProcessesExplained(
  refinedArr = refinedPhosphoMatrix,
  descript = completeDescription,
  chosenPValue = 0.05,
  pValuesVector = allPValsPhospho[["cluster4"]],
  logFCVector = allLog2FCPhospho[["cluster4"]]
)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster4,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 4)",
                 processesExplainedList = phosphocluster4APE,
                 obviate.log2FC = TRUE)

phosphocluster5APE <- allProcessesExplained(
  refinedArr = refinedPhosphoMatrix,
  descript = completeDescription,
  chosenPValue = 0.05,
  pValuesVector = allPValsPhospho[["cluster5"]],
  logFCVector = allLog2FCPhospho[["cluster5"]]
)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster5,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 5)",
                 processesExplainedList = phosphocluster5APE,
                 obviate.log2FC = TRUE)

phosphocluster6APE <- allProcessesExplained(
  refinedArr = refinedPhosphoMatrix,
  descript = completeDescription,
  chosenPValue = 0.05,
  pValuesVector = allPValsPhospho[["cluster6"]],
  logFCVector = allLog2FCPhospho[["cluster6"]]
)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster6,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Cluster 6)",
                 processesExplainedList = phosphocluster6APE,
                 obviate.log2FC = TRUE)

# Combined clusters

phosphocluster13APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["clusters13"]],
                                             logFCVector = allLog2FCPhospho[["clusters13"]])

plotProcessBars3(processesArr = allPhosphoProcesses$clusters13,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Clusters 1 and 3)",
                 processesExplainedList = phosphocluster13APE,
                 obviate.log2FC = TRUE)

phosphocluster25APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["clusters25"]],
                                             logFCVector = allLog2FCPhospho[["clusters25"]])

plotProcessBars3(processesArr = allPhosphoProcesses$clusters25,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Clusters 2 and 5)",
                 processesExplainedList = phosphocluster25APE,
                 obviate.log2FC = TRUE)

phosphocluster36APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["clusters36"]],
                                             logFCVector = allLog2FCPhospho[["clusters36"]])

plotProcessBars3(processesArr = allPhosphoProcesses$clusters36,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Clusters 3 and 6)",
                 processesExplainedList = phosphocluster36APE,
                 obviate.log2FC = TRUE)


phosphocluster45APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["clusters45"]],
                                             logFCVector = allLog2FCPhospho[["clusters45"]])

plotProcessBars3(processesArr = allPhosphoProcesses$clusters45,
                 howMany = 10,
                 main2 = "Top GO Biological Processes (Clusters 4 and 5)",
                 processesExplainedList = phosphocluster45APE,
                 obviate.log2FC = TRUE)
