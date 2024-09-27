library(xlsx)
library(readxl)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")

load(file = "C:/Users/usuario/Downloads/fatPaper/fatima/rda/refinedMatrix.rda")

clustersExcel <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/chosenT/mock/clusters.xlsx"

clustersExcel <- "C:/Users/usuario/Downloads/fatPaper/fatima/7_Mfuzz/chosenT/mock/clusters.xlsx"

nSheets <- excel_sheets(clustersExcel)
proteinsPerCluster <- list()
for(i in 1:length(nSheets)){
  proteinsPerCluster[[nSheets[i]]]  <- read.xlsx2(file = clustersExcel, sheetIndex = i)[,1]
}

save(x = proteinsPerCluster,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")

allClustersProcesses <- newAlteredProcesses(refinedArr = refinedMatrix,
                                            proteinVec = unlist(proteinsPerCluster))
cluster1Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                         proteinVec = proteinsPerCluster[["Cluster1"]])
cluster2Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                         proteinVec = proteinsPerCluster[["Cluster2"]])
cluster3Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                         proteinVec = proteinsPerCluster[["Cluster3"]])
cluster4Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                         proteinVec = proteinsPerCluster[["Cluster4"]])
combined2_4Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                            proteinVec = c(proteinsPerCluster[["Cluster2"]],
                                                           proteinsPerCluster[["Cluster4"]]))
listClustersProcesses <- list("allClusters" = allClustersProcesses,
                              "cluster1Processes"= cluster1Processes,
                              "cluster2Processes" = cluster2Processes,
                              "cluster3Processes" = cluster3Processes,
                              "cluster4Processes" = cluster4Processes,
                              "combined2_4Processes" = combined2_4Processes)
save(x = listClustersProcesses,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")

################################################################################

# Calculate log2(anyTime vs MOCK)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAbundances.rda")
PDAbundances[PDAbundances == 0] <- 1
PDMean <- matrix(data = c(apply(X = PDAbundances[,1:3], MARGIN = 1, FUN = mean),
                          apply(X = PDAbundances[,4:6], MARGIN = 1, FUN = mean),
                          apply(X = PDAbundances[,7:9], MARGIN = 1, FUN = mean),
                          apply(X = PDAbundances[,10:12], MARGIN = 1, FUN = mean),
                          apply(X = PDAbundances[,13:15], MARGIN = 1, FUN = mean)),
                 ncol = 5)
colnames(PDMean) <- c("MOCK.t0",
                      "SARS.CoV2.t3",
                      "SARS.CoV2.t6",
                      "SARS.CoV2.t9",
                      "SARS.CoV2.t16")
rownames(PDMean) <- rownames(PDAbundances)
logMatrix <- matrix(data = c(log2(PDMean[,2]/PDMean[,1]),
                             log2(PDMean[,3]/PDMean[,1]),
                             log2(PDMean[,4]/PDMean[,1]),
                             log2(PDMean[,5]/PDMean[,1])),
                    ncol = 4)
colnames(logMatrix) <- c("log2FC.t3.t0",
                         "log2FC.t6.t0",
                         "log2FC.t9.t0",
                         "log2FC.t16.t0")
rownames(logMatrix) <- rownames(PDAbundances)
save(x = logMatrix,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")

################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDMockSig.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/ttestPD.rda")

# Order newAlteredProcesses function results.

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")

completeProcessesTable <- matrix(data = NA,
                                 nrow = nrow(listClustersProcesses$allClusters),
                                 ncol = 1 + (ncol(listClustersProcesses[[1]]) - 1) * length(listClustersProcesses))
cols <- colnames(listClustersProcesses[[1]])[-1]
d <- rep(cols, each = 6)
e <- c(".allClusters", ".cluster1", ".cluster2", ".cluster3", ".cluster4", ".clusters2And4")
colnames(completeProcessesTable) <- c("total.proteins", sprintf(paste(d, "%s", sep = ""), e))

r <- rownames(listClustersProcesses[["allClusters"]])
rownames(completeProcessesTable) <- r

k <- 2
for(i in 2:ncol(listClustersProcesses[[1]])){
  for(j in 1:length(listClustersProcesses)){
    completeProcessesTable[,k] <- listClustersProcesses[[j]][r,i]
    k <- k + 1
  }
}
completeProcessesTable[,1] <- listClustersProcesses[[1]][r,1]

save(x = completeProcessesTable,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/completeProcessesTable.rda")

subTable1AlteredProteins <- completeProcessesTable[,1:7]
subTable2ProteinRates <- completeProcessesTable[,c(1,8:13)]
subTable3PValues <- completeProcessesTable[,14:19]
subTable4QValues <- completeProcessesTable[,20:25]

write.xlsx2(x = completeProcessesTable,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/completeProcessesTable.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allInformation",
            append = FALSE)

write.xlsx2(x = subTable1AlteredProteins,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/completeProcessesTable.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "alteredProteins",
            append = TRUE)

write.xlsx2(x = subTable2ProteinRates,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/subTable2ProteinRates.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "proteinRates",
            append = FALSE)

write.xlsx2(x = subTable3PValues,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/subTable3PValues.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "pValues",
            append = FALSE)

write.xlsx2(x = subTable4QValues,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/subTable4QValues.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "qValues",
            append = FALSE)

################################################################################

# APE: "All Processes Explained"

# Que para el cluster, el resto de prote?nas para un contraste de p-valor
# elegido tengan un p-valor no significativo.

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/sarsData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")
colnames(description) <- c("accession", "description", "gene.ID")

# Cluster 1: p-value (MOCK vs. t16)

cluster1PValues <- rep(NA, nrow(refinedMatrix))
names(cluster1PValues) <- rownames(refinedMatrix)
cluster1PValues[proteinsPerCluster[["Cluster1"]]] <- 0.01
cluster1PValues[is.na(cluster1PValues)] <- 1

cluster1APE <- allProcessesExplained(refinedArr = refinedMatrix,
                                     descript = description,
                                     chosenPValue = 0.05,
                                     pValuesVector = cluster1PValues,
                                     logFCVector = logMatrix[,"log2FC.t16.t0"])

printTopProcesses(processesMatrix = listClustersProcesses$cluster1Processes,
                  processesExplained = cluster1APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/topProcessesC1.xlsx",
                  howMany = 10)

write.xlsx2(x = cluster1Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/cluster1Processes.xlsx",
            sheetName = "cluster1",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

write.xlsx2(x = cluster2Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/cluster2Processes.xlsx",
            sheetName = "cluster2",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = cluster4Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/cluster4Processes.xlsx",
            sheetName = "cluster4",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = combined2_4Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/clusters2And4Processes.xlsx",
            sheetName = "clusters2And4",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

write.xlsx2(x = cluster3Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/cluster3Processes.xlsx",
            sheetName = "cluster3",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

plotProcessBars2(processesArr = listClustersProcesses[["cluster1Processes"]],
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz cluster 1 proteins)",
                 processesExplainedList = cluster1APE)

plotProcessBars3(processesArr = listClustersProcesses[["cluster1Processes"]],
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz cluster 1 proteins)",
                 processesExplainedList = cluster1APE)

subCluster1APE <- cluster1APE[rownames(listClustersProcesses[["cluster1Processes"]][1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/GO/stringnets/"

for(i in c(2:5,7:10)){
  
  variablePathPart <- paste("process", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subCluster1APE[[i]]$q.value,
            log2FC = subCluster1APE[[i]]$log2.FC,
            proteinNames = subCluster1APE[[i]]$protein.code,
            descript = description,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath)
  
}

################################################################################

cluster1STRINGResults <- list()
isExcel2 <- FALSE

for(i in c(2,3,7)){
  
  variablePathPart <- paste("process", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  cluster1STRINGResults[[i]] <- stringnet(qvalue = subCluster1APE[[i]]$q.value,
                                          log2FC = subCluster1APE[[i]]$log2.FC,
                                          proteinNames = subCluster1APE[[i]]$protein.code,
                                          descript = description,
                                          clusteringAlgorithm = "fastgreedy",
                                          isSmallNet = TRUE,
                                          STRINGVer = "11.5",
                                          selectedSpecies = 9606, # human
                                          scoreThres = 400,
                                          chosenPVal = 1,
                                          isExcel = isExcel2,
                                          filesPath = completePath)
  
}

names(cluster1STRINGResults) <- sprintf("process%s", c(2,3,7))

bigClusters <- list()
for(i in 1:length(cluster1STRINGResults)){
  iclusters <- list()
  imatrixProcess <- cluster1STRINGResults[[i]]$clustersMatrix
  for(j in 1:ncol(imatrixProcess)){
    if(colSums(imatrixProcess)[j] > 1){
      ibool <- imatrixProcess[,j] == 1
      iclusters[[j]] <- rownames(imatrixProcess)[ibool]
    }
  }
  names(iclusters) <- sprintf("subcluster%s", 1:length(iclusters))
  bigClusters[[i]] <- iclusters
}
names(bigClusters) <- names(cluster1STRINGResults)

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/GO/stringnets/"

for(i in 1:length(bigClusters)){
  for(j in 1:length(bigClusters[[i]])){
    print(paste("Process ", as.numeric(i), ", subprocess ", as.numeric(j)))
    proteinSet <- bigClusters[[i]][[j]]
    stringnet(qvalue = subCluster1APE[[i]]$q.value[proteinSet,],
              log2FC = subCluster1APE[[i]]$log2.FC[proteinSet,],
              proteinNames = proteinSet,
              descript = description,
              clusteringAlgorithm = "fastgreedy",
              isSmallNet = TRUE,
              STRINGVer = "11.5",
              selectedSpecies = 9606, # human
              scoreThres = 400,
              chosenPVal = 1,
              isExcel = isExcel2,
              filesPath = completePath)
  }
}

################################################################################

# KEGG (cluster 1)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/sarsData.rda")

cluster1Log2FC <- logMatrix[,"log2FC.t16.t0"]
cluster1Bool <- names(cluster1Log2FC) %in% names(which(cluster1PValues <= 0.5))
cluster1Bool <- !(cluster1Bool)
cluster1Log2FC[cluster1Bool] <- 0

cluster1Data <- sarsData[,c(1:3,13:15)]
cluster1Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
cluster1KEGG <- keggProcesses(rawData = cluster1Data,
                              classesVector = cluster1Classes,
                              keyType = "UNIPROT",
                              is.pval.and.log = TRUE,
                              pval = cluster1PValues,
                              qval = cluster1PValues,
                              log2FC = cluster1Log2FC)

View(cluster1KEGG$bioProcesses$greater)
cluster1KEGG$log2FC
View(cluster1KEGG$hsaList$`hsa04914 Progesterone-mediated oocyte maturation`)

setworkspace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/KEGG/"
keggNetworks(keggres = cluster1KEGG$bioProcesses,
             foldchanges = cluster1KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setworkspace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs)

excelKeggRes <- cluster1KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("cluster1.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- cluster1KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(cluster1KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(cluster1KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/KEGG/allKEGGProcesses.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

cluster1TopProcesses <- cluster1KEGG$hsaList[rownames(cluster1KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = cluster1TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = cluster1TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster1/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
  
}

# Clusters 2 and 4: p-value (MOCK vs. t16)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/sarsData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")

colnames(description)[c(1, 2)] <- c("accesion", "description")

clusters2And4PValues <- rep(NA, nrow(refinedMatrix))
names(clusters2And4PValues) <- rownames(refinedMatrix)
proteinsClusters2And4 <- c(proteinsPerCluster[["Cluster2"]], proteinsPerCluster[["Cluster4"]])
clusters2And4PValues[proteinsClusters2And4] <- 0.01
clusters2And4PValues[is.na(clusters2And4PValues)] <- 1

clusters2And4APE <- allProcessesExplained(refinedArr = refinedMatrix,
                                          descript = description,
                                          chosenPValue = 0.05,
                                          pValuesVector = clusters2And4PValues,
                                          logFCVector = logMatrix[,"log2FC.t16.t0"])

printTopProcesses(processesMatrix = listClustersProcesses$combined2_4Processes,
                  processesExplained = clusters2And4APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/clusters2And4TopProcesses.xlsx",
                  howMany = 10)

write.xlsx2(x = listClustersProcesses$combined2_4Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/clusters2AndProcesses.xlsx",
            sheetName = "clusters2And4",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = listClustersProcesses$cluster2Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/cluster2Processes.xlsx",
            sheetName = "cluster2",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = listClustersProcesses$cluster4Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/cluster4Processes.xlsx",
            sheetName = "cluster4",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

plotProcessBars2(processesArr = listClustersProcesses$combined2_4Processes,
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz clusters 2 and 4 proteins combined)",
                 processesExplainedList = clusters2And4APE)

plotProcessBars3(processesArr = listClustersProcesses$combined2_4Processes,
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz clusters 2 and 4 proteins combined)",
                 processesExplainedList = clusters2And4APE)

subClusters2And4APE <- clusters2And4APE[rownames(listClustersProcesses[["combined2_4Processes"]][1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/GO/stringnets/"
isExcel2 <- FALSE

for(i in c(1:4, 6:10)){
  
  variablePathPart <- paste("process", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subClusters2And4APE[[i]]$q.value,
            log2FC = subClusters2And4APE[[i]]$log2.FC,
            proteinNames = subClusters2And4APE[[i]]$protein.code,
            descript = description,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606,
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = isExcel2,
            filesPath = completePath)
  
}

# KEGG (clusters 2 and 4)

clusters2And4Log2FC <- logMatrix[,"log2FC.t16.t0"]
clusters2And4Bool <- names(clusters2And4Log2FC) %in% names(which(clusters2And4PValues <= 0.5))
clusters2And4Bool <- !(clusters2And4Bool)
clusters2And4Log2FC[clusters2And4Bool] <- 0

clusters2And4Data <- sarsData[,c(1:3,13:15)]
clusters2And4Classes <- c(rep("MOCK.t0",3), rep("SARS.CoV2.t16", 3))
clusters2And4KEGG <- keggProcesses(rawData = clusters2And4Data,
                                   classesVector = clusters2And4Classes,
                                   keyType = "UNIPROT",
                                   is.pval.and.log = TRUE,
                                   pval = clusters2And4PValues,
                                   qval = clusters2And4PValues,
                                   log2FC = clusters2And4Log2FC)

View(clusters2And4KEGG$bioProcesses$greater)
clusters2And4KEGG$log2FC
View(clusters2And4KEGG$hsaList$`hsa05012 Parkinson's disease`)

setworkspace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/KEGG/images/"
keggNetworks(keggres = clusters2And4KEGG$bioProcesses,
             foldchanges = clusters2And4KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setworkspace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs)

excelKeggRes <- clusters2And4KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("clusters2And4.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- clusters2And4KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(clusters2And4KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(clusters2And4KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/KEGG/allKEGGProcesses.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

clusters2And4TopProcesses <- clusters2And4KEGG$hsaList[rownames(clusters2And4KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = clusters2And4TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = clusters2And4TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/clusters2And4/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
  
}

################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/sarsData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")

# Cluster 3: p-value (MOCK vs. t3)

cluster3PValues <- rep(NA, nrow(refinedMatrix))
names(cluster3PValues) <- rownames(refinedMatrix)
cluster3PValues[proteinsPerCluster[["Cluster3"]]] <- 0.01
cluster3PValues[is.na(cluster3PValues)] <- 1

cluster3APE <- allProcessesExplained(refinedArr = refinedMatrix,
                                     descript = description,
                                     chosenPValue = 0.05,
                                     pValuesVector = cluster3PValues,
                                     logFCVector = logMatrix[,"log2FC.t3.t0"])

printTopProcesses(processesMatrix = listClustersProcesses$cluster3Processes,
                  processesExplained = cluster3APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/cluster3TopProcesses.xlsx",
                  howMany = 10)

write.xlsx2(x = cluster3Processes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/cluster3Processes.xlsx",
            sheetName = "cluster3",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

plotProcessBars2(processesArr = listClustersProcesses[["cluster3Processes"]],
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz cluster 3 proteins)",
                 processesExplainedList = cluster3APE)

plotProcessBars3(processesArr = listClustersProcesses[["cluster3Processes"]],
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz cluster 3 proteins)",
                 processesExplainedList = cluster3APE)

subCluster3APE <- cluster3APE[rownames(listClustersProcesses[["cluster3Processes"]][1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/GO/stringnets/"
boolVec <- c(rep(TRUE, 3), FALSE, rep(TRUE, 6))

for(i in 1:10){
  
  variablePathPart <- paste("process", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subCluster3APE[[i]]$q.value,
            log2FC = subCluster3APE[[i]]$log2.FC,
            proteinNames = subCluster3APE[[i]]$protein.code,
            descript = description,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = boolVec[i],
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath)
  
}

# KEGG (cluster 3)

cluster3Log2FC <- logMatrix[,"log2FC.t3.t0"]
# cluster3Bool <- names(cluster3Log2FC) %in% names(which(cluster3PValues <= 0.5))
# cluster3Bool <- !(cluster3Bool)
# cluster3Log2FC[cluster3Bool] <- 0

cluster3Data <- sarsData[,c(1:3,4:6)]
cluster3Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t3", 3))
cluster3KEGG <- keggProcesses(rawData = cluster3Data,
                              classesVector = cluster3Classes,
                              keyType = "UNIPROT",
                              is.pval.and.log = TRUE,
                              pval = cluster3PValues,
                              qval = cluster3PValues,
                              log2FC = cluster3Log2FC)

View(cluster3KEGG$bioProcesses$greater)
cluster3KEGG$log2FC
# View(cluster3KEGG$hsaList$)

setworkspace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/KEGG/images/"
keggNetworks(keggres = cluster3KEGG$bioProcesses,
             foldchanges = cluster3KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setworkspace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs)

excelKeggRes <- cluster3KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("cluster3.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- cluster3KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(cluster3KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(cluster3KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/KEGG/allKEGGProcesses.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

cluster3TopProcesses <- cluster3KEGG$hsaList[rownames(cluster3KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = cluster3TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = cluster3TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/cluster3/KEGG/topKEGGProcesses.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
  
}


