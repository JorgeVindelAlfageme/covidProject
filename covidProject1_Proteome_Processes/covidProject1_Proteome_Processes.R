# 1) Loading packages, data, and separating proteomics values from protein
# descriptions

# sarsPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/lessData.xlsx"
sarsPath <- "C:/Users/jvalf/Desktop/Archivos/trabajo2/partOfFornax/fatima/lessData.xlsx"
library(xlsx)
library(limma)
sarsData <- read.xlsx2(file = sarsPath, sheetIndex = 1)
rownames(sarsData) <- sarsData[,1]
description <- sarsData[, c(1,2)]
onlyDescription <- description[, 2]

# 2) Adding new data to description matrix

description[, 3] <- retrieveGeneID(description[,2])
colnames(description)[3] <- "gene.ID"
anyNA(description[,3])

# 3) Deleting protein names and descriptions from the originally loaded data
# matrix, shortening its column names and reordering these columns as the user
# pleased

sarsData <- sarsData[,-c(1, 2)]
cols <- colnames(sarsData)
cols <- gsub(pattern = "Abundance..F2..", replacement = "", x = cols, fixed = TRUE)
cols <- gsub(pattern = "..Sample..", replacement = ".", x = cols, fixed = TRUE)
cols <- gsub(pattern = "..0", replacement = ".t0", x = cols, fixed = TRUE)
cols <- gsub(pattern = "CoV2..", replacement = "CoV2.t", x = cols, fixed = TRUE)
colnames(sarsData) <- cols
sarsData <- sarsData[, c(1:3,7:9,10:15,4:6)]
cols <- colnames(sarsData)
rows <- rownames(sarsData)
sarsData <- apply(X = sarsData, MARGIN = 2, FUN = as.numeric)
rownames(sarsData) <- rows
colnames(sarsData)[1:3] <- c("126.MOCK.t0", "129N.MOCK.t0", "131C.MOCK.t0")
cols <- colnames(sarsData)

# 4) Normalizing data using the quantile data standardization (not recommended)

sarsData[sarsData == 0] <- 1
normSARS <- normalizeBetweenArrays(object = log2(na.omit(sarsData)), method = "quantile")

# 5) Generating a vector of the different classes the data set represents

uniqueGroups <- c("MOCK.t0",
                  "SARS.CoV2.t3",
                  "SARS.CoV2.t6",
                  "SARS.CoV2.t9",
                  "SARS.CoV2.t16")
groups <- rep(x = uniqueGroups, each = 3)

# 6) ANOVA calculations

ANOVASARS <- multiANOVA3(dat = normCopy, classes = as.factor(groups))
ANOVASARS <- ANOVASARS[order(ANOVASARS[,2], decreasing = FALSE),]

adjANOVASARS <- adjANOVA(ANOVAres = ANOVASARS)

# 7) Analyzing altered processes considering significant proteins

allProteins <- rownames(normSARS)
pvec <- rep(0.01, length(allProteins))

# "UniProtR" was a library upon which the user-defined function "bioProcesses2"
# was defined. It can be possible that this function can't be run due to the
# changes that happened in the library. The result was saved in a variable that
# can be loaded, even if the original function can't be executed.

referenceMatrix <- bioProcesses2(uniProtAccessionVector = allProteins,
                                 pvalVec = pvec)

t16Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,3] <= 0.05]
t9Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,6] <= 0.05]
t6Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,5] <= 0.05]
t3Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,4] <= 0.05]

refinedMatrix <- referenceMatrix[rowSums(referenceMatrix) > 0,]

################################################################################

load(file = "C:/Users/usuario/Downloads/fatPaper/fatima/rda/sarsData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/logMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/proteinsPerCluster.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDMockSig.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/ttestPD.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/listClustersProcesses.rda")
load(file = "C:/Users/usuario/Downloads/fatPaper/fatima/rda/descriptionSARS.rda")

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

# APE: "All Processes Explained"

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

plotProcessBars3(processesArr = listClustersProcesses[["cluster1Processes"]],
                 howMany = 10,
                 main2 = "Top GO Biological processes (Cluster 1)",
                 processesExplainedList = cluster1APE)

# Clusters 2 and 4: p-value (MOCK vs. t16)

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

plotProcessBars3(processesArr = listClustersProcesses$combined2_4Processes,
                 howMany = 10,
                 main2 = "Top GO Biological processes (Clusters 2 and 4)",
                 processesExplainedList = clusters2And4APE)

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

plotProcessBars3(processesArr = listClustersProcesses[["cluster3Processes"]],
                 howMany = 10,
                 main2 = "Top GO Biological processes (Cluster 3)",
                 processesExplainedList = cluster3APE)
