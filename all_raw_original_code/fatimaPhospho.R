load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normCopy.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")

library(xlsx)
library(readxl)
sigRClusters <- list()
excelFilePath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_Sig/clustersRSig.xlsx"
nSheets <- length(excel_sheets(excelFilePath))

for(i in 1:nSheets){
  sigRClusters[[i]] <- read.xlsx2(file = excelFilePath, sheetIndex = i)[,1]
}
# sigRClusters[[10]] <- read.xlsx2(file = excelFilePath, sheetIndex = 10)[,1]
names(sigRClusters) <- sprintf("Cluster%s", 1:10)

pvalMatrix <- adjANOVASARS[,3:12]
newTitles <- reducedDescription(fullDescription = description)

for(i in 1:length(sigRClusters)){
  iprots <- sigRClusters[[i]]
  idir <- names(sigRClusters)[i]
  photoBox2(is.already = FALSE,
            nameDir = idir,
            where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_sig/boxAndDesc/",
            dat = normCopy[iprots,],
            pvalues = pvalMatrix[iprots,],
            groups2 = groups,
            is.ANOVA = TRUE,
            main3 = newTitles[iprots])
}

################################################################################

# An?lisis de fosfoprote?mica:

library(xlsx)

phosphoExcel <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoReduced.xlsx"

phosphoExcel <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/20220801_SARS_Cells_TMT16plex_Fosfo_Time_t-test_All comparison_Listofpeptides.xlsx"

phosphoDescription <- read.xlsx2(file = phosphoExcel, sheetIndex = 1)
phosphoRaw <- read.xlsx2(file = phosphoExcel, sheetIndex = 2)
phosphoRaw <- phosphoRaw[, -ncol(phosphoRaw)]

phosphoAll <- cbind(phosphoRaw[, c(1,2)],
                    phosphoDescription[, c(2, 4:ncol(phosphoDescription))],
                    phosphoRaw[, 3:ncol(phosphoRaw)])

phosphoCols <- colnames(phosphoAll)
phosphoCols <- gsub(pattern = "...", replacement = ".", x = phosphoCols, fixed = TRUE)
phosphoCols <- gsub(pattern = "..", replacement = ".", x = phosphoCols, fixed = TRUE)
phosphoCols <- gsub(pattern = "X", replacement = "", x = phosphoCols, fixed = TRUE)
phosphoCols <- gsub(pattern = "MOCK", replacement = "MOCK.0h", x = phosphoCols, fixed = TRUE)
colnames(phosphoAll) <- phosphoCols

chosenData <- phosphoAll

phosphoConverter <- function(chosenData){
  
  # Make sure the first column indicates the peptide, and the second column,
  # the proteins that have that peptide. The rest of the column have to contain
  # the data that is going to be replicated.
  
  backUpCols <- c(colnames(chosenData)[1:2],
                  "Is.shared", "Sharing.proteins",
                  colnames(chosenData)[3:ncol(chosenData)])
  
  emptyShell <- matrix(data = NA, nrow = nrow(chosenData), ncol = ncol(chosenData) + 2)
  
  i <- 1
  k <- 1
  
  while((nrow(chosenData) + 1) != i){
    
    iproteins <- unlist(strsplit(x = chosenData[i,2], split = "; ", fixed = TRUE))
    
    if(length(iproteins) > 1){
      addToShell <- matrix(data = NA,
                           nrow = length(iproteins) - 1,
                           ncol = ncol(emptyShell))
      emptyShell <- rbind(emptyShell, addToShell)
      for(j in 1:length(iproteins)){
        emptyShell[(k + j - 1), 1] <- chosenData[i, 1]
        emptyShell[(k + j - 1), 2] <- iproteins[j]
        emptyShell[(k + j - 1), 3] <- TRUE
        emptyShell[(k + j - 1), 4] <- chosenData[i,2]
        emptyShell[(k + j - 1), 5:ncol(emptyShell)] <- as.character(chosenData[i, 3:ncol(chosenData)])
      }
      k <- k + j
    }else{
      emptyShell[k,] <- c(as.character(chosenData[i,c(1,2)]), as.character(FALSE), as.character(chosenData[i,2:ncol(chosenData)]))
      k <- k + 1
    }
    
    i <- i + 1
    
  }
  
  shellRowNames <- sprintf(paste(emptyShell[,1], "..", "%s", sep = ""),
                           emptyShell[,2])
  shellRowNames <- gsub(pattern = "[",
                        replacement = "",
                        x = shellRowNames,
                        fixed = TRUE)
  shellRowNames <- gsub(pattern = "]",
                        replacement = "",
                        x = shellRowNames,
                        fixed = TRUE)
  rownames(emptyShell) <- shellRowNames
  colnames(emptyShell) <- backUpCols
  
  return(emptyShell)
  
}

newPhosphoAll <- phosphoConverter(chosenData = chosenData)
newPhosphoAll <- phosphoConverter(chosenData = phosphoAll)
save(x = newPhosphoAll,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/newPhosphoAll.rda")

load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/newPhosphoAll.rda")

# Ordenar el resultado

newPhosphoRows <- rownames(newPhosphoAll)

phosphoPVals <- newPhosphoAll[ ,6:9]
phosphoPVals <- apply(X = phosphoPVals, MARGIN = 2, FUN = as.numeric)
rownames(phosphoPVals) <- newPhosphoRows
save(x = phosphoPVals,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoPVals.rda")

phosphoData <- newPhosphoAll[,10:24]
phosphoData <- apply(X = phosphoData, MARGIN = 2, FUN = as.numeric)
rownames(phosphoData) <- newPhosphoRows
save(x = phosphoData,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoData.rda")

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")

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
save(x = phosphoDescription,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoDescription.rda")

################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoDescription.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoPVals.rda")

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

save(x = phosphoDescription,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoDescUpdated.rda")

################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoData.rda")
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

save(x = aggSum,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/aggSum.rda")

################################################################################

library(xlsx)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoDescUpdated.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/aggSum.rda")

allPhosphoDataPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/allR1.xlsx"
allPhosphoDescPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/allR2.xlsx"

allPhosphoDataPath <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/allR1.xlsx"
allPhosphoDescPath <- "C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/allR2.xlsx"

allRawPhosphoData <- read.xlsx2(file = allPhosphoDataPath, sheetIndex = 1)
descAllPhosphoData <- read.xlsx2(file = allPhosphoDescPath, sheetIndex = 1)

# Solo fosfop?ptidos, con prote?na de origen reconocida y sin compartir su
# secuencia con otra prote?na... ?Algo m?s?

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

save(x = filteredData,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/filteredData.rda")
save(x = filteredDesc,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/filteredDesc.rda")

################################################################################

sigPhosphoPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoReduced.xlsx"

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

save(x = sigData,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/sigData.rda")
save(x = sigDesc,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/sigDesc.rda")
save(x = aggSum,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/aggSum.rda")

################################################################################

phosphoFuzz <- mfuzzing(aggregatedResult = aggSum, chosenRange = 2:10)
save(x = phosphoFuzz,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/phosphoFuzz.rda")

resPhosphoFuzz <- clAndPlot(mfuzzRes = phosphoFuzz,
                            chosenC = 6,
                            chosenMF = c(3, 2),
                            chosenTimes = c(0, 3, 6, 9, 16),
                            xlab2 = "Time (h)")

load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/resPhosphoFuzz.rda")

save(x = resPhosphoFuzz,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/resPhosphoFuzz.rda")

sigPhosphoProteinsClusters <- list()
for(i in 1:max(resPhosphoFuzz$cluster)){
  inames <- names(resPhosphoFuzz$cluster)[resPhosphoFuzz$cluster == i]
  iprots <- sigDesc[inames, 3]
  clusterName <- paste("cluster", as.numeric(i), sep = "")
  sigPhosphoProteinsClusters[[clusterName]] <- iprots
}
save(x = sigPhosphoProteinsClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/sigPhosphoProteinsClusters.rda")

# Retrieving significant proteome elements

proteomeClustersProteins <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/chosenT/mock/clusters.xlsx"
library(readxl)
library(xlsx)
proteomeClusters <- list()
for(i in 1:length(excel_sheets(path = proteomeClustersProteins))){
  
  iproteins  <- read.xlsx2(file = proteomeClustersProteins, sheetIndex = i)[,1]
  clusterName <- paste("cluster", as.numeric(i), sep = "")
  proteomeClusters[[clusterName]] <- iproteins
  
}
save(x = proteomeClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/proteomeClusters.rda")

# Boxplots; excel tables

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
descriptionSARS <- description
rm(description)
newSigDesc <- cbind.data.frame(sigDesc,descriptionSARS[sigDesc[,3],c(2,3)])
save(x = newSigDesc,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/newSigDesc.rda")

phosphoClustersExcel <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phospho_mfuzz/phosphoClustersDesc.xlsx"

for(i in 1:max(resPhosphoFuzz$cluster)){
  
  sigIDs <- names(resPhosphoFuzz$cluster)[resPhosphoFuzz$cluster == i]
  isheet <- paste("cluster", as.numeric(i), sep = "")
  if(i == 1){
    write.xlsx2(x = newSigDesc[sigIDs,],
                file = phosphoClustersExcel,
                row.names = TRUE,
                col.names = TRUE,
                sheetName = isheet,
                append = FALSE)
  }else{
    write.xlsx2(x = newSigDesc[sigIDs,],
                file = phosphoClustersExcel,
                row.names = TRUE,
                col.names = TRUE,
                sheetName = isheet,
                append = TRUE)
  }
  
}

phosphoBoxplotsTitles <- paste(newSigDesc[,"Modifications"],
                               " (",
                               newSigDesc[,3],
                               ")",
                               sep = "")

nGroups <- ncol(aggSum)
nContrasts <- nGroups * ((nGroups - 1)/2)
falsePMatrix <- matrix(data = 1, nrow = nrow(sigData), ncol = nContrasts)
rownames(falsePMatrix) <- rownames(sigData)

sameGroups <- c(rep("MOCK.t0", 3),
                rep("SARS.CoV2.t3", 3),
                rep("SARS.CoV2.t6", 3),
                rep("SARS.CoV2.t9", 3),
                rep("SARS.CoV2.t16", 3))


names(phosphoBoxplotsTitles) <- rownames(newSigDesc)
wherePhosphoBoxplots <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phospho_mfuzz/phosphoBoxplots/"

for(i in 1:max(resPhosphoFuzz$cluster)){
  iphospho <- names(resPhosphoFuzz$cluster)[resPhosphoFuzz$cluster == i]
  photoBox2(is.already = FALSE,
            nameDir = paste("cluster", as.numeric(i), sep = ""),
            where = wherePhosphoBoxplots,
            dat = log2(sigData[iphospho,]),
            pvalues = falsePMatrix[iphospho,],
            groups2 = sameGroups,
            is.ANOVA = TRUE,
            main3 = phosphoBoxplotsTitles[iphospho])
}

c("sig.ID.48", "sig.ID.326")

sigData[sigData == 0.0] <- 1 # !!!

multiBox2(protValues = log2(sigData["sig.ID.48",]),
          protName = "sig.ID.48",
          groups = sameGroups,
          pvals = rep(0, nContrasts),
          onlySig = TRUE,
          main2 = phosphoBoxplotsTitles["sig.ID.48"])

multiBox2(protValues = log2(sigData["sig.ID.326",]),
          protName = "sig.ID.326",
          groups = sameGroups,
          pvals = rep(0, nContrasts),
          onlySig = TRUE,
          main2 = phosphoBoxplotsTitles["sig.ID.326"])

################################################################################

# Venn

# rownames(PDMockSig)
# unique(unlist(sigPhosphoProteinsClusters))

library(VennDiagram)

draw.pairwise.venn(area1 = length(rownames(PDMockSig)),
                   area2 = length(unique(unlist(sigPhosphoProteinsClusters))),
                   cross.area = length(intersect(x = rownames(PDMockSig),
                                                 y = unique(unlist(sigPhosphoProteinsClusters)))),
                   category = c("Proteome","Phospho"),
                   fill = c("forestgreen","skyblue"))

################################################################################

allProteinsInPhospho <- unique(filteredDesc[,3])

refinedPhosphoMatrix <- bioProcesses2(uniProtAccessionVector = allProteinsInPhospho,
                                      pvalVec = rep(0.01, length(allProteinsInPhospho)))
save(x = refinedPhosphoMatrix,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")

load("C:/Users/usuario/Downloads/fatPaper/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")

allSigPhospho <- unique(unlist(sigPhosphoProteinsClusters))
phosphoAlteredProcesses <- newAlteredProcesses(refinedPhosphoMatrix, allSigPhospho)

c1.3 <- unique(c(sigPhosphoProteinsClusters$cluster1,
                 sigPhosphoProteinsClusters$cluster3))
c2.5 <- unique(c(sigPhosphoProteinsClusters$cluster2,
                 sigPhosphoProteinsClusters$cluster5))
# c2.4.5 <- unique(c(sigPhosphoProteinsClusters$cluster2,
#                    sigPhosphoProteinsClusters$cluster4,
#                    sigPhosphoProteinsClusters$cluster5))
c4.5 <- unique(c(sigPhosphoProteinsClusters$cluster4,
                 sigPhosphoProteinsClusters$cluster5))
c3.6 <- unique(c(sigPhosphoProteinsClusters$cluster3,
                 sigPhosphoProteinsClusters$cluster6))

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

# a <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 2)],
#                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC245[rownames(phosphoAlteredProcesses),c(2)],
#                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(2)])
# 
# colnames(a) <- c("total.proteins",
#                  sprintf("altered.proteins.%s",
#                          c("all",
#                            sprintf("cluster.%s", 1:6),
#                            sprintf("clusters.%s", c("2.5", "2.4.5", "3.6")))))

a <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 2)],
                      phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC13[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC45[rownames(phosphoAlteredProcesses),c(2)],
                      phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(2)])

colnames(a) <- c("total.proteins",
                 sprintf("altered.proteins.%s",
                         c("all",
                           sprintf("cluster.%s", 1:6),
                           sprintf("clusters.%s", c("1.3", "2.5", "4.5", "3.6")))))

View(a)

# b <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 3)],
#                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC245[rownames(phosphoAlteredProcesses),c(3)],
#                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(3)])
# 
# colnames(b) <- c("total.proteins",
#                  sprintf("altered.ratio.%s",
#                          c("all",
#                            sprintf("cluster.%s", 1:6),
#                            sprintf("clusters.%s", c("2.5", "2.4.5", "3.6")))))

b <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 3)],
                      phosphoAlteredC1[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC2[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC3[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC4[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC5[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC6[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC13[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC25[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC45[rownames(phosphoAlteredProcesses), c(3)],
                      phosphoAlteredC36[rownames(phosphoAlteredProcesses), c(3)])

colnames(b) <- c("total.proteins",
                 sprintf("altered.ratio.%s",
                         c("all",
                           sprintf("cluster.%s", 1:6),
                           sprintf("clusters.%s", c("1.3", "2.5", "4.5", "3.6")))))

# c <- cbind.data.frame(phosphoAlteredProcesses[,c(4)],
#                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC245[rownames(phosphoAlteredProcesses),c(4)],
#                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(4)])
# 
# colnames(c) <- c(sprintf("fisher.pvalue.%s",
#                          c("all",
#                            sprintf("cluster.%s", 1:6),
#                            sprintf("clusters.%s", c("2.5", "2.4.5", "3.6")))))

c <- cbind.data.frame(phosphoAlteredProcesses[,c(4)],
                      phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC13[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC45[rownames(phosphoAlteredProcesses),c(4)],
                      phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(4)])

colnames(c) <- c(sprintf("fisher.pvalue.%s",
                         c("all",
                           sprintf("cluster.%s", 1:6),
                           sprintf("clusters.%s", c("1.3", "2.5", "4.5", "3.6")))))

# d <- cbind.data.frame(phosphoAlteredProcesses[,c(5)],
#                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC245[rownames(phosphoAlteredProcesses),c(5)],
#                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(5)])
# 
# colnames(d) <- c(sprintf("fisher.qvalue.%s",
#                          c("all",
#                            sprintf("cluster.%s", 1:6),
#                            sprintf("clusters.%s", c("2.5", "2.4.5", "3.6")))))

d <- cbind.data.frame(phosphoAlteredProcesses[,c(5)],
                      phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC2[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC3[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC5[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC6[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC13[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC45[rownames(phosphoAlteredProcesses),c(5)],
                      phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(5)])

colnames(d) <- c(sprintf("fisher.qvalue.%s",
                         c("all",
                           sprintf("cluster.%s", 1:6),
                           sprintf("clusters.%s", c("1.3", "2.5", "4.5", "3.6")))))

library(xlsx)

newFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/phosphoProcessesAllInfo.xlsx"

write.xlsx2(x = a,
            file = newFile,
            sheetName = "numberOfProteins",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = b,
            file = newFile,
            sheetName = "proteinRates",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)
write.xlsx2(x = c,
            file = newFile,
            sheetName = "pvalue",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)
write.xlsx2(x = d,
            file = newFile,
            sheetName = "qvalue",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)

a2 <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 2)],
                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(2)],
                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(2)],
                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(2)],
                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(2)])

colnames(a2) <- c("total.proteins",
                  sprintf("altered.proteins.%s",
                          c("all", "cluster.1", "clusters.2.5", "clusters.3.6", "cluster.4")))

b2 <- cbind.data.frame(phosphoAlteredProcesses[,c(1, 3)],
                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(3)],
                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(3)],
                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(3)],
                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(3)])

colnames(b2) <- c("total.proteins",
                  sprintf("protein.rate.%s",
                          c("all", "cluster.1", "clusters.2.5", "clusters.3.6", "cluster.4")))

c2 <- cbind.data.frame(phosphoAlteredProcesses[,c(4)],
                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(4)],
                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(4)],
                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(4)],
                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(4)])

colnames(c2) <- c(sprintf("fisher.pvalue.%s",
                          c("all", "cluster.1", "clusters.2.5", "clusters.3.6", "cluster.4")))

d2 <- cbind.data.frame(phosphoAlteredProcesses[,c(5)],
                       phosphoAlteredC1[rownames(phosphoAlteredProcesses),c(5)],
                       phosphoAlteredC25[rownames(phosphoAlteredProcesses),c(5)],
                       phosphoAlteredC36[rownames(phosphoAlteredProcesses),c(5)],
                       phosphoAlteredC4[rownames(phosphoAlteredProcesses),c(5)])

colnames(d2) <- c(sprintf("fisher.qvalue.%s",
                          c("all", "cluster.1", "clusters.2.5", "clusters.3.6", "cluster.4")))

newFile2 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/perChosenClusters/allChosenClusters.xlsx"

write.xlsx2(x = a2,
            file = newFile2,
            sheetName = "numberOfProteins",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = b2,
            file = newFile2,
            sheetName = "proteinRates",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)
write.xlsx2(x = c2,
            file = newFile2,
            sheetName = "pvalue",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)
write.xlsx2(x = d2,
            file = newFile2,
            sheetName = "qvalue",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)

newFile3 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/perChosenClusters/phosphoCluster1.xlsx"
newFile4 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/perChosenClusters/phosphoClusters25.xlsx"
newFile5 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/perChosenClusters/phosphoCluster36.xlsx"
newFile6 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/perChosenClusters/phosphoCluster4.xlsx"

write.xlsx2(x = phosphoAlteredC1,
            file = newFile3,
            sheetName = "Hoja1",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = phosphoAlteredC25,
            file = newFile4,
            sheetName = "Hoja1",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = phosphoAlteredC36,
            file = newFile5,
            sheetName = "Hoja1",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)
write.xlsx2(x = phosphoAlteredC4,
            file = newFile6,
            sheetName = "Hoja1",
            col.names = TRUE,
            row.names = TRUE,
            append = FALSE)

newFile7 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/tables/perChosenClusters/proteinDescriptionC1.xlsx"
newFile8 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/tables/perChosenClusters/proteinDescriptionC25.xlsx"
newFile9 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/tables/perChosenClusters/proteinDescriptionC36.xlsx"
newFile10 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/tables/perChosenClusters/proteinDescriptionC4.xlsx"

# NECESITO UNA DESCRIPCI?N COMPLETA CON TODAS LAS PROTE?NAS

write.xlsx2(x = completeDescription[sigPhosphoProteinsClusters$cluster1,],
            file = newFile7,
            sheetName = "cluster1",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)
write.xlsx2(x = completeDescription[c2.5,],
            file = newFile8,
            sheetName = "clusters25",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)
write.xlsx2(x = completeDescription[c3.6,],
            file = newFile9,
            sheetName = "clusters36",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)
write.xlsx2(x = completeDescription[sigPhosphoProteinsClusters$cluster4,],
            file = newFile10,
            sheetName = "cluster4",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)

################################################################################

lackingProteins <- setdiff(allProteinsInPhospho, rownames(descriptionSARS))
library(xlsx)
lackingAll <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/lackingProteins.xlsx",
                         sheetIndex = 1)
lackingAll <- lackingAll[,-1]
colnames(lackingAll) <- c("accession", "description", "gene.ID")
completeDescriptionRows <- c(rownames(descriptionSARS), rownames(lackingAll))
completeDescription <- rbind(descriptionSARS, lackingAll)

save(x = completeDescription,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/completeDescription.xlsx")

################################################################################

save(x = sigPhosphoProteinsClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/sigPhosphoProteinsClusters.rda")

allPhosphoProcesses <- list("cluster1" = phosphoAlteredC1,
                            "cluster2" = phosphoAlteredC2,
                            "cluster3" = phosphoAlteredC3,
                            "cluster4" = phosphoAlteredC4,
                            "cluster5" = phosphoAlteredC5,
                            "cluster6" = phosphoAlteredC6,
                            "cluster25" = phosphoAlteredC25,
                            "cluster36" = phosphoAlteredC36)

save(x = allPhosphoProcesses,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPhosphoProcesses.rda")

allPValsPhospho <- list("cluster1" = c(),
                        "cluster25" = c(),
                        "cluster36" = c(),
                        "cluster4" = c())
allLog2FCPhospho <- allPValsPhospho

relevantClusters <- list("cluster1" = sigPhosphoProteinsClusters$cluster1,
                         "clusters25" = c2.5,
                         "clusters36" = c3.6,
                         "cluster4" = sigPhosphoProteinsClusters$cluster4)

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

save(x = allPValsPhospho,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPValsPhospho.rda")
save(x = allLog2FCPhospho,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allLog2FCPhospho.rda")

################################################################################

# Phosphocluster 1

phosphocluster1APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                            descript = completeDescription,
                                            chosenPValue = 0.05,
                                            pValuesVector = allPValsPhospho[["cluster1"]],
                                            logFCVector = allLog2FCPhospho[["cluster1"]])

printTopProcesses(processesMatrix = allPhosphoProcesses$cluster1,
                  processesExplained = phosphocluster1APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/GO/topPhosphoCluster1.xlsx",
                  howMany = 10)

plotProcessBars2(processesArr = allPhosphoProcesses$cluster1,
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz cluster 1 proteins)",
                 processesExplainedList = phosphocluster1APE,
                 obviate.log2FC = TRUE)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster1,
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz cluster 1 proteins)",
                 processesExplainedList = phosphocluster1APE,
                 obviate.log2FC = TRUE)

subPhospho1APE <- phosphocluster1APE[rownames(allPhosphoProcesses$cluster1[1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/GO/STRINGnets/"

for(i in c(2:7)){
  
  variablePathPart <- paste("phospho", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subPhospho1APE[[i]]$q.value,
            log2FC = subPhospho1APE[[i]]$log2.FC,
            proteinNames = subPhospho1APE[[i]]$protein.code,
            descript = completeDescription,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath,
            colSet = "greyblue",
            is.pval.title = TRUE)
  
}

rm(qvalue, log2FC, proteinNames, descript, clusteringAlgorithm, isSmallNet,
   STRINGVer, selectedSpecies, scoreThres, chosenPVal, isExcel, filesPath,
   colSet, is.pval.title, string_db, idata, mapped, translatorTable,
   stringClusters, clustersMatrix, matrixStartPoint, interactionsShell,
   tables, hits, mapped_pval, payload_id, insetLegendX, posGradientX1,
   posGradientX2, posGradient, textPosX, textPosY, enterYes, colorGradient,
   legend_image, rangeValue, op)

i <- 1

stringnet(qvalue = subPhospho1APE[[i]]$q.value[subPhospho1APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho1APE[[i]]$log2.FC[subPhospho1APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho1APE[[i]]$protein.code[subPhospho1APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/GO/STRINGnets/process1/phospho1.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 4

stringnet(qvalue = subPhospho1APE[[i]]$q.value[subPhospho1APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho1APE[[i]]$log2.FC[subPhospho1APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho1APE[[i]]$protein.code[subPhospho1APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/GO/STRINGnets/process4/phospho4.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)


save(x = relevantClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")
save(x = allProteinsInPhospho,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allProteinsInPhospho.rda")
save(x = completeDescription,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")

# KEGG

phospho1Data <- matrix(data = abs(rnorm(n = (length(allLog2FCPhospho$cluster1) * 6),
                                      mean = 1,
                                      sd = 1)),
                       nrow = length(allLog2FCPhospho$cluster1),
                       ncol = 6)
phospho1Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(phospho1Data) <- names(allLog2FCPhospho$cluster1)
colnames(phospho1Data) <- phospho1Classes

phospho1KEGG <- keggProcesses(rawData = phospho1Data,
                              classesVector = phospho1Classes,
                              keyType = "UNIPROT",
                              is.pval.and.log = TRUE,
                              pval = allPValsPhospho$cluster1,
                              qval = allPValsPhospho$cluster1,
                              log2FC = allLog2FCPhospho$cluster1)

View(phospho1KEGG$bioProcesses$greater)

setWorkSpace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/KEGG/"
keggNetworks(keggres = phospho1KEGG$bioProcesses,
             foldchanges = phospho1KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setWorkSpace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs,
             is.log2FC.absent = TRUE)

excelKeggRes <- phospho1KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("cluster1.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- phospho1KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(phospho1KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(phospho1KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/KEGG/phospho1KEGG.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

phospho1TopProcesses <- phospho1KEGG$hsaList[rownames(phospho1KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = phospho1TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/KEGG/phospho1TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = phospho1TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster1/KEGG/phospho1TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
}

################################################################################

# Clusters 2 and 5

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allProteinsInPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPValsPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allLog2FCPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPhosphoProcesses.rda")

phosphocluster25APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["cluster25"]],
                                             logFCVector = allLog2FCPhospho[["cluster25"]])

printTopProcesses(processesMatrix = allPhosphoProcesses$cluster25,
                  processesExplained = phosphocluster25APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/GO/topPhosphoClusters25.xlsx",
                  howMany = 10)

plotProcessBars2(processesArr = allPhosphoProcesses$cluster25,
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz clusters 2 and 5 proteins)",
                 processesExplainedList = phosphocluster25APE,
                 obviate.log2FC = TRUE)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster25,
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz clusters 2 and 5 proteins)",
                 processesExplainedList = phosphocluster25APE,
                 obviate.log2FC = TRUE)

subPhospho25APE <- phosphocluster25APE[rownames(allPhosphoProcesses$cluster25[1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/GO/STRINGnets/"

for(i in c(1:10)){
  
  variablePathPart <- paste("phospho", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subPhospho25APE[[i]]$q.value,
            log2FC = subPhospho25APE[[i]]$log2.FC,
            proteinNames = subPhospho25APE[[i]]$protein.code,
            descript = completeDescription,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath,
            colSet = "greyblue",
            is.pval.title = TRUE)
  
}

i <- 1

stringnet(qvalue = subPhospho25APE[[i]]$q.value[subPhospho25APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho25APE[[i]]$log2.FC[subPhospho25APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho25APE[[i]]$protein.code[subPhospho25APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/GO/STRINGnets/process1/phospho1.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 2

stringnet(qvalue = subPhospho25APE[[i]]$q.value[subPhospho25APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho25APE[[i]]$log2.FC[subPhospho25APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho25APE[[i]]$protein.code[subPhospho25APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/GO/STRINGnets/process2/phospho2.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 4

stringnet(qvalue = subPhospho25APE[[i]]$q.value[subPhospho25APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho25APE[[i]]$log2.FC[subPhospho25APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho25APE[[i]]$protein.code[subPhospho25APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/GO/STRINGnets/process4/phospho4.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

# KEGG processes (clusters 2 and 5)

phospho25Data <- matrix(data = abs(rnorm(n = (length(allLog2FCPhospho$cluster25) * 6),
                                         mean = 1,
                                         sd = 1)),
                       nrow = length(allLog2FCPhospho$cluster25),
                       ncol = 6)
phospho25Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(phospho25Data) <- names(allLog2FCPhospho$cluster25)
colnames(phospho25Data) <- phospho25Classes

phospho25KEGG <- keggProcesses(rawData = phospho25Data,
                               classesVector = phospho25Classes,
                               keyType = "UNIPROT",
                               is.pval.and.log = TRUE,
                               pval = allPValsPhospho$cluster25,
                               qval = allPValsPhospho$cluster25,
                               log2FC = allLog2FCPhospho$cluster25)

View(phospho25KEGG$bioProcesses$greater)

setWorkSpace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/KEGG/"
keggNetworks(keggres = phospho25KEGG$bioProcesses,
             foldchanges = phospho25KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setWorkSpace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs,
             is.log2FC.absent = TRUE)

excelKeggRes <- phospho25KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("clusters2.5.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- phospho25KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(phospho25KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(phospho25KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/KEGG/phospho25KEGG.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

phospho25TopProcesses <- phospho25KEGG$hsaList[rownames(phospho25KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = phospho25TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/KEGG/phospho25TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = phospho25TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters25/KEGG/phospho25TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
}

################################################################################

# Clusters 3 and 6

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allProteinsInPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPValsPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allLog2FCPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPhosphoProcesses.rda")

phosphocluster36APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                             descript = completeDescription,
                                             chosenPValue = 0.05,
                                             pValuesVector = allPValsPhospho[["cluster36"]],
                                             logFCVector = allLog2FCPhospho[["cluster36"]])

printTopProcesses(processesMatrix = allPhosphoProcesses$cluster36,
                  processesExplained = phosphocluster36APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/GO/topPhosphoClusters36.xlsx",
                  howMany = 10)

plotProcessBars2(processesArr = allPhosphoProcesses$cluster36,
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz clusters 3 and 6 proteins)",
                 processesExplainedList = phosphocluster36APE,
                 obviate.log2FC = TRUE)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster36,
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz clusters 3 and 6 proteins)",
                 processesExplainedList = phosphocluster36APE,
                 obviate.log2FC = TRUE)

subPhospho36APE <- phosphocluster36APE[rownames(allPhosphoProcesses$cluster36[1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/GO/STRINGnets/"

for(i in c(1:10)){
  
  variablePathPart <- paste("phospho", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subPhospho36APE[[i]]$q.value,
            log2FC = subPhospho36APE[[i]]$log2.FC,
            proteinNames = subPhospho36APE[[i]]$protein.code,
            descript = completeDescription,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath,
            colSet = "greyblue",
            is.pval.title = TRUE)
  
}

i <- 2

stringnet(qvalue = subPhospho36APE[[i]]$q.value[subPhospho36APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho36APE[[i]]$log2.FC[subPhospho36APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho36APE[[i]]$protein.code[subPhospho36APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = FALSE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/GO/STRINGnets/process2/phospho2.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 4

stringnet(qvalue = subPhospho36APE[[i]]$q.value[subPhospho36APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho36APE[[i]]$log2.FC[subPhospho36APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho36APE[[i]]$protein.code[subPhospho36APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/GO/STRINGnets/process4/phospho4.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 8

stringnet(qvalue = subPhospho36APE[[i]]$q.value[subPhospho36APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho36APE[[i]]$log2.FC[subPhospho36APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho36APE[[i]]$protein.code[subPhospho36APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/GO/STRINGnets/process8/phospho8.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

# KEGG processes (clusters 3 and 6)

phospho36Data <- matrix(data = abs(rnorm(n = (length(allLog2FCPhospho$cluster36) * 6),
                                         mean = 1,
                                         sd = 1)),
                        nrow = length(allLog2FCPhospho$cluster36),
                        ncol = 6)
phospho36Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(phospho36Data) <- names(allLog2FCPhospho$cluster36)
colnames(phospho36Data) <- phospho36Classes

phospho36KEGG <- keggProcesses(rawData = phospho36Data,
                               classesVector = phospho36Classes,
                               keyType = "UNIPROT",
                               is.pval.and.log = TRUE,
                               pval = allPValsPhospho$cluster36,
                               qval = allPValsPhospho$cluster36,
                               log2FC = allLog2FCPhospho$cluster36)

View(phospho36KEGG$bioProcesses$greater)

setWorkSpace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/KEGG/"
keggNetworks(keggres = phospho36KEGG$bioProcesses,
             foldchanges = phospho36KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setWorkSpace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs,
             is.log2FC.absent = TRUE)

excelKeggRes <- phospho36KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("clusters3.6.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- phospho36KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(phospho36KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(phospho36KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}


write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/KEGG/phospho36KEGG.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

phospho36TopProcesses <- phospho36KEGG$hsaList[rownames(phospho36KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = phospho36TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/KEGG/phospho36TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = phospho36TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/clusters36/KEGG/phospho36TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
}

################################################################################

# Cluster 4

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/refinedPhosphoMatrix.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allProteinsInPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPValsPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allLog2FCPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allPhosphoProcesses.rda")

phosphocluster4APE <- allProcessesExplained(refinedArr = refinedPhosphoMatrix,
                                            descript = completeDescription,
                                            chosenPValue = 0.05,
                                            pValuesVector = allPValsPhospho[["cluster4"]],
                                            logFCVector = allLog2FCPhospho[["cluster4"]])

printTopProcesses(processesMatrix = allPhosphoProcesses$cluster4,
                  processesExplained = phosphocluster4APE,
                  filePath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/GO/topPhosphoCluster4.xlsx",
                  howMany = 10)

plotProcessBars2(processesArr = allPhosphoProcesses$cluster4,
                 howMany = 10,
                 main2 = "-log10(p- and q-value) of altered GO processes (Mfuzz cluster 4 proteins)",
                 processesExplainedList = phosphocluster4APE,
                 obviate.log2FC = TRUE)

plotProcessBars3(processesArr = allPhosphoProcesses$cluster4,
                 howMany = 10,
                 main2 = "-log10(p-value) of altered GO processes (Mfuzz cluster 4 proteins)",
                 processesExplainedList = phosphocluster4APE,
                 obviate.log2FC = TRUE)

subPhospho4APE <- phosphocluster4APE[rownames(allPhosphoProcesses$cluster4[1:10,])]

constantPathPart <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/GO/STRINGnets/"

for(i in c(1:10)){
  
  variablePathPart <- paste("phospho", as.character(i), "Interactions.xlsx", sep = "")
  completePath <- paste(constantPathPart, variablePathPart, sep = "")
  
  stringnet(qvalue = subPhospho4APE[[i]]$q.value,
            log2FC = subPhospho4APE[[i]]$log2.FC,
            proteinNames = subPhospho4APE[[i]]$protein.code,
            descript = completeDescription,
            clusteringAlgorithm = "fastgreedy",
            isSmallNet = TRUE,
            STRINGVer = "11.5",
            selectedSpecies = 9606, # human
            scoreThres = 400,
            chosenPVal = 1,
            isExcel = TRUE,
            filesPath = completePath,
            colSet = "greyblue",
            is.pval.title = TRUE)
  
}

i <- 2

stringnet(qvalue = subPhospho4APE[[i]]$q.value[subPhospho4APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho4APE[[i]]$log2.FC[subPhospho4APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho4APE[[i]]$protein.code[subPhospho4APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = TRUE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/GO/STRINGnets/process2/process2.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

i <- 5

stringnet(qvalue = subPhospho4APE[[i]]$q.value[subPhospho4APE[[i]]$q.value <= 0.05],
          log2FC = subPhospho4APE[[i]]$log2.FC[subPhospho4APE[[i]]$q.value <= 0.05],
          proteinNames = subPhospho4APE[[i]]$protein.code[subPhospho4APE[[i]]$q.value <= 0.05],
          descript = completeDescription,
          clusteringAlgorithm = "fastgreedy",
          isSmallNet = TRUE,
          STRINGVer = "11.5",
          selectedSpecies = 9606, # human
          scoreThres = 400,
          chosenPVal = 1,
          isExcel = FALSE,
          filesPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/GO/STRINGnets/process5/process5.1Interactions.xlsx",
          colSet = "greyblue",
          is.pval.title = TRUE)

# KEGG processes (cluster 4)

phospho4Data <- matrix(data = abs(rnorm(n = (length(allLog2FCPhospho$cluster4) * 6),
                                        mean = 1,
                                        sd = 1)),
                        nrow = length(allLog2FCPhospho$cluster4),
                        ncol = 6)
phospho4Classes <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(phospho4Data) <- names(allLog2FCPhospho$cluster4)
colnames(phospho4Data) <- phospho4Classes

phospho4KEGG <- keggProcesses(rawData = phospho4Data,
                              classesVector = phospho4Classes,
                              keyType = "UNIPROT",
                              is.pval.and.log = TRUE,
                              pval = allPValsPhospho$cluster4,
                              qval = allPValsPhospho$cluster4,
                              log2FC = allLog2FCPhospho$cluster4)

View(phospho4KEGG$bioProcesses$greater)

setWorkSpace <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/KEGG/"
keggNetworks(keggres = phospho4KEGG$bioProcesses,
             foldchanges = phospho4KEGG$log2FC,
             selectSpecies = "hsa",
             nPathways = 10,
             setWorkSpace = setWorkSpace,
             createDir = FALSE,
             keggSetsHS = kegg.sets.hs,
             is.log2FC.absent = TRUE)

excelKeggRes <- phospho4KEGG$bioProcesses$greater
excelKeggRes <- cbind(excelKeggRes, rep(NA, nrow(excelKeggRes)), rep(NA, nrow(excelKeggRes)))
colnames(excelKeggRes)[c((ncol(excelKeggRes)-1), ncol(excelKeggRes))] <- c("cluster4.proteins", "rate")

for(i in 1:nrow(excelKeggRes)){
  
  itotal <- phospho4KEGG$bioProcesses$greater[i,"set.size"]
  iname <- rownames(phospho4KEGG$bioProcesses$greater)[i]
  isig <- sum(na.omit(phospho4KEGG$hsaList[[iname]][,"is.sig.0.05"]))
  irate <- isig/itotal
  excelKeggRes[i,(ncol(excelKeggRes)-1)] <- isig
  excelKeggRes[i,ncol(excelKeggRes)] <- irate
  
}

write.xlsx2(x = excelKeggRes,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/KEGG/phospho4KEGG.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "allKEGGProcesses",
            append = FALSE)

phospho4TopProcesses <- phospho4KEGG$hsaList[rownames(phospho4KEGG$bioProcesses$greater)[1:10]]
for(i in 1:10){
  
  isheetName <- paste("process", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = phospho4TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/KEGG/phospho4TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = FALSE)
  }else{
    write.xlsx2(x = phospho4TopProcesses[[i]],
                file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/clusters_1by1/cluster4/KEGG/phospho4TopKEGG.xlsx",
                row.names = FALSE,
                col.names = TRUE,
                sheetName = isheetName,
                append = TRUE)
  }
  
}

################################################################################

library(xlsx)

allProteomePath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/9_funcAnalysisMock/allGroupsTogether_GO.xlsx"
allPhosphoPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoProcesses/tables/allPhosphoTogether_GO.xlsx"

allProteomeGO <- read.xlsx2(file = allProteomePath, sheetIndex = 1)
allPhosphoGO <- read.xlsx2(file = allPhosphoPath, sheetIndex = 1)

proteomeNames <- allProteomeGO[,1]
phosphoNames <- allPhosphoGO[,1]

allProteomeGO <- allProteomeGO[,-1]
allProteomeGO <- apply(X = allProteomeGO, MARGIN = 2, FUN = as.numeric)
rownames(allProteomeGO) <- proteomeNames

allPhosphoGO <- allPhosphoGO[,-1]
allPhosphoGO <- apply(X = allPhosphoGO, MARGIN = 2, FUN = as.numeric)
rownames(allPhosphoGO) <- phosphoNames

intersectRownames <- intersect(proteomeNames, phosphoNames)
logProteo <- -log10(allProteomeGO[,"fisher.pvalue.allClusters"])
logPhospho <- -log10(allPhosphoGO[,"fisher.pvalue.all"])

logProteo <- logProteo[intersectRownames]
logPhospho <- logPhospho[intersectRownames]

plot(x = logProteo,
     y = logPhospho,
     main = "Proteome vs phosphoproteome -log10(GO altered processes p-value)",
     xlab = "-log10(proteome GO altered processes p-value)",
     ylab = "-log10(phosphoproteome GO altered processes p-value)",
     xlim = range(pretty(logProteo)),
     ylim = range(pretty(logPhospho)),
     pch = 20,
     cex = 1,
     col = "#377EB8")

pearson1 <- cor.test(x = logProteo,
                     y = logPhospho,
                     method = "pearson",
                     conf.level = 0.95)

lineCoeff1 <- lm(logPhospho ~ logProteo)
abline(lineCoeff1, col = "#984EA3", lwd = 2)

legend(x = 6,
       y = 2.5,
       legend = paste("Pearson correlation coefficient (R): ",
                      as.character(round(pearson1$estimate, 2)),
                      sep = ""))

sigProteo <- rownames(allProteomeGO)[allProteomeGO[,"fisher.pvalue.allClusters"] <= 0.05]
sigPhospho <- rownames(allPhosphoGO)[allPhosphoGO[,"fisher.pvalue.all"] <= 0.05]
sigIntersect <- intersect(sigProteo, sigPhospho)

plot(x = logProteo[sigIntersect],
     y = logPhospho[sigIntersect],
     main = "Proteome vs phosphoproteome -log10(GO significantly altered processes p-value)",
     xlab = "-log10(proteome GO significantly altered processes p-value)",
     ylab = "-log10(phosphoproteome GO significantly altered processes p-value)",
     xlim = range(pretty(c(0, logProteo[sigIntersect]))),
     ylim = range(pretty(c(0, logPhospho[sigIntersect]))),
     pch = 20,
     cex = 1,
     col = "#377EB8")

pearson2 <- cor.test(x = logProteo[sigIntersect],
                     y = logPhospho[sigIntersect],
                     method = "pearson",
                     conf.level = 0.95)

lineCoeff2 <- lm(logPhospho[sigIntersect] ~ logProteo[sigIntersect])
abline(lineCoeff2, col = "#984EA3", lwd = 2)

legend(x = 3.5,
       y = 0.5,
       legend = paste("Pearson correlation coefficient (R): ",
                      as.character(round(pearson2$estimate, 2)),
                      sep = ""))

# Diagrama de Venn con el n?mero de procesos significativamente alterados del
# proteoma y del fosfoproteoma.

library(VennDiagram)

venn.diagram(list("Proteome" = sigProteo,
                  "Phos." = sigPhospho),     
             filename = "sigGO_ProteoVSPhospho.png",
             imagetype = "png",
             fill = c("lightgreen", "lightblue"))

################################################################################

# Comparaci?n entre procesos de prote?mica y fosfoprote?mica: KEGG

# Prote?mica:

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAgg.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDMockSig.rda")

allProteome <- rownames(PDAgg)
sigProteome <- rownames(PDMockSig)

proteomeLog2FC <- rep(0, length(allProteome))
names(proteomeLog2FC) <- allProteome
proteomeLog2FC[sigProteome] <- 1

proteomePVal <- rep(1, length(allProteome))
names(proteomeLog2FC) <- allProteome
proteomePVal[sigProteome] <- 0.01

proteoSampleData <- matrix(data = abs(rnorm(n = (length(proteomeLog2FC) * 6),
                                            mean = 1,
                                            sd = 1)),
                            nrow = length(proteomeLog2FC),
                            ncol = 6)
proteoSampleClasses <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(proteoSampleData) <- names(proteomeLog2FC)
colnames(proteoSampleData) <- proteoSampleClasses

proteoKEGG <- keggProcesses(rawData = proteoSampleData,
                            classesVector = proteoSampleClasses,
                            keyType = "UNIPROT",
                            is.pval.and.log = TRUE,
                            pval = proteomePVal,
                            qval = proteomePVal,
                            log2FC = proteomeLog2FC)

rownames(proteoKEGG$bioProcesses$greater)[1:30]

printKEGGresults(keggProcessesResults = proteoKEGG,
                 wanna.Excel = TRUE,
                 excelPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/comparison/KEGG/proteomeKEGGprocesses.xlsx")

# Fosfoprote?mica:

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/allLog2FCPhospho.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")

combinedLog2FC <- rep(0, length(allLog2FCPhospho[[1]]))
names(combinedLog2FC) <- names(allLog2FCPhospho[[1]])

combinedPVal <- rep(1, length(allLog2FCPhospho[[1]]))
names(combinedPVal) <- names(allLog2FCPhospho[[1]])

allSigPhospho <- unique(unlist(relevantClusters))
combinedLog2FC[allSigPhospho] <- 1
combinedPVal[allSigPhospho] <- 0.01

phosphoSampleData <- matrix(data = abs(rnorm(n = (length(allLog2FCPhospho[[1]]) * 6),
                                             mean = 1,
                                             sd = 1)),
                            nrow = length(allLog2FCPhospho[[1]]),
                            ncol = 6)
phosphoSampleClasses <- c(rep("MOCK.t0", 3), rep("SARS.CoV2.t16", 3))
rownames(phosphoSampleData) <- names(allLog2FCPhospho[[1]])
colnames(phosphoSampleData) <- phosphoSampleClasses

phosphoKEGG <- keggProcesses(rawData = phosphoSampleData,
                             classesVector = phosphoSampleClasses,
                             keyType = "UNIPROT",
                             is.pval.and.log = TRUE,
                             pval = combinedPVal,
                             qval = combinedPVal,
                             log2FC = combinedLog2FC)

rownames(phosphoKEGG$bioProcesses$greater)[1:30]

intersect(rownames(proteoKEGG$bioProcesses$greater)[1:30],
          rownames(phosphoKEGG$bioProcesses$greater)[1:30])

library(ggvenn)

ggvenn(
  list("proteome.KEGG" = rownames(proteoKEGG$bioProcesses$greater)[1:30],
       "phosph.KEGG" = rownames(phosphoKEGG$bioProcesses$greater)[1:30]), 
  fill_color = c("forestgreen", "skyblue"),
  stroke_size = 0.5,
  set_name_size = 4
)

printKEGGresults(keggProcessesResults = phosphoKEGG,
                 wanna.Excel = TRUE,
                 excelPath = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/comparison/KEGG/phosphoKEGGprocesses.xlsx")

################################################################################

# Mfuzz: Tabla

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/sigDesc.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/resPhosphoFuzz.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/completeDescription.rda")

clustersWithRep <- list()
for(i in 1:max(resPhosphoFuzz$cluster)){
  
  peptideNames <- names(resPhosphoFuzz$cluster)[resPhosphoFuzz$cluster == i]
  proteinNames <- sigDesc[peptideNames, "Protein"]
  clustersWithRep[[i]] <- proteinNames
  
}

c2.5 <- c(clustersWithRep[[2]], clustersWithRep[[5]])
c3.6 <- c(clustersWithRep[[3]], clustersWithRep[[6]])

updatedWithRep <- list("1" = clustersWithRep[[1]],
                       "2" = c2.5,
                       "3" = c3.6,
                       "4" = clustersWithRep[[4]])

allSigPhospho <- unique(unlist(updatedWithRep))
mFuzzTable <- matrix(data = NA, ncol = 6, nrow = length(allSigPhospho))
rownames(mFuzzTable) <- allSigPhospho
colnames(mFuzzTable) <- c("all.peptides",
                          "cluster.1.peptides",
                          "clusters.2.5.peptides",
                          "clusters.3.6.peptides",
                          "cluster.4.peptides",
                          "how.many.clusters")

for(i in 1:nrow(mFuzzTable)){
  
  iprot <- rownames(mFuzzTable)[i]
  
  for(j in 1:length(updatedWithRep)){
    
    sumBool <- sum(updatedWithRep[[j]] %in% iprot)
    mFuzzTable[iprot, (1 + j)] <- sumBool
    
  }
  
  sumClusters <- sum(mFuzzTable[iprot,2:(1+j)])
  mFuzzTable[iprot, "all.peptides"] <- sumClusters
  howManyClusters <- sum(mFuzzTable[iprot,2:(1+j)] > 0)
  mFuzzTable[iprot, "how.many.clusters"] <- howManyClusters
  
}

mFuzzTable <- mFuzzTable[order(mFuzzTable[,"how.many.clusters"], decreasing = TRUE),]
v <- c()

for(i in max(mFuzzTable[,"how.many.clusters"]):1){
  
  subTable <- mFuzzTable[mFuzzTable[,"how.many.clusters"] == i,]
  subTable <- subTable[order(subTable[,"all.peptides"], decreasing = TRUE),]
  v <- append(x = v, values = rownames(subTable))
  
}

mFuzzTable <- mFuzzTable[v,]

write.xlsx2(x = cbind.data.frame(completeDescription[rownames(mFuzzTable),],
                                 mFuzzTable),
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/comparison/mFuzzTable.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Hoja1")

library(ggvenn)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/fatima2/phospho/phosphoRDA/relevantClusters.rda")

ggvenn(
  relevantClusters, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
