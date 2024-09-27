# 1) Loading packages, data, and separating proteomics values from protein
# descriptions.

# sarsPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/lessData.xlsx"
sarsPath <- "C:/Users/jvalf/Desktop/Archivos/trabajo2/partOfFornax/fatima/lessData.xlsx"
library(xlsx)
library(limma)
sarsData <- read.xlsx2(file = sarsPath, sheetIndex = 1)
rownames(sarsData) <- sarsData[,1]
description <- sarsData[, c(1,2)]
onlyDescription <- description[, 2]

# 2) Adding new data to description matrix. 

retrieveGeneID <- function(desc){
  geneID <- c()
  for(string in onlyDescription){
    re <- regexpr(pattern = "GN=\\S* ", text = string)
    geneCode <- substr(x = string,
                       start = re[1] + 3,
                       stop = re[1] + attr(x = re, which = "match.length") - 2)
    if(geneCode == "" | is.na(geneCode)){
      geneCode <- ""
    }
    geneID <- append(x = geneID, values = geneCode)
  }
  return(geneID)
}

description[, 3] <- retrieveGeneID(description[,2])
colnames(description)[3] <- "gene.ID"
anyNA(description[,3])

# 3) Deleting protein names and descriptions from the originally loaded data
# matrix, shortening its column names and reordering these columns as the user
# pleased.

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

save(x = description,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/descriptionSARS.rda")
save(x = sarsData,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/sarsData.rda")

################################################################################

# 4) Normalizing data using the quantile data standardization (not recommended). 

sarsData[sarsData == 0] <- 1
normSARS <- normalizeBetweenArrays(object = log2(na.omit(sarsData)), method = "quantile")

save(x = normSARS,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/normSARS.rda")

################################################################################

# 5) Generating a vector of the different classes the data set represents.

uniqueGroups <- c("MOCK.t0",
                  "SARS.CoV2.t3",
                  "SARS.CoV2.t6",
                  "SARS.CoV2.t9",
                  "SARS.CoV2.t16")
groups <- rep(x = uniqueGroups, each = 3)

save(x = groups,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/SARSgroups.rda")

################################################################################

# 6) Calculating PCA, after altering some proteomics values to make them
# slightly different between each other because they are theoretically
# identical, which difficults some posterior calculations.

sarsPCA <- PCACal(scaledData = normSARS)

set.seed(1)

normCopy <- normSARS

normCopy <- justChangeIt(dataset = normSARS, classes = groups)

save(x = normCopy,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normCopy.rda")

sarsPCA <- PCACal(scaledData = normCopy)
library(RColorBrewer)
SARSColors <- brewer.pal(n = length(unique(groups)), name = "Set1")

# PCAPlot(PCAObject = sarsPCA,
#         is.chosen = TRUE,
#         firstPC = 1,
#         secPC = 2,
#         groups = unique(groups),
#         colVec = SARSColors,
#         legPos = "bottomright",
#         main2 = "Mock infection and SARS infected samples over time PCA")

PCAPlot(PCAObject = sarsPCA,
        is.chosen = TRUE,
        firstPC = 1,
        secPC = 2,
        groups = unique(groups),
        colVec = SARSColors,
        legPos = "bottomright",
        main2 = "Mock infection and SARS infected samples over time PCA",
        cexLeg = 1.5,
        cexPoints = 1.75,
        cexMain = 2,
        cexLabels = 1.55,
        cexAxes = 1.75)

# 7) Calculating the results of k-means on these data and representing its
# results onto PCA plot.

# kmcluster(data = normCopy,
#           chCol = brewer.pal(n = 8, name = "Set1")[6:8],
#           legPos = "bottomright")

kmcluster(data = normCopy,
          maxClus = 10,
          maxIter = 1000,
          chCol = c("lawngreen", "aquamarine4", "skyblue",
                    "palevioletred2", "lemonchiffon4", "orange", "brown2",
                    "darkorchid", "yellow", "navyblue"),
          main2 = "default",
          legPos = "bottomright",
          is.chosen = TRUE,
          firstPC = 1,
          secPC = 2,
          wanna.scale2 = FALSE,
          cexLeg = 1.75,
          cexPoints = 1.75,
          cexMain = 2,
          cexLabels = 1.55,
          cexAxes = 1.75)

# 8) Heat map representation of the similarity of protein expression among
# samples.

savecols <- colnames(normCopy)
heatcols <- gsub(pattern = "SARS.CoV2",
                 replacement = "CoV",
                 x = savecols,
                 fixed = TRUE)
colnames(normCopy) <- heatcols

# heat(dataF = normCopy,
#      xcex = 0.7,
#      colD = FALSE,
#      ylab = FALSE)

heat(dataF = normCopy,
     distance = "correlation",
     xcex = 0.75,
     ycex = 0.75,
     ylab = FALSE,
     rowD = TRUE,
     colD = FALSE,
     colSet = "greenred",
     rotateXLabs = NULL,
     rotateYLabs = NULL,
     ADJCOL = c(NA, 0),
     ADJROW = c(0, NA))

colnames(normCopy) <- savecols

################################################################################

# 9) ANOVA calculations.

ANOVASARS <- multiANOVA3(dat = normCopy, classes = as.factor(groups))
ANOVASARS <- ANOVASARS[order(ANOVASARS[,2], decreasing = FALSE),]

save(x = ANOVASARS,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/ANOVASARS.rda")

difProtSARS <- rownames(ANOVASARS)[ANOVASARS[,2] <= 0.05]
ANOVAout <- cbind.data.frame(description[rownames(ANOVASARS),], ANOVASARS)

library(xlsx)

write.xlsx2(x = ANOVAout,
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/ANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "allProteins")

adjANOVASARS <- adjANOVA(ANOVAres = ANOVASARS)

save(x = adjANOVASARS,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/adjANOVASARS.rda")

adjANOVAout <- cbind.data.frame(description[rownames(adjANOVASARS),], adjANOVASARS)

write.xlsx2(x = adjANOVAout,
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/adjANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "allProteins")

# [1) Cleaning the table from those proteins that are from the virus,
# 2) Choosing only those proteins that have an significant adjusted post-hoc
# p-value, which reduces the number of significant proteins].

humanInd <- grep("OS=Homo sapiens", ANOVAout$Description)
SARSInd <- setdiff(1:nrow(ANOVAout), humanInd) # Not only does it actually
# contain proteins from SARS-CoV2, but from animals such as cow, pig and rabbit
# (these are contaminants from the proteomics experiment method).

# ANOVAout[SARSInd, ]

# SARSInd <- grep("OS=Severe acute respiratory syndrome coronavirus 2",
#                 ANOVAout$Description)

write.xlsx2(x = ANOVAout[humanInd,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/ANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "humanProts",
            append = TRUE)

write.xlsx2(x = ANOVAout[SARSInd,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/ANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "SARSProts",
            append = TRUE)

write.xlsx2(x = adjANOVAout[humanInd,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/adjANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "humanProts",
            append = TRUE)

write.xlsx2(x = adjANOVAout[SARSInd,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/adjANOVASARS.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "SARSProts",
            append = TRUE)

pvalsANOVA <- adjANOVASARS[humanInd, 3:12]

adjANOVAProts <- rescueAdj(pvalsANOVA) # checks if any comparison that
# generated a p-value resulted in a significant p-value. If no p-values from
# a row are significant, that protein would not be considered further on.

write.xlsx2(x = adjANOVAout[adjANOVAProts,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/enrichmentProts.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "humanProts",
            append = FALSE)

pvalsANOVA <- adjANOVASARS[SARSInd, 3:12]
adjANOVAProts <- rescueAdj(pvalsANOVA)

write.xlsx2(x = adjANOVAout[adjANOVAProts,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/enrichmentProts.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "SARSProts",
            append = TRUE)

pvalsANOVA <- adjANOVASARS[humanInd, 3:12]
adjANOVAProts <- rescueAdj(pvalsANOVA)

subsetPvals <- pvalsANOVA[adjANOVAProts, 1:4]

infectionSpecific <- rescueAdj(pvalsArr = subsetPvals)

write.xlsx2(x = adjANOVAout[infectionSpecific,],
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/stateSpecific.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "infectionSpecific",
            append = FALSE)

subsetPvals <- adjANOVASARS[infectionSpecific, 7:12]
uniqueGroups2 <- uniqueGroups[-1]
listT <- list()

for(i in 1:length(uniqueGroups2)){
  
  ti <- grep(pattern = uniqueGroups2[i], x = colnames(subsetPvals))
  arri <- subsetPvals[,ti]
  listT[[i]] <- specificAdj(pvalsArr = arri, zerosAndOnes = c(1,1,1))

}

names(listT) <- uniqueGroups2

View(subsetPvals[, c(2, 4, 6)])

save(x = listT,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/listT.rda")

loopNames <- c("t3.VS.tAny.Specific",
               "t6.VS.tAny.Specific",
               "t9.VS.tAny.Specific",
               "t16.VS.tAny.Specific")

for(i in 1:length(listT)){
  
  write.xlsx2(x = adjANOVAout[listT[[i]],],
              file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/2_statisticalAnalysis/stateSpecific.xlsx",
              row.names = FALSE,
              col.names = TRUE,
              sheetName = loopNames[i],
              append = TRUE)

}

worthMeasuring <- Reduce(intersect, listT)
mockVSt3 <- pvalsANOVA[,"adj.pval.SARS.CoV2.t3-MOCK.t0"]
worthMeasuring <- append(x = worthMeasuring,
                         values = names(mockVSt3[mockVSt3 <= 0.05]))
worthData <- normCopy[worthMeasuring,]

save(x = worthData,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/worthData.rda")

# 10) Using ANOVA results to produce curated figures.

bestData <- normCopy[adjANOVAProts,]
saveCols <- colnames(bestData)
colnames(bestData) <- gsub(pattern = "SARS.CoV2",
                           replacement = "CoV",
                           x = colnames(bestData))
heat(dataF = bestData, ylab = FALSE, colD = FALSE, xcex = 0.7)
colnames(bestData) <- saveCols

PCAPlot(PCAObject = PCACal(scaledData = bestData),
        is.chosen = TRUE,
        firstPC = 1,
        secPC = 2,
        groups = uniqueGroups,
        colVec = SARSColors,
        legPos = "bottomright",
        main2 = "SARS-CoV2 infection samples PCA with significant proteins selected from ANOVA",
        cexLeg = 1.5,
        cexPoints = 1.75,
        cexMain = 1.45,
        cexLabels = 1.5,
        cexAxes = 1.5)

kmcluster(data = bestData,
          chCol = RColorBrewer::brewer.pal(n = 8, name = "Set1"),
          legPos = "bottomright",
          maxClus = 10,
          maxIter = 1000,
          main2 = "SARS-CoV2 infection samples PCA with significant proteins selected from ANOVA (k-means clustering)",
          is.chosen = TRUE,
          firstPC = 1,
          secPC = 2,
          wanna.scale2 = FALSE,
          cexLeg = 1.5,
          cexPoints = 1.75,
          cexMain = 1.1,
          cexLabels = 1.5,
          cexAxes = 1.5)

# 11) Random forest classifier building

save(x = bestData,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/bestData.rda")

# rrfSARS <- recursiveRandomForest(dataset = t(bestData),
#                                  groups = as.factor(groups),
#                                  treeNumber = 1000,
#                                  retentionProportion = 0.81,
#                                  chosenIterationsNumber = 50,
#                                  withNormalization = FALSE,
#                                  topNumber = 50,
#                                  iterationsBest = 1000,
#                                  searchForBest = TRUE,
#                                 minimization = TRUE)

set.seed(1)
rrfSARS <- recursiveRandomForest3(dataset = t(bestData),
                                  groups = as.factor(groups),
                                  treeNumber = 1000,
                                  retentionProportion = 0.81,
                                  chosenIterationsNumber = 50,
                                  withNormalization = FALSE,
                                  topNumber = 50,
                                  iterationsBest = 1000,
                                  searchForBest = TRUE,
                                  howManyFeats = 8,
                                  minimization = TRUE)

save(x = rrfSARS,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/rrfSARS.rda")

head(rrfSARS$importance.matrix, n = 20)
rrfSARS$chosen.features
rrfSARS$best.RF

set.seed(1)
rfSARS <- randomForest(x = t(worthData),
                       y = as.factor(groups),
                       importance = TRUE,
                       proximity = TRUE)

save(x = rfSARS,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/rfSARS.rda")

# MDS2(rfObj = rfSARS,
#      classCol = groups,
#      main2 = " 2D Multidimensional scaling (MDS) plot after using random forest",
#      legPos = "topleft")

MDS2(rfObj = rfSARS,
     classCol = groups,
     areChosenColors = TRUE,
     chosenColors = RColorBrewer::brewer.pal(n = length(unique(groups)),
                                             name = "Set1"),
     main2 = "Multidimensional scaling (MDS) plot (chosen proteins from covid study into random forest)",
     legPos = "topleft",
     pointCex = 1.75,
     legSize = 1.5,
     cexMain = 1.4,
     cexLabels = 1.5,
     cexAxes = 1.5)

savecols <- colnames(worthData)
colnames(worthData) <- gsub(pattern = "SARS.CoV2",
                            replacement = "CoV",
                            x = colnames(worthData),
                            fixed = TRUE)
heat(dataF = worthData, xcex = 0.7, ycex = 1, colD = FALSE)
colnames(worthData) <- savecols

worthOut <- adjANOVASARS[rownames(worthData),]
worthOut <- worthOut[order(worthOut[,2], decreasing = FALSE),]
write.xlsx2(x = data.frame(description[rownames(worthOut),], worthOut),
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/4_randomForest/RFproteins.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Hoja1",
            append = TRUE)

# multiROC(groupsReal = groups,
#          groupsPredicted = rfSARS$predicted,
#          predictionProb = rfSARS$votes)

multiROC2(realClasses = groups,
          probabilities = rfSARS$votes)

# 12) Other figures generated from covid study data

pvalsMat <- adjANOVASARS[, 3:12]

# c("t16-MOCK", "t3-MOCK", "t6-MOCK", "t9-MOCK", "t3-t16",
# "t6-t16", "t9-t16", "t6-t3", "t9-t3", "t9-t6")

t6Peaks <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(0,0,1,0,0,1,0,1,0,1))
multiBox(protValues = bestData[t6Peaks,],
         protName = t6Peaks,
         groups = groups,
         pvals = pvalsMat[t6Peaks,],
         onlySig = TRUE)

t9Peaks <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(0,0,0,1,0,0,1,0,1,1))
multiBox(protValues = bestData[t9Peaks,],
         protName = t9Peaks,
         groups = groups,
         pvals = pvalsMat[t9Peaks,],
         onlySig = TRUE)

t16Peaks <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(1,0,0,0,1,1,1,0,0,0))
photoBox(is.already = FALSE,
         nameDir = "t16Peaks",
         where = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/5_furtherStatistics/",
         dat = bestData[t16Peaks,],
         pvalues = pvalsMat[t16Peaks,],
         groups2 = groups,
         is.ANOVA = TRUE)

t3Peaks <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(0,1,0,0,1,0,0,1,1,0))
multiBox(protValues = bestData[t3Peaks,],
         protName = t3Peaks,
         groups = groups,
         pvals = pvalsMat[t3Peaks,],
         onlySig = TRUE)

mockPeaks1 <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(1,0,1,1,0,0,0,0,0,0))
mockPeaks2 <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(1,0,0,1,0,0,0,0,0,0))
mockPeaks3 <- specificAdj(pvalsArr = pvalsMat, zerosAndOnes = c(1,0,0,0,0,0,0,0,0,0))
mockPeaks <- c(mockPeaks1, mockPeaks2, mockPeaks3)
photoBox(is.already = FALSE,
         nameDir = "mockPeaks",
         where = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/5_furtherStatistics/",
         dat = bestData[mockPeaks,],
         pvalues = pvalsMat[mockPeaks,],
         groups2 = groups,
         is.ANOVA = TRUE)

mockData <- adjANOVASARS[mockPeaks,]
mockData <- mockData[order(mockData[,2], decreasing = FALSE),]

write.xlsx2(x = cbind.data.frame(description[rownames(mockData),], mockData),
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/5_furtherStatistics/peaks.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "mock",
            append = FALSE)
write.xlsx2(x = data.frame(description[t16Peaks,], adjANOVASARS[t16Peaks,]),
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/5_furtherStatistics/peaks.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "t16",
            append = TRUE)

# 13) Analyzing altered processes considering significant proteins

allProteins <- rownames(normCopy)
pvec <- rep(0.01, length(allProteins))

# referenceMatrix <- bioProcesses(uniProtAccessionVector = allProteins,
#                                 pvalVec = pvec)

referenceMatrix <- bioProcesses2(uniProtAccessionVector = allProteins,
                                 pvalVec = pvec)

t16Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,3] <= 0.05]
t9Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,6] <= 0.05]
t6Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,5] <= 0.05]
t3Proteins <- rownames(adjANOVASARS)[adjANOVASARS[,4] <= 0.05]

refinedMatrix <- referenceMatrix[rowSums(referenceMatrix) > 0,]

save(x = refinedMatrix,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/refinedMatrix.rda")

truet16Proteins <- intersect(x = rownames(refinedMatrix), t16Proteins)

library(xlsx)

load(file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/refinedMatrix.rda")
load(file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/adjANOVASARS.rda")

strangeProteins <- c("P0DTD1",
                     "P0DTC3",
                     "P0DTC7",
                     "P0DTC8",
                     "P0DTC5",
                     "P0DTC2",
                     "P0DTD2",
                     "P0DTC9",
                     "P00761",
                     "P00883",
                     "P02769")

humanProteins <- setdiff(rownames(adjANOVASARS), strangeProteins)
# humanProteins <- setdiff(rownames(refinedMatrix), strangeProteins)
refinedMatrix <- refinedMatrix[humanProteins,]
description <- description[humanProteins,]
adjANOVASARS <- adjANOVASARS[humanProteins,]

t6Proteins <- names(adjANOVASARS[,5])[adjANOVASARS[,5] <= 0.05]
t9Proteins <- names(adjANOVASARS[,6])[adjANOVASARS[,6] <= 0.05]
t16Proteins <- names(adjANOVASARS[,3])[adjANOVASARS[,3] <= 0.05]

t6Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                   proteinVec = t6Proteins)
t9Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                    proteinVec = t9Proteins)
t16Processes <- newAlteredProcesses(refinedArr = refinedMatrix,
                                    proteinVec = t16Proteins)

write.xlsx2(x = t6Processes,
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/processes.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "t6Processes",
            append = FALSE)
write.xlsx2(x = t9Processes,
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/processes.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "t9Processes",
            append = TRUE)
write.xlsx2(x = t16Processes,
            file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/processes.xlsx",
            row.names = TRUE,
            col.names = TRUE,
            sheetName = "t16Processes",
            append = TRUE)

t6ProcessesExplained <- allProcessesExplained(refinedArr = refinedMatrix,
                                              descript = description,
                                              chosenPValue = 0.05,
                                              pValuesVector = adjANOVASARS[,5],
                                              logFCVector = -(adjANOVASARS[,15]))

t9ProcessesExplained <- allProcessesExplained(refinedArr = refinedMatrix,
                                              descript = description,
                                              chosenPValue = 0.05,
                                              pValuesVector = adjANOVASARS[,6],
                                              logFCVector = -(adjANOVASARS[,16]))

t16ProcessesExplained <- allProcessesExplained(refinedArr = refinedMatrix,
                                               descript = description,
                                               chosenPValue = 0.05,
                                               pValuesVector = adjANOVASARS[,3],
                                               logFCVector = -(adjANOVASARS[,13]))

allProcessesExplainedList <- list(t6ProcessesExplained,
                                  t9ProcessesExplained,
                                  t16ProcessesExplained)
names(allProcessesExplainedList) <- c("t6", "t9", "t16")
save(x = allProcessesExplainedList,
     file = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/rda/allProcessesExplainedList.rda")

printTopProcesses(processesMatrix = t6Processes,
                  processesExplained = t6ProcessesExplained,
                  filePath = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/topt6.xlsx",
                  howMany = 10)
printTopProcesses(processesMatrix = t9Processes,
                  processesExplained = t9ProcessesExplained,
                  filePath = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/topt9.xlsx",
                  howMany = 10)
printTopProcesses(processesMatrix = t16Processes,
                  processesExplained = t16ProcessesExplained,
                  filePath = "C:/Users/jvalf/Desktop/Archivos/Trabajo/AlYFer/datos/fatima/6_alteredProcesses/topt16.xlsx",
                  howMany = 10)

plotProcessBars(processesArr = t6Processes,
                howMany = 10,
                main2 = "-log10(p- and q-values) of altered GO processes (t6)")
plotProcessBars(processesArr = t9Processes,
                howMany = 10,
                main2 = "-log10(p- and q-values) of altered GO processes (t9)")
plotProcessBars(processesArr = t16Processes,
                howMany = 10,
                main2 = "-log10(p- and q-values) of altered GO processes (t16)")

################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/allProcesses.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/allProcessesExplainedList.rda")

t6Processes <- allProcesses$t6
t9Processes <- allProcesses$t9
t16Processes <- allProcesses$t16

t6List <- allProcessesExplainedList$t6
t9List <- allProcessesExplainedList$t9
t16List <- allProcessesExplainedList$t16

plotProcessBars2(processesArr = t6Processes,
                 howMany = 10,
                 main2 = "-log10(p- and q-values) of altered GO processes (t6)",
                 processesExplainedList = t6List)

plotProcessBars2(processesArr = t9Processes,
                 howMany = 10,
                 main2 = "-log10(p- and q-values) of altered GO processes (t9)",
                 processesExplainedList = t9List)

plotProcessBars2(processesArr = t16Processes,
                 howMany = 10,
                 main2 = "-log10(p- and q-values) of altered GO processes (t16)",
                 processesExplainedList = t16List)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/bestData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
adjANOVASARS <- adjANOVASARS[rownames(bestData),]

pvalt3 <- adjANOVASARS[,4]
pvalt6 <- adjANOVASARS[,5]
pvalt9 <- adjANOVASARS[,6]
pvalt16 <- adjANOVASARS[,3]

sum(pvalt3 <= 0.05)
sum(pvalt6 <= 0.05)
sum(pvalt9 <= 0.05)
sum(pvalt16 <= 0.05)

library("ggvenn")
library(xlsx)
excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/2_statisticalAnalysis/stateSpecific.xlsx"
t3Specific <- read.xlsx2(file = excelFile, sheetName = "t3.VS.tAny.Specific")[,1]
t6Specific <- read.xlsx2(file = excelFile, sheetName = "t6.VS.tAny.Specific")[,1]
t9Specific <- read.xlsx2(file = excelFile, sheetName = "t9.VS.tAny.Specific")[,1]
t16Specific <- read.xlsx2(file = excelFile, sheetName = "t16.VS.tAny.Specific")[,1]

ggvenn(list("t3" = t3Specific,
            "t6" = t6Specific,
            "t9" = t9Specific,
            "t16" = t16Specific))

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/worthData.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
View(description[rownames(worthData),])

photoBox(is.already = FALSE,
         nameDir = "box", 
         where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/4_randomForest/",
         dat = worthData,
         pvalues = adjANOVASARS[rownames(worthData),
                                grep(pattern = "adj.pval.",
                                     x = colnames(adjANOVASARS),
                                     fixed = TRUE)],
         groups2 = groups,
         is.ANOVA = TRUE)

View(description["Q14CX7",])
