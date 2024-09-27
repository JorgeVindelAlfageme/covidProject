load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")

library(xlsx)

contaminants <- c("P00761", "P00883", "P02769")
allHumanProts <- setdiff(x = rownames(adjANOVASARS), y = contaminants)
adjANOVASARS <- adjANOVASARS[allHumanProts,]
ANOVASig <- adjANOVASARS[adjANOVASARS[,"qvalue.Ftest"] <= 0.05,]
adjTable <- adjANOVASARS[allHumanProts,]
adjTable <- cbind.data.frame(description[allHumanProts,], adjTable)
aggSig <- aggRes[rownames(ANOVASig),]

mfuzzSig <- mfuzzing(aggregatedResult = aggSig, chosenRange = 2:50)
save(x = mfuzzSig,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzSig.rda")
# load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzSig.rda")
cSig <- 10
clSig <- mfuzz(eset = mfuzzSig[["standarized.data"]],
               c = cSig,
               m = mfuzzSig[["m.fuzzification"]])
mfuzz.plot2(eset = mfuzzSig[["standarized.data"]],
            cl = clSig,
            mfrow = c(3, 2),
            time.labels = c(0, 3, 6, 9, 16),
            xlab = "Time (h)")

clustersSig <- clSig$cluster
excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_sig/clustersSig.xlsx"
for(i in 1:max(clustersSig)){
  
  iprots <- names(clustersSig)[clustersSig == i]
  itable <- adjTable[iprots,]
  itable <- itable[order(itable[,"pvalue.Ftest"], decreasing = FALSE),]
  isheet <- paste("Cluster", as.character(i), sep = "")
  
  if(i == 1){
    write.xlsx2(x = itable,
                file = excelFile,
                sheetName = isheet,
                row.names = FALSE,
                col.names = TRUE,
                append = FALSE)
  }else{
    write.xlsx2(x = itable,
                file = excelFile,
                sheetName = isheet,
                row.names = FALSE,
                col.names = TRUE,
                append = TRUE)
  }
  
  
}

pvalMatrix <- adjANOVASARS[,3:12]

for(i in 1:max(clustersSig)){
  
  iprots <- names(clustersSig)[clustersSig == i]
  idir <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = idir,
           where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_sig/clusterBoxplots/",
           dat = normCopy[iprots,],
           pvalues = pvalMatrix[iprots,],
           groups2 = groups,
           is.ANOVA = TRUE)
  
}

fileLocation <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm/RNorm_All/RNorm_All_Clusters2.xlsx"
for(i in 7:length(mfuzzClusters)){
  isheet <- paste("Cluster", as.character(i), sep = "")
  iprots <- mfuzzClusters[[i]]
  clusterData <- adjTable[iprots,]
  clsterData <- clusterData[order(clusterData[,"pvalue.Ftest"], decreasing = FALSE),]
  if(i == 1){
    write.xlsx2(x = clusterData,
                file = fileLocation,
                sheetName = isheet,
                row.names = FALSE,
                col.names = TRUE,
                append = FALSE)
  }else{
    write.xlsx2(x = clusterData,
                file = fileLocation,
                sheetName = isheet,
                row.names = FALSE,
                col.names = TRUE,
                append = TRUE)
  }
}

sigSARS <- rownames(adjANOVASARS)[adjANOVASARS[,"qvalue.Ftest"] <= 0.05]
sigOnlyHuman <- setdiff(x = sigSARS, nonHuman)

################################################################################


