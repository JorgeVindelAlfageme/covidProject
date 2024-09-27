load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAbundances.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzPDSig.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/ttestPD.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/allTTestPD.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAgg.rda")

# Calcular PDAbundances2 (con la tabla de abundancias de la tabla de los 4
# t-tests).

set.seed(1)

lastData <- PDAbundances[rownames(allTTestPD),]
lastAgg <- PDAgg[rownames(allTTestPD),]
mfuzzLast <- mfuzzing(aggregatedResult = lastAgg, chosenRange = 2:50)
save(x = mfuzzLast,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzLast.rda")

timeSARS <- c(0, 3, 6, 9, 16)
cPDLast <- 4
clPDLast <- clAndPlot(mfuzzRes = mfuzzLast,
                      chosenC = cPDLast,
                      chosenMF = c(2, 2),
                      chosenTimes = timeSARS,
                      xlab2 = "Time (h)")
save(x = clPDLast,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/clPDLast.rda")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/chosenT/all/clusters.xlsx"
for(i in 1:max(clPDLast$cluster)){
  
  iprots <- names(clPDLast$cluster)[clPDLast$cluster == i]
  itable <- cbind.data.frame(description[iprots,], allTTestPD[iprots,])
  icluster <- paste("Cluster", as.character(i), sep = "")
  if(i == 1){
    write.xlsx2(x = itable,
                file = excelFile,
                sheetName = icluster,
                col.names = TRUE,
                row.names = FALSE,
                append = FALSE)
  }else{
    write.xlsx2(x = itable,
                file = excelFile,
                sheetName = icluster,
                col.names = TRUE,
                row.names = FALSE,
                append = TRUE)
  }
}

boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/chosenT/mock/boxplots/"
for(i in 1:max(clPDSig$cluster)){
  
  iprots <- names(clPDSig$cluster)[clPDSig$cluster == i]
  idata <- PDAbundances[iprots,]
  ipval <- ttestPD[iprots,]
  icluster <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = icluster,
           where = boxDir,
           dat = idata,
           pvalues = ipval,
           groups2 = groups,
           is.ANOVA = TRUE)
  
}


