library(xlsx)
PDFile <- "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_2/input/PD_mock.xlsx"
PDAbundances <- read.xlsx2(file = PDFile, sheetName = "abundances")
acc <- PDAbundances[,1]
PDAbundances <- PDAbundances[,-1]
PDAbundances <- apply(FUN = as.numeric, X = PDAbundances, MARGIN = 2)
rownames(PDAbundances) <- acc
PDAbundances <- PDAbundances[,-ncol(PDAbundances)]
PDAgg <- cbind(apply(X = PDAbundances[,1:3], MARGIN = 1, FUN = sum)/3,
               apply(X = PDAbundances[,4:6], MARGIN = 1, FUN = sum)/3,
               apply(X = PDAbundances[,7:9], MARGIN = 1, FUN = sum)/3,
               apply(X = PDAbundances[,10:12], MARGIN = 1, FUN = sum)/3,
               apply(X = PDAbundances[,13:15], MARGIN = 1, FUN = sum)/3)
rownames(PDAgg) <- acc
colnames(PDAgg) <- c("MOCK.t0",
                     "SARS.CoV2.t3",
                     "SARS.CoV2.t6",
                     "SARS.CoV2.t9",
                     "SARS.CoV2.t16")
PDAgg <- log2(PDAgg)
contaminants <- c("P00761", "P00883", "P02769")
PDAgg <- PDAgg[-which(rownames(PDAgg) %in% contaminants),]

save(x = PDAgg,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAgg.rda")

###

sarsPath <- "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/input/lessData.xlsx"
sarsData <- read.xlsx2(file = sarsPath, sheetIndex = 1)
rownames(sarsData) <- sarsData[, 1]
description <- sarsData[, c(1, 2)]
onlyDescription <- description[, 2]

# 2) Adding new data to description matrix

description[, 3] <- retrieveGeneID(description[,2])
colnames(description)[3] <- "gene.ID"

# 3) Generating a vector of the different classes the data set represents

uniqueGroups <- c("MOCK.t0",
                  "SARS.CoV2.t3",
                  "SARS.CoV2.t6",
                  "SARS.CoV2.t9",
                  "SARS.CoV2.t16")
groups <- rep(x = uniqueGroups, each = 3)

###

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")

# All significant proteins for t-test comparisons considering MOCK.t0 group:
mockT3 <- read.xlsx2(file = PDFile, sheetName = "t3")
mockT6 <- read.xlsx2(file = PDFile, sheetName = "t6")
mockT9 <- read.xlsx2(file = PDFile, sheetName = "t9")
mockT16 <- read.xlsx2(file = PDFile, sheetName = "t16")
rownames(mockT3) <- mockT3[,1]
rownames(mockT6) <- mockT6[,1]
rownames(mockT9) <- mockT9[,1]
rownames(mockT16) <- mockT16[,1]

PDMockSig <- unique(c(mockT3[,1], mockT6[,1], mockT9[,1], mockT16[,1]))
PDMockSig <- PDAgg[PDMockSig,]
save(x = PDMockSig,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDMockSig.rda")

ttestPD <- matrix(
  data = NA,
  nrow = 6709, # Number of original proteins in raw data tables.
  ncol = 10 # Number of possible paired comparisons between groups.
)
# ttestPD <- matrix(data = NA, nrow = nrow(pvalMat), ncol = ncol(pvalMat))

rownames(ttestPD) <- rownames(PDAgg)
colnames(ttestPD) <- c(
  "adj.pval.SARS.CoV2.t16-MOCK.t0", "adj.pval.SARS.CoV2.t3-MOCK.t0",
  "adj.pval.SARS.CoV2.t6-MOCK.t0", "adj.pval.SARS.CoV2.t9-MOCK.t0",
  "adj.pval.SARS.CoV2.t3-SARS.CoV2.t16", "adj.pval.SARS.CoV2.t6-SARS.CoV2.t16",
  "adj.pval.SARS.CoV2.t9-SARS.CoV2.t16", "adj.pval.SARS.CoV2.t6-SARS.CoV2.t3",
  "adj.pval.SARS.CoV2.t9-SARS.CoV2.t3", "adj.pval.SARS.CoV2.t9-SARS.CoV2.t6"
)

# t3 -> 2
# t16 -> 1
# t6 -> 3
# t9 -> 4

for(i in 1:nrow(ttestPD)){
  iprot <- rownames(ttestPD)[i]
  if(!(is.na(mockT3[iprot,2])) & mockT3[iprot,2] <= 0.05){
    ttestPD[i,2] <- mockT3[iprot,2]
  }
  if(!(is.na(mockT6[iprot,2])) & mockT6[iprot,2] <= 0.05){
    ttestPD[i,3] <- mockT6[iprot,2]
  }
  if(!(is.na(mockT9[iprot,2])) & mockT9[iprot,2] <= 0.05){
    ttestPD[i,4] <- mockT9[iprot,2]
  }
  if(!(is.na(mockT16[iprot,2])) & mockT16[iprot,2] <= 0.05){
    ttestPD[i,1] <- mockT16[iprot,2]
  }
}
ttestPD[is.na(ttestPD)] <- 1

save(x = ttestPD,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/ttestPD.rda")

# Let's start

infoSARS2 <- cbind.data.frame(description[rownames(ttestPD),], ttestPD)

timeSARS <- c(0, 3, 6, 9, 16)

mfuzzPDAll <- mfuzzing(aggregatedResult = PDAgg, chosenRange = 2:10)
save(x = mfuzzPDAll,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzPDAll.rda")
cPDAll <- 4
clPDAll <- clAndPlot(mfuzzRes = mfuzzPDAll,
                     chosenC = cPDAll,
                     chosenMF = c(3, 2),
                     chosenTimes = timeSARS,
                     xlab2 = "Time (h)")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/tTestPD_All/clustersPDAll.xlsx"
for(i in 1:max(clPDAll$cluster)){
  
  iprots <- names(clPDAll$cluster)[clPDAll$cluster == i]
  itable <- infoSARS2[iprots,]
  # itable <- itable[order(itable[,"adj.pval.MOCK.t0-SARS.CoV2.t3"], decreasing = FALSE),]
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

PDAbundances <- PDAbundances[,-ncol(PDAbundances)]
boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/tTestPD_All/clusterBoxplots/"
for(i in 1:max(clPDAll$cluster)){
  
  iprots <- names(clPDAll$cluster)[clPDAll$cluster == i]
  idata <- log2(PDAbundances[iprots,])
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
