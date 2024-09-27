# Functions:

multiBox <- function(protValues, protName, groups,
                     pvals, onlySig = TRUE, chosenPVals = FALSE, whatPVals,
                     textCex = 1.5){
  
  protDataF <- data.frame(protValues, groups)
  colnames(protDataF) <- c(protName, "classes")
  
  rights <- 1:(length(unique(groups))) + 1/(length(unique(groups))*4)
  rights[length(rights)] <- length(rights) - 1/(length(rights)*4)
  lefts <- 1:(length(unique(groups))) - 1/(length(unique(groups))*4)
  lefts[1] <- 1 + 1/(length(lefts)*4)
  
  firstPretty <- range(pretty(c(min(protValues), max(protValues))))[2]
  secondPretty <- range(pretty(c(min(protValues), max(protValues))))[1]
  difPretty <- firstPretty - secondPretty
  
  indComb <- list()
  namesIndComb <- c()
  k <- 1
  for(i in 1:length(unique(groups))){
    if(i < length(unique(groups))){
      for(j in (i+1):length(unique(groups))){
        indComb[[k]] <- c(i, j)
        tempPName <- paste(unique(groups)[i], unique(groups)[j], sep = "-")
        namesIndComb <- append(x = namesIndComb, values = tempPName)
        k <- k + 1
      }
    }
  }
  names(indComb) <- namesIndComb
  
  sigHeights <- c()
  textHeights <- c()
  
  if(onlySig){
    whichPVals <- which(pvals <= 0.05)
  }
  if(chosenPVals){
    whichPVals <- whatPVals
  }
  
  p <- 1
  q <- 2
  for(i in 1:length(whichPVals)){
    newSig <- max(protValues) + (p*(difPretty/50))
    sigHeights <- append(x = sigHeights, values = newSig)
    newText <- max(protValues) + (q*(difPretty/50))
    textHeights <- append(x = textHeights, values = newText)
    p <- p + 2
    q <- q + 2
  }
  
  finalPretty <- max(textHeights)
  
  protDataF$classes <- factor(protDataF$classes , levels = unique(groups)) ########### !!!
  boxplot(protDataF[,1] ~ protDataF[,2],
          data = protDataF,
          pch = 21,
          bg = "grey",
          xaxt = "n",
          xlab = "Groups",
          ylab = paste(protName, " normalized levels", sep = ""),
          main = paste(protName, " normalized levels across groups", sep = ""),
          ylim = range(pretty(c(min(protValues), finalPretty))))
  
  axis(side = 1, at = 1:length(unique(groups)), labels = unique(groups))
  
  if(length(whichPVals) >= 1){
    for(i in 1:length(whichPVals)){
      ########################################################################## !!!
      
      if(length(unique(groups)) == 2){
        tempStart <- 1
        tempEnd <- 2
      }
      ti <- whichPVals[i]
      newnames <- gsub(pattern = "adj.pval.", replacement = "", x = names(ti), fixed = TRUE)
      newnames <- gsub(pattern = "pval.", replacement = "", x = newnames, fixed = TRUE)
      group1 <- strsplit(x = newnames, split = "-", fixed = TRUE)[[1]][1]
      group2 <- strsplit(x = newnames, split = "-", fixed = TRUE)[[1]][2]
      orderMe <- c()
      orderMe <- append(x = orderMe, values = which(unique(groups) == group1))
      orderMe <- append(x = orderMe, values = which(unique(groups) == group2))
      orderMe <- orderMe[order(orderMe)]
      
      tempStart <- orderMe[1]
      tempEnd <- orderMe[2]
      ########################################################################## !!!
      # tempStart <- indComb[[ti]][1]
      # tempEnd <- indComb[[ti]][2]
      
      tempStart <- rights[tempStart]
      tempEnd <- lefts[tempEnd]
      tempMean <- (tempStart + tempEnd)/2
      littleSum <- difPretty/100
      segments(x0 = tempStart,
               x1 = tempEnd,
               y0 = sigHeights[i])
      segments(x0 = tempStart,
               y0 = sigHeights[i],
               y1 = sigHeights[i] - littleSum)
      segments(x0 = tempEnd,
               y0 = sigHeights[i],
               y1 = sigHeights[i] - littleSum)
      tempPval <- pvals[ti]
      
      if(tempPval > 0.05){
        tempText <- "ns"
      }
      if(tempPval <= 0.05){
        tempText <- "*"
      }
      if(tempPval <= 0.01){
        tempText <- "**"
      }
      if(tempPval <= 0.001){
        tempText <- "***"
      }
      if(tempPval <= 0.0001){
        tempText <- "****"
      }
      
      text(x = tempMean, y = textHeights[i], labels = tempText, cex = textCex)
    }
  }
  
}

photoBox <- function(is.already = TRUE,
                     nameDir,
                     where,
                     dat,
                     chosenQuality = 100,
                     chosenWidth = 600,
                     chosenHeight = 480,
                     pvalues,
                     groups2,
                     is.ANOVA = TRUE){
  
  if(!(is.already)){
    where <- paste(where, nameDir, "/", sep = "")
    dir.create(path = where)
  }
  
  for(i in 1:nrow(dat)){
    fileName <- paste(rownames(dat)[i], ".jpeg", sep = "")
    filePath <- paste(where, fileName, sep = "")
    jpeg(file = filePath,
         quality = chosenQuality,
         width = chosenWidth,
         height = chosenHeight)
    if(is.ANOVA){
      multiBox(protValues = dat[rownames(dat)[i],],
               protName = rownames(dat)[i],
               groups = groups2,
               pvals = pvalues[i,],
               onlySig = TRUE,
               chosenPVals = FALSE)
    }else{
      multiBox(protValues = dat[rownames(dat)[i],],
               protName = rownames(dat)[i],
               groups = groups2,
               pvals = pvalues[i],
               onlySig = TRUE,
               chosenPVals = FALSE)
    }
    
    dev.off()
  }
  
}

mfuzzing <- function(aggregatedResult,
                     chosenRange = 2:50){
  
  if(!("package:Mfuzz" %in% search())){
    library(Mfuzz)
  } 
  
  allMyResults <- list()
  
  aggSet <- ExpressionSet(assayData = aggregatedResult)
  aggTmp <- filter.std(aggSet, min.std = 0)
  
  originalProteins <- rownames(aggSet)
  filteredProteins <- rownames(aggTmp@assayData$exprs)
  filteredOut <- setdiff(x = originalProteins, y = filteredProteins)
  allMyResults[[1]] <- originalProteins
  allMyResults[[2]] <- filteredProteins
  allMyResults[[3]] <- filteredOut
  
  standAgg <- standardise(aggTmp)
  allMyResults[[4]] <- standAgg
  mAgg <- mestimate(standAgg)
  allMyResults[[5]] <- mAgg
  dAgg <- Dmin(eset = standAgg, m = mAgg, crange = chosenRange)
  allMyResults[[6]] <- dAgg
  
  plot(x = chosenRange,
       y = dAgg,
       pch = 16,
       cex = 0.8,
       col = "black",
       xlab = "Number of clusters",
       ylab = "Minimum distance 'D min' between cluster centroid",
       main = "Minimum distance between cluster centroid vs. Number of clusters",
       xlim = range(pretty(chosenRange)),
       ylim = range(pretty(dAgg)))
  points(x = chosenRange,
         y = dAgg,
         col = "black",
         type = "l",
         cex = 1)
  
  names(allMyResults) <- c("original.features",
                           "filtered.in",
                           "filtered.out",
                           "standarized.data",
                           "m.fuzzification",
                           "d.parameter")
  
  return(allMyResults)
  
}

clAndPlot <- function(mfuzzRes,
                      chosenC,
                      chosenMF = c(3, 2),
                      chosenTimes,
                      xlab2 = "Time (h)"){
  
  if(!("package:Mfuzz" %in% search())){
    library(Mfuzz)
  }
  
  clObj <- mfuzz(eset = mfuzzRes[["standarized.data"]],
                 centers = chosenC,
                 m = mfuzzRes[["m.fuzzification"]])
  mfuzz.plot2(eset = mfuzzRes[["standarized.data"]],
              cl = clObj,
              mfrow = chosenMF,
              time.labels = chosenTimes,
              xlab = xlab2)
  return(clObj)
  
}

library(xlsx)

# Load necessary R data

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normCopy.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/aggRes.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/infoSARS.rda")

timeSARS <- c(0, 3, 6, 9, 16)
mfuzzRAll <- mfuzzing(aggregatedResult = aggRes, chosenRange = 2:50)
save(x = mfuzzRAll,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzRAll.rda")
clRAll <- clAndPlot(mfuzzRes = mfuzzRAll,
                    chosenC = 20,
                    chosenMF = c(3, 2),
                    chosenTimes = timeSARS,
                    xlab2 = "Time (h)")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_All/clustersRAll.xlsx"
for(i in 1:max(clRAll$cluster)){
  
  iprots <- names(clRAll$cluster)[clRAll$cluster == i]
  itable <- infoSARS[iprots,]
  itable <- itable[order(itable[,"pvalue.Ftest"], decreasing = FALSE),]
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

# Else:

partialPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_All/"
for(i in 1:max(clRAll$cluster)){
  
  iprots <- names(clRAll$cluster)[clRAll$cluster == i]
  itable <- infoSARS[iprots,]
  itable <- itable[order(itable[,"pvalue.Ftest"], decreasing = FALSE),]
  icluster <- paste("Cluster", as.character(i), sep = "")
  completePath <- paste(partialPath, icluster, sep = "")
  
  write.xlsx2(x = itable,
              file = completePath,
              sheetName = icluster,
              col.names = TRUE,
              row.names = FALSE,
              append = FALSE)
  
}

pvalMat <- adjANOVASARS[,3:12]

boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_All/clusterBoxplots/"
for(i in 1:max(clRAll$cluster)){
  
  iprots <- names(clRAll$cluster)[clRAll$cluster == i]
  idata <- normCopy[iprots,]
  ipval <- pvalMat[iprots,]
  icluster <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = icluster,
           where = boxDir,
           dat = idata,
           pvalues = ipval,
           groups2 = groups,
           is.ANOVA = TRUE)
  
}

################################################################################

# Significant proteins obtained with R:

sigProt <- rownames(adjANOVASARS)[adjANOVASARS[,"qvalue.Ftest"] <= 0.05]
aggSig <- aggRes[sigProt,]
mfuzzRSig <- mfuzzing(aggregatedResult = aggSig, chosenRange = 2:50)
save(x = mfuzzRSig,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzRSig.rda")
clRSig <- clAndPlot(mfuzzRes = mfuzzRSig,
                    chosenC = 10,
                    chosenMF = c(3, 2),
                    chosenTimes = timeSARS,
                    xlab2 = "Time (h)")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_Sig/clustersRSig.xlsx"
for(i in 1:max(clRSig$cluster)){
  
  iprots <- names(clRSig$cluster)[clRSig$cluster == i]
  itable <- infoSARS[iprots,]
  itable <- itable[order(itable[,"pvalue.Ftest"], decreasing = FALSE),]
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

boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_Sig/clusterBoxplots/"
for(i in 1:max(clRSig$cluster)){
  
  iprots <- names(clRSig$cluster)[clRSig$cluster == i]
  idata <- normCopy[iprots,]
  ipval <- pvalMat[iprots,]
  icluster <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = icluster,
           where = boxDir,
           dat = idata,
           pvalues = ipval,
           groups2 = groups,
           is.ANOVA = TRUE)
  
}

################################################################################

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

#-------------------------------------------------------------------------------

infoSARS2 <- description
mfuzzPDSig <- mfuzzing(aggregatedResult = PDMockSig, chosenRange = 2:10)
save(x = mfuzzPDSig,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzPDSig.rda")
cPDSig <- 4
clPDSig <- clAndPlot(mfuzzRes = mfuzzPDSig,
                     chosenC = cPDSig,
                     chosenMF = c(2, 2),
                     chosenTimes = timeSARS,
                     xlab2 = "Time (h)")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/tTestPD_Sig/clustersPDSig.xlsx"
for(i in 1:max(clPDSig$cluster)){
  
  iprots <- names(clPDSig$cluster)[clPDSig$cluster == i]
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

boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/tTestPD_Sig/clusterBoxplots/"
for(i in 1:max(clPDSig$cluster)){
  
  iprots <- names(clPDSig$cluster)[clPDSig$cluster == i]
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

################################################################################

file2 <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/allTTestPD.xlsx"
PDAbund2 <- read.xlsx2(file = file2, sheetName = "abundances")
tAll <- read.xlsx2(file = file2, sheetName = "ttest")

acc2 <- PDAbund2[,1]
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/PDAgg.rda")
PDAgg2 <- PDAgg[acc2,]
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
acc3 <- rownames(tAll)
tAll <- tAll[,-1]
colnames(tAll) <- colnames(adjANOVASARS[,3:12])
tAll <- apply(X = tAll, MARGIN = 2, FUN = as.numeric)
rownames(tAll) <- acc3

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
infoSARS3 <- cbind.data.frame(description[rownames(tAll),], tAll)

mfuzzPDLast <- mfuzzing(aggregatedResult = PDAgg2, chosenRange = 2:50)
save(x = mfuzzPDLast,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzPDLast.rda")
timeSARS <- c(0, 3, 6, 9, 16)
cPDLast <- 4
clPDLast<- clAndPlot(mfuzzRes = mfuzzPDLast,
                     chosenC = cPDLast,
                     chosenMF = c(2, 2),
                     chosenTimes = timeSARS,
                     xlab2 = "Time (h)")

excelFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/newTTestPD_sig/clustersPD.xlsx"

for(i in 1:max(clPDLast$cluster)){
  
  iprots <- names(clPDLast$cluster)[clPDLast$cluster == i]
  itable <- infoSARS3[iprots,]
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

boxDir <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/newTTestPD_sig/clusterBoxplots/"

for(i in 1:max(clPDLast$cluster)){
  
  iprots <- names(clPDLast$cluster)[clPDLast$cluster == i]
  idata <- log2(PDAbundances[iprots,])
  ipval <- tAll[iprots,]
  icluster <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = icluster,
           where = boxDir,
           dat = idata,
           pvalues = ipval,
           groups2 = groups,
           is.ANOVA = TRUE)
  
}

################################################################################

torf <- normCopy[rownames(adjANOVASARS)[adjANOVASARS[,"qvalue.Ftest"] <= 0.05],]
SARSProteins <- c("P0DTD1",
                  "P0DTC3",
                  "P0DTC7",
                  "P0DTC8",
                  "P0DTC5",
                  "P0DTC2",
                  "P0DTD2",
                  "P0DTC9")

torf <- torf[-which(rownames(torf) %in% SARSProteins),]

library(xlsx)

peptideAtlasList <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/plasma.txt"
peptideAtlasMatrix <- read.delim(file = peptideAtlasList)
plasmaProteins <- peptideAtlasMatrix[,1]
torf <- torf[intersect(rownames(torf), plasmaProteins),]
torf <- t(torf)

save(x = torf,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/torf.rda")

plasmaRF <- recursiveRandomForest3(dataset = torf,
                                   groups = as.factor(groups),
                                   treeNumber = 1000,
                                   retentionProportion = howManyLoops(nfeat = ncol(torf),
                                                                      wantedLoops = 30)[1],
                                   chosenIterationsNumber = 50,
                                   withNormalization = FALSE,
                                   topNumber = 50,
                                   iterationsBest = 1000,
                                   searchForBest = TRUE,
                                   minimization = TRUE)

View(plasmaRF[["importance.matrix"]])
plasmaRF[["chosen.features"]]
plasmaRF[["best.RF"]]

plasmaRF2 <- recursiveRandomForest3(dataset = torf,
                                    groups = as.factor(groups),
                                    treeNumber = 1000,
                                    retentionProportion = howManyLoops(nfeat = ncol(torf),
                                                                       wantedLoops = 30)[1],
                                    chosenIterationsNumber = 50,
                                    withNormalization = FALSE,
                                    topNumber = 50,
                                    iterationsBest = 1000,
                                    searchForBest = FALSE,
                                    howManyFeats = 15)