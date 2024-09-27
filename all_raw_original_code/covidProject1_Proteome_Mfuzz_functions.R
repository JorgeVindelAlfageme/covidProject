## retrieveGeneID():

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

## mfuzzing()

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

## clAndPlot():

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

# "UniProtR" was a library upon which the user-defined function "bioProcesses2"
# was defined. It can be possible that this function can't be run due to the
# changes that happened in the library. The result was saved in a variable that
# can be loaded, even if the original function can't be executed.

## bioProcesses2():

bioProcesses2 <- function(uniProtAccessionVector, pvalVec){
  
  if(!("package:UniprotR" %in% base::search())){
    library(UniprotR)
  }
  
  GOlist <- list()
  
  proteinProcesses <- GetProteinGOInfo(ProteinAccList = uniProtAccessionVector,
                                       directorypath = NULL)[,3]
  
  for(i in 1:length(proteinProcesses)){
    GOlist[[i]] <- unlist(strsplit(x = proteinProcesses[i], split = "; "))
  }
  names(GOlist) <- uniProtAccessionVector
  
  processesColumns<- unique(na.omit(unname(unlist(GOlist))))
  
  numberRowsProcesses <- length(uniProtAccessionVector)
  numberColumnsProcesses <- length(processesColumns)
  
  processesMatrix <- matrix(data = rep(0, numberRowsProcesses * numberColumnsProcesses),
                            nrow = numberRowsProcesses, byrow = TRUE)
  rownames(processesMatrix) <- uniProtAccessionVector
  colnames(processesMatrix) <- processesColumns
  
  for(i in 1:length(GOlist)){
    if(pvalVec[i] <= 0.05){
      processesMatrix[i,]  <- as.numeric(processesColumns %in% GOlist[[i]])
    }
  }
  
  processesMatrix <- processesMatrix[rowSums(processesMatrix) > 0,]
  
  return(processesMatrix)
  
}

## newAlteredProcesses():

# 

newAlteredProcesses <- function(refinedArr,
                                proteinVec){
  
  zero <- matrix(data = 0, nrow = nrow(refinedArr), ncol = ncol(refinedArr))
  rownames(zero) <- rownames(refinedArr)
  colnames(zero) <- colnames(refinedArr)
  
  trueProteins <- intersect(proteinVec, rownames(refinedArr))
  zero[trueProteins,] <- refinedArr[trueProteins,]
  
  refinedSum <- apply(X = refinedArr, MARGIN = 2, FUN = sum)
  zeroSum <- apply(X = zero, MARGIN = 2, FUN = sum)
  
  processesResults <- matrix(data = NA, nrow = length(refinedSum), ncol = 5)
  rownames(processesResults) <- names(refinedSum)
  colnames(processesResults) <- c("total.proteins",
                                  "altered.proteins",
                                  "protein.rate",
                                  "fisher.pvalue",
                                  "fisher.qvalue")
  processesResults[,1] <- refinedSum
  processesResults[,2] <- zeroSum
  processesResults[,3] <- zeroSum / refinedSum
  
  totalProteins <- nrow(refinedArr)
  alteredProteins <- length(trueProteins)
  sig <- alteredProteins
  noSig <- totalProteins - sig
  
  fisher.vector <- c()
  
  for(i in 1:length(refinedSum)){
    
    ni <- as.numeric(refinedSum[i])
    isig <- as.numeric(zeroSum[i])
    inosig <- ni - isig
    
    trueNoSig <- noSig - inosig
    trueSig <- sig - isig
    
    fisher.data <- data.frame("non.specific" = c(trueNoSig, trueSig),
                              "specific" = c(inosig, isig))
    rownames(fisher.data) <- c("non.significant", "significant")
    
    ipvalue <- as.numeric(fisher.test(fisher.data, alternative = "greater")$p.value)
    fisher.vector <- append(x = fisher.vector, values = ipvalue)
    
  }
  
  fisher.adjusted <- p.adjust(p = fisher.vector, "fdr")
  
  processesResults[,4] <- fisher.vector
  processesResults[,5] <- fisher.adjusted
  
  processesResults <- processesResults[order(processesResults[,4], decreasing = FALSE),]
  
  return(processesResults)
  
}

## allProcessesExplained():

allProcessesExplained <- function(refinedArr,
                                  descript,
                                  chosenPValue = 0.05,
                                  pValuesVector,
                                  logFCVector){
  
  zero <- matrix(data = 0, nrow = nrow(refinedArr), ncol = ncol(refinedArr))
  rownames(zero) <- rownames(refinedArr)
  colnames(zero) <- colnames(refinedArr)
  selectedProteins <- names(pValuesVector)[pValuesVector <= chosenPValue]
  trueSelectedProteins <- intersect(x = rownames(refinedArr), selectedProteins)
  zero[trueSelectedProteins,] <- refinedArr[trueSelectedProteins,]
  
  matList <- list()
  for(i in 1:ncol(refinedArr)){
    
    icol <- refinedArr[,i]
    iproteins <- names(icol)[icol == 1]
    imat <- data.frame(matrix(data = NA, nrow = length(iproteins), ncol = 8))
    imat[,1] <- colnames(refinedArr)[i]
    imat[,2] <- iproteins
    imat[,3] <- descript[iproteins,"description"]
    imat[,4] <- descript[iproteins,"gene.ID"]
    imat[,5] <- pValuesVector[iproteins]
    imat[,6] <- pValuesVector[iproteins] <= 0.05
    imat[,7] <- logFCVector[iproteins]
    isUpOrDown <- logFCVector[iproteins] >= 0
    isUpOrDown[isUpOrDown == TRUE] <- "Up"
    isUpOrDown[isUpOrDown == "FALSE"] <- "Down"
    imat[,8] <- isUpOrDown
    colnames(imat) <- c("biological.process",
                        "protein.code",
                        "description",
                        "gene.ID",
                        "q.value",
                        "is.significant.0.05",
                        "log2.FC",
                        "is.Up.Or.Down")
    imat <- imat[order(imat[,5], decreasing = FALSE),]
    matList[[i]] <- imat
    
  }
  names(matList) <- colnames(refinedArr)
  
  return(matList)
  
}