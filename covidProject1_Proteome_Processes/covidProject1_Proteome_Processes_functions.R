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

## multiANOVA3():

multiANOVA3 <- function(dat, classes, adjMet = "BH", posthoc = "tukey"){
  
  if(!("package:DescTools" %in% base::search())){
    library(DescTools)
  }
  
  n <- nlevels(classes)
  nrowShellp <- nrow(dat)
  ncolShellp <- 1 + (n * ((n - 1)/2))
  ANOVAshellp <- matrix(data = rep(NA, (nrowShellp * ncolShellp)),
                        nrow = nrowShellp)
  rownames(ANOVAshellp) <- rownames(dat)
  
  nrowShellFC <- nrow(dat)
  ncolShellFC <- n * (n - 1)
  ANOVAshellFC <- matrix(data = rep(NA, (nrowShellFC * ncolShellFC)),
                         nrow = nrowShellFC)
  rownames(ANOVAshellFC) <- rownames(dat)
  allFCNames <- c()
  
  for(i in 1:nrowShellp){
    
    aovDataFrame <- cbind.data.frame(dat[i,], classes)
    firstColName <- rownames(dat)[i]
    colnames(aovDataFrame)[1] <- firstColName
    
    ANOVA <- aov(formula = aovDataFrame[,1] ~ aovDataFrame[,2], data = aovDataFrame)
    ANOVApval <- summary(ANOVA)[[1]][[5]][1]
    if(posthoc == "tukey"){
      tukey <- TukeyHSD(ANOVA)
      posthoc.pvals <- tukey[[1]][,4]
      allPVals <- unname(c(ANOVApval, posthoc.pvals))
    }else{
      posthocResults <- PostHocTest(x = ANOVA, method = posthoc)
      posthoc.pvals <- posthocResults[[1]][,4]
      allPVals <- unname(c(ANOVApval, posthoc.pvals))
    }
    
    if(i == 1){
      colnames(ANOVAshellp) <- c("pval.Ftest", paste("pval.", names(posthoc.pvals), sep = ""))
    }
    
    ANOVAshellp[i,] <- allPVals
    
    allFC <- c()
    aovMeans <- aggregate(aovDataFrame[,1], list(aovDataFrame[,2]), mean)
    for(j in 1:nrow(aovMeans)){
      for(k in 1:nrow(aovMeans)){
        if(j != k){
          if(i == 1){
            FCname <-  paste("log2.ratio.",
                             as.character(aovMeans[j,1]),
                             "-",
                             as.character(aovMeans[k,1]),
                             sep = "")
            allFCNames <- append(x = allFCNames, values = FCname)
          }
          FC <- aovMeans[j,2] - aovMeans[k,2]
          allFC <- append(x = allFC, values = FC)
        }
      }
    }
    
    if(i == 1){
      colnames(ANOVAshellFC) <- allFCNames
    }
    
    ANOVAshellFC[i,] <- allFC
  }
  ANOVAq <- p.adjust(p = ANOVAshellp[,1], method = adjMet)
  
  result <- cbind(ANOVAshellp[,1], ANOVAq, ANOVAshellp[,-1], ANOVAshellFC)
  colnames(result)[1:2] <- c("pvalue.Ftest", "qvalue.Ftest")
  return(result)
}

## adjANOVA():

adjANOVA <- function(ANOVAres, adjMethod = "fdr"){
  
  postHocInd <- grep(pattern = "pval", x = colnames(ANOVAres))[-1]
  postHoc <- ANOVAres[,postHocInd]
  postHoc <- apply(X = postHoc, MARGIN = 2, FUN = as.numeric)
  rownames(postHoc) <- rownames(ANOVAres)
  postHoc <- apply(X = postHoc, MARGIN = 2, FUN = function(x) p.adjust(p = x, adjMethod))
  colnames(postHoc) <- paste0("adj.", colnames(postHoc), sep = "")
  ANOVAres[,postHocInd] <- postHoc
  colnames(ANOVAres)[postHocInd] <- paste0("adj.", colnames(ANOVAres)[postHocInd], sep = "")
  
  return(ANOVAres)
  
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

## plotProcessBars3():

plotProcessBars3 <- function(processesArr,
                             howMany = 10,
                             main2 = "default",
                             processesExplainedList,
                             obviate.log2FC = FALSE,
                             factorPosition = 4/7){
  
  if(obviate.log2FC){
    red <- "#377EB8"
    lawn <- "#377EB8"
    blue <- "#E41A1C"
  }else{
    red <- "#E41A1C"
    lawn <- "#7CFC00"
    blue <- "#377EB8"
  }
  
  if(main2 == "default"){
    main2 <- "-log10(p-value) of altered GO processes"
  }
  
  processesArr <- processesArr[,"fisher.pvalue"]
  processesArr <- -log10(processesArr)
  processesArr <- processesArr[howMany:1]
  processesNames1 <- names(processesArr)
  processesNames2 <- gsub(pattern = "\\s\\[.*]",
                          replacement = "",
                          x = processesNames1)
  
  space1 <- 1.5
  width <- 1
  nSpaces <- howMany
  firstText <- width + space1
  textPositions <- c(firstText)
  saveValue <- firstText
  for(i in 1:(nSpaces-1)){
    saveValue <- saveValue + firstText
    textPositions <- append(x = textPositions, values = saveValue)
  }
  textExtra <- 0.5
  textPositions <- textPositions + textExtra
  
  specificExplained <- processesExplainedList[processesNames1]
  upOrDownVec <- c()
  for(i in 1:length(specificExplained)){
    logSum <- sum(specificExplained[[i]][,"log2.FC"])
    if(logSum < 0){
      upOrDownVec <- append(x = upOrDownVec, values = lawn)
    }
    if(logSum >= 0){
      upOrDownVec <- append(x = upOrDownVec, values = red)
    }
  }
  
  barplot(height = processesArr,
          horiz = TRUE,
          beside = TRUE,
          col = upOrDownVec,
          xlim = range(pretty(c(0, max(processesArr)))),
          main = main2,
          yaxt = "n",
          space = space1,
          xlab = "-log10(p-value)",
          ylab = "Biological processes")
  
  meanPoint <- range(pretty(c(0, max(processesArr))))[2] * factorPosition
  
  abline(v = -log10(0.05), lty = 2, col = blue, lwd = 2)
  text(x = meanPoint, y = textPositions, labels = processesNames2)
  
}
