# ) retrieveGeneID():

# 

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

## 4) PCACal():

# It computes the PCA object from the given data. To do so, the R base function
# employed is "prcomp()".

PCACal <- function(scaledData, wanna.scale = FALSE){
  
  PCA <- prcomp(x = t(na.omit(scaledData)), scale. = wanna.scale)
  return(PCA)
  
}

PCA <- PCACal(scaledData = scData)

# ) justChangeIt():

# 

justChangeIt <- function(dataset, classes){
  
  for(i in 1:nrow(dataset)){
    
    for(j in 1:length(classes)){
      
      jind <- classes == classes[j]
      
      if(length(unique(dataset[i, jind])) == 1){
        
        dataset[i,jind] <- abs(rnorm(n = length(dataset[i,jind]),
                                     mean = dataset[i,1],
                                     sd = dataset[i,1]/100))
      
      }
      
    }
    
  }
  
  return(dataset)

}

## 5) elbowM():

# Elbow method is a criterion used in machine learning after graphical
# representation. It is a rule followed after the value that produces the most
# drastic reduction in data variance for choosing the number of clusters or
# others for machine learning models. This is also a method utilized for PCA.
# When choosing the two principal components for an scatterplot, elbow method
# can be used to select those two principal components by finding the one that
# implies the most relative variance reduction.

elbowM <- function(vec, namedVec = FALSE){
  
  ratesVec <- c()
  diffsVec <- c()
  len <- length(vec)
  for(i in 2:len){
    if(i <= (len-1)){
      rate1 <- unname(vec[i-1]) - unname(vec[i])
      rate2 <- unname(vec[i]) - unname(vec[i+1])
      rate <- round(rate1/rate2, 5)
      ratesVec <- append(x = ratesVec, values = rate)
      diff <- vec[i-1] - vec[i]
      diffsVec <- append(x = diffsVec, values = diff)
    }else{
      diff <- vec[i-1] - vec[i]
      diffsVec <- append(x = diffsVec, values = diff)
    }
  }
  ratesVec <- c(NA, ratesVec, NA)
  copyDiffs <- diffsVec
  diffsVec <- c(NA, diffsVec)
  dataF <- rbind(vec, ratesVec, diffsVec)
  rownames(dataF) <- c("Original vector",
                       "Continuous difference rate",
                       "Continuous difference")
  if(namedVec == TRUE){
    colnames(dataF) <- names(vec)
  }else{
    colnames(dataF) <- 1:len
  }
  print(dataF)
  orderedVec <- copyDiffs[order(copyDiffs, decreasing = TRUE)]
  if(all(as.numeric(copyDiffs) == as.numeric(orderedVec))){
    return(paste0("Recommended plotting components: ",
                  as.character(colnames(dataF)[which.max(ratesVec)-1]),
                  ", ",
                  as.character(colnames(dataF)[which.max(ratesVec)]),
                  sep = ""))
  }else{
    return(paste0("Recommended plotting components: ",
                  as.character(colnames(dataF)[which.max(diffsVec)-1]),
                  ", ",
                  as.character(colnames(dataF)[which.max(diffsVec)]),
                  sep = ""))
  }
}

## 9) PCAPlot():

# It plots a PCA representation. It requires a PCA object as input.

PCAPlot <- function(PCAObject,
                    is.chosen = FALSE,
                    firstPC,
                    secPC,
                    groups,
                    colVec,
                    legPos = "topright",
                    main2 = "default",
                    cexLeg = 1,
                    cexPoints = 1,
                    cexMain = 1,
                    cexLabels = 1,
                    cexAxes = 1
                    ){
  if(is.chosen == FALSE){
    
    elbow <- elbowM(summary(PCAObject)$importance[2,], namedVec = FALSE)
    
    firstPC <- as.numeric(substr(x = elbow,
                                 start = nchar(elbow)-3,
                                 stop = nchar(elbow)-3))
    secPC <- as.numeric(substr(x = elbow,
                               start = nchar(elbow),
                               stop = nchar(elbow)))
  }
  
  if(main2 == "default"){
    
    main2 <- paste0("PC",
                    as.character(firstPC),
                    " vs. PC",
                    as.character(secPC),
                    sep = "")
    
  }
  
  plot(x = PCAObject$x[,firstPC], y = PCAObject$x[,secPC],
       xlim = range(pretty(c(min(PCAObject$x[,firstPC]), max(PCAObject$x[,firstPC])))),
       ylim = range(pretty(c(min(PCAObject$x[,secPC]), max(PCAObject$x[,secPC])))),
       main = main2,
       xlab = paste0("PC", as.character(firstPC), " (",
                     as.character(round(100*summary(PCAObject)$importance[2,firstPC],1)),
                     "%)",
                     sep = ""),
       ylab = paste0("PC", as.character(secPC), " (",
                     as.character(round(100*summary(PCAObject)$importance[2,secPC],1)),
                     "%)",
                     sep = ""),
       cex = cexPoints,
       cex.main = cexMain,
       cex.lab = cexLabels,
       cex.axis = cexAxes)
  
  for(i in 1:length(groups)){
    
    idx <- which(grepl(pattern = groups[i], x = rownames(PCAObject$x), fixed = TRUE))
    
    points(x = PCAObject$x[idx,firstPC],
           y = PCAObject$x[idx,secPC],
           pch = 21,
           col = "black",
           bg = colVec[i],
           cex = cexPoints)
    
  }
  
  legend(x = legPos,
         legend = groups,
         pch = 21,
         col = "black",
         pt.bg = colVec,
         cex = cexLeg)
  
}

## 10) kmcluster():

# It executes k-means and plots its results in a PCA plot, so sample
# classification can be observed in a tow-dimensional plot. It requires a
# normalized data matrix.

kmcluster <- function(data,
                      maxClus = 10,
                      maxIter = 100,
                      chCol = c("lawngreen", "aquamarine4", "skyblue",
                                "palevioletred2", "lemonchiffon4", "orange", "brown2",
                                "darkorchid", "yellow", "navyblue"),
                      main2 = "default",
                      legPos = "bottomright",
                      is.chosen = TRUE,
                      firstPC = 1,
                      secPC = 2,
                      wanna.scale2 = FALSE,
                      cexLeg = 1,
                      cexPoints = 1,
                      cexMain = 1,
                      cexLabels = 1,
                      cexAxes = 1){
  
  tWith = rep(NA, maxClus)
  for(i in 1:maxClus){
    KM = kmeans(x = t(data), center = i, nstart = maxIter) # transposes data
    tWith[i] = KM$tot.withinss
  }
  barplot(height = tWith, names.arg = 1:length(tWith),
          main = "Total within-cluster sum of squares per number of clusters",
          xlab = "Number of clusters",
          ylab = "Total within-cluster sum of squares",
          ylim = range(pretty(c(0,max(tWith)))))
  eRes <- elbowM(vec = tWith, namedVec = FALSE)
  chosenK <- as.numeric(substr(x = eRes,
                               start = nchar(eRes),
                               stop = nchar(eRes)))
  KM = kmeans(x = t(data), center = chosenK, nstart = maxIter)
  KClasses <- KM$cluster
  KMPalette <- sample(x = chCol, size = chosenK, replace = FALSE)
  
  PCAObj <- PCACal(scaledData = data, wanna.scale = wanna.scale2)
  if(!is.chosen){
    eRes <- elbowM(summary(PCAObj)$importance[2,], namedVec = FALSE)
    PCx <- as.numeric(substr(x = eRes, start = nchar(eRes)-3, stop = nchar(eRes)-3))
    PCy <- as.numeric(substr(x = eRes, start = nchar(eRes), stop = nchar(eRes)))
  }else{
    PCx <- firstPC
    PCy <- secPC
  }
  
  
  if(main2 == "default"){
    main2 = paste0("PC", as.character(PCx),
                   " vs. PC", as.character(PCy),
                   sep = "")
  }
  
  plot(x = PCAObj$x[,PCx], y = PCAObj$x[,PCy],
       xlim = range(pretty(c(min(PCAObj$x[,PCx]), max(PCAObj$x[,PCx])))),
       ylim = range(pretty(c(min(PCAObj$x[,PCy]), max(PCAObj$x[,PCy])))),
       main = main2,
       xlab = paste0("PC", as.character(PCx), " (",
                     as.character(round(100*summary(PCAObj)$importance[2,PCx],1)),
                     "%)",
                     sep = ""),
       ylab = paste0("PC", as.character(PCy), " (",
                     as.character(round(100*summary(PCAObj)$importance[2,PCy],1)),
                     "%)",
                     sep = ""),
       cex = cexPoints,
       cex.main = cexMain,
       cex.lab = cexLabels,
       cex.axis = cexAxes)
  
  clus <- c()
  for(i in 1:chosenK){
    bool <- unname(i == KClasses)
    points(x = PCAObj$x[bool,PCx],
           y = PCAObj$x[bool,PCy],
           pch = 21,
           col = "black",
           bg = KMPalette[i],
           cex = cexPoints)
    clus <- append(x = clus, values = paste0("Cluster ", as.character(i), sep = ""))
  }
  
  legend(x = legPos,
         legend = clus,
         pch = 21,
         col = "black",
         pt.bg = KMPalette,
         cex = cexLeg)
  
}

## 3) heat() # requires "gplots" package

# This function requires the "gplots" package. An input containing only the data
# in the form of a matrix or data frame is needed, as well as the distance
# method calculation as a string. The data only have to contain the samples and
# protein levels in their columns and rows, respectively .Possible distance
# calculations are: "euclidean", "manhattan" and "correlation". The key index
# shows the deviation from a mean value in terms of Z-score (number of standard
# deviations) calculated using the protein levels from all samples. In each
# axis, a dendrogram is plotted. While the distance function is specified,
# clustering method is not.
# heat() function does require NA values omission.

heat <- function(dataF,
                 distance = "correlation",
                 xcex = 0.75,
                 ycex = 0.75,
                 ylab = TRUE,
                 rowD = TRUE,
                 colD = TRUE,
                 colSet = "greenred",
                 rotateXLabs = NULL,
                 rotateYLabs = NULL,
                 ADJCOL = c(NA, 0),
                 ADJROW = c(0, NA)){
  if(!("package:gplots" %in% base::search())){
    library(gplots)
  }
  
  if(colSet == "greenred"){
    colsUsed <- colorRampPalette(colors = c("green", "black", "red"))(75)
  }
  if(colSet == "bluered"){
    colsUsed <- colorRampPalette(colors = c("blue", "black", "red"))(75)
  }
  if(colSet == "redgreen"){
    colsUsed <- colorRampPalette(colors = c("red", "black", "green"))(75)
  }
  
  if(ylab){
    if(distance == "euclidean" | distance == "manhattan"){
      heatmap.2(x = dataF,
                distfun = function(x) {dist(x, method = distance)},
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
    if(distance == "correlation"){
      heatmap.2(x = dataF,
                distfun = function(x)
                  as.dist((1 - cor(t(x), use = "pairwise.complete.obs"))/2),
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
  }else{
    if(distance == "euclidean" | distance == "manhattan"){
      heatmap.2(x = dataF, distfun = function(x) {dist(x, method = distance)},
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                labRow = ylab,
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
    if(distance == "correlation"){
      heatmap.2(x = dataF,
                distfun = function(x)
                  as.dist((1 - cor(t(x), use = "pairwise.complete.obs"))/2),
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                labRow = ylab,
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
  }
}

## ) multiANOVA3():

# 

multiANOVA3 <- function(dat, classes, adjMet = "BH", posthoc = "tukey"){
  
  if(!("package:DescTools" %in% search())){
    library(DescTools)
  }
  
  n <- nlevels(classes)
  nrowShellp <- nrow(dat)
  ncolShellp <- 1 + (n * ((n - 1)/2))
  ANOVAshellp <- matrix(data = rep(NA, (nrowShellp * ncolShellp)),
                        nrow = nrowShellp)
  rownames(ANOVAshellp) <- rownames(dat)
  
  # Preparing p-values shell
  
  nrowShellFC <- nrow(dat)
  ncolShellFC <- n * (n - 1)
  ANOVAshellFC <- matrix(data = rep(NA, (nrowShellFC * ncolShellFC)),
                         nrow = nrowShellFC)
  rownames(ANOVAshellFC) <- rownames(dat)
  allFCNames <- c()
  
  for(i in 1:nrowShellp){
    
    # p-value shell is filled with p-values 
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
    # Name p-values shell columns
    if(i == 1){
      colnames(ANOVAshellp) <- c("pval.Ftest", paste("pval.", names(posthoc.pvals), sep = ""))
    }
    
    ANOVAshellp[i,] <- allPVals
    
    # Preparing fold changes shell
    allFC <- c()
    aovMeans <- aggregate(aovDataFrame[,1], list(aovDataFrame[,2]), mean)
    for(j in 1:nrow(aovMeans)){
      for(k in 1:nrow(aovMeans)){
        if(j != k){
          # Name fold changes shell columns.
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

## ) adjANOVA():

# 

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


## ) rescueAdj():

#

rescueAdj <- function(pvalsArr, chosen.pval = 0.05){
  
  chANOVA <- c()
  
  for(i in 1:nrow(pvalsArr)){
    
    if(sum(pvalsArr[i,] <= chosen.pval) > 0){
      chANOVA <- append(x = chANOVA, values = rownames(pvalsArr)[i])
    }
    
  }
  
  return(chANOVA)

}

## ) specificAdj():

# It only extracts those proteins that are considered significant in the
# comparisons indicated by the user, as a binary or boolean vector, considering
# the structure the ANOVA matrix results has.

specificAdj <- function(pvalsArr, zerosAndOnes){
  
  chANOVA <- c()
  
  for(i in 1:nrow(pvalsArr)){
    
    if(sum(as.numeric(pvalsArr[i,] <= 0.05) == zerosAndOnes) == length(zerosAndOnes)){
      chANOVA <- append(x = chANOVA, values = rownames(pvalsArr)[i])
    }
    
  }
  
  return(chANOVA)

}

## ) multiROC2():

# 

multiROC2 <- function(realClasses,
                      probabilities){
  
  if(!("package:verification" %in% search())){
    library(verification)
  }
  
  allClasses <- unique(realClasses)
  allClasses <- allClasses[order(allClasses)]
  
  if(length(allClasses) == 2){
    
    controlClass <- allClasses[1]
    casesClass <- allClasses[2]
    
    newX <- as.numeric(realClasses == casesClass)
    predProb <- probabilities[,casesClass] # should have lower values associated
    # to control class.
    
    finallyAGoodROC <- verification::roc.plot(x = newX, # a binary observation
                                              # (coded{0, 1})
                                              pred = predProb, # a probability
                                              # prediction
                                              xlab =  "1 - specificity",
                                              ylab = "Sensitivity"
    )
    
    return(finallyAGoodROC)
    
  }else{ # I'm not sure if this will need a value conversion for probabilities. 
    
    listOfROCs <- list()
    
    for(i in 1:length(allClasses)){
      
      iclass <- allClasses[i]
      irealClasses <- as.numeric(realClasses != iclass)
      ipredProb <- probabilities[,iclass]
      iROC <- verification::roc.plot(x = irealClasses, # a binary observation
                                     # (coded{0, 1})
                                     pred = ipredProb, # a probability
                                     # prediction
                                     xlab =  "1 - specificity",
                                     ylab = "Sensitivity"
      )
      listOfROCs[[iclass]] <- iROC
      
    }
    
    return(listOfROCs)
    
  }
  
}

## 29) multiBox():

# 

multiBox <- function(protValues,
                     protName,
                     groups,
                     pvals,
                     onlySig = TRUE,
                     chosenPVals = FALSE,
                     whatPVals,
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

## ) bioProcesses():

# 

bioProcesses<- function(uniProtAccessionVector, pvalVec){
  
  if(!("package:UniprotR" %in% search())){
    library(UniprotR)
  }
  
  GOlist <- list()
  
  proteinProcesses <- GetProteinGOInfo(ProteinAccList = uniProtAccessionVector, directorypath = NULL)[,3]
  
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
  
  return(processesMatrix)
  
}

## 21) bioProcesses2():

# 

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
  zero[trueProteins,] <- refinedMatrix[trueProteins,]
  
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

## ) allProcessesExplained():

# 

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
  zero[trueSelectedProteins,] <- refinedMatrix[trueSelectedProteins,]
  
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

## ) printTopProcesses():

# 

printTopProcesses <- function(processesMatrix,
                              processesExplained,
                              filePath,
                              howMany = 10){
  
  if(!("package:xlsx" %in% search())){
    library(xlsx)
  }
  
  topProcesses <- rownames(processesMatrix)[1:howMany]
  subList <- processesExplained[topProcesses]
  topProcesses <- sprintf("process%s", 1:howMany)
  
  for(i in 1:length(subList)){
    if(i == 1){
      write.xlsx2(x = subList[[i]],
                  file = filePath,
                  sheetName = topProcesses[i],
                  col.names = TRUE,
                  row.names = FALSE,
                  append = FALSE)
    }else{
      write.xlsx2(x = subList[[i]],
                  file = filePath,
                  sheetName = topProcesses[i],
                  col.names = TRUE,
                  row.names = FALSE,
                  append = TRUE)
    }
  }
  
}

## plotProcessBars():

# 

plotProcessBars <- function(processesArr,
                            howMany = 10,
                            main2 = "default",
                            processesExplainedList){
  
  red <- "#E41A1C"
  blue <- "#377EB8"
  green <- "#4DAF4A"
  
  if(main2 == "default"){
    main2 <- "-log10(p- and q-values) of altered GO processes"
  }
  
  processesArr <- processesArr[,c("fisher.qvalue","fisher.pvalue")]
  processesArr <- -log10(processesArr)
  processesArr <- processesArr[howMany:1,]
  processesArr <- t(processesArr)
  processesNames <- colnames(processesArr)
  processesNames <- gsub(pattern = "\\s\\[.*]",
                         replacement = "",
                         x = processesNames)
  
  space1 <- howMany / 50
  space2 <- howMany * 0.2
  nSpaces <- howMany - 1
  firstSpace <- space2 + 1 + space1 + 1
  spacesPosition <- c(firstSpace)
  saveValue <- firstSpace
  for(i in 1:nSpaces){
    saveValue <- saveValue + firstSpace
    spacesPosition <- append(x = spacesPosition, values = saveValue)
  }
  spacesPosition <- spacesPosition + (space1 * 4)
  
  barplot(height = processesArr,
          horiz = TRUE,
          beside = TRUE,
          col = c(green, blue),
          xlim = range(pretty(processesArr)),
          main = main2,
          yaxt = "n",
          space = c(space1, space2),
          xlab = "-log10(p-/q-value)",
          ylab = "Biological processes")
  
  meanPoint <- range(pretty(processesArr))[2] * 4/7
  
  abline(v = -log10(0.05), lty = 2, col = red, lwd = 2)
  text(x = meanPoint, y = spacesPosition, labels = processesNames)
  
}

## ) plotProcessesBar2():

# 

plotProcessBars2 <- function(processesArr,
                             howMany = 10,
                             main2 = "default",
                             processesExplainedList){
  
  orange <- "#FF8C00"
  red <- "#E41A1C"
  lawn <- "#7CFC00"
  sea <- "#20B2AA"
  blue <- "#377EB8"
  
  downColors <- c(red, orange)
  upColors <- c(sea, lawn)
  
  if(main2 == "default"){
    main2 <- "-log10(p- and q-values) of altered GO processes"
  }
  
  processesArr <- processesArr[,c("fisher.qvalue","fisher.pvalue")]
  processesArr <- -log10(processesArr)
  processesArr <- processesArr[howMany:1,]
  processesArr <- t(processesArr)
  processesNames1 <- colnames(processesArr)
  processesNames2 <- gsub(pattern = "\\s\\[.*]",
                          replacement = "",
                          x = processesNames1)
  
  space1 <- howMany / 50
  space2 <- howMany * 0.2
  nSpaces <- howMany - 1
  firstSpace <- space2 + 1 + space1 + 1
  spacesPosition <- c(firstSpace)
  saveValue <- firstSpace
  for(i in 1:nSpaces){
    saveValue <- saveValue + firstSpace
    spacesPosition <- append(x = spacesPosition, values = saveValue)
  }
  spacesPosition <- spacesPosition + (space1 * 4)
  
  specificExplained <- processesExplainedList[processesNames1]
  upOrDownVec <- c()
  for(i in 1:length(specificExplained)){
    
    upOrDownCol <- specificExplained[[i]][,ncol(specificExplained[[i]])]
    
    if(length(upOrDownCol) == 1){
      mostCommon <- upOrDownCol
    }else{
      aggregateDF <- data.frame("log2.FC" = specificExplained[[i]][,"log2.FC"],
                                "is.Up.Or.Down" = upOrDownCol)
      rownames(aggregateDF) <- specificExplained[[i]][,"protein.code"]
      aggregateResults <- aggregate(aggregateDF$log2.FC,
                                    by = list(Category = aggregateDF$is.Up.Or.Down),
                                    FUN = sum)
      print(aggregateResults)
      mostCommon <- aggregateResults[which.max(abs(aggregateResults$x)),1]
    }
    
    if(mostCommon == "Down"){
      upOrDownVec <- append(x = upOrDownVec, values = downColors)
    }
    if(mostCommon == "Up"){
      upOrDownVec <- append(x = upOrDownVec, values = upColors)
    }
    
  }
  
  barplot(height = processesArr,
          horiz = TRUE,
          beside = TRUE,
          col = upOrDownVec,
          xlim = range(pretty(processesArr)),
          main = main2,
          yaxt = "n",
          space = c(space1, space2),
          xlab = "-log10(p-/q-value)",
          ylab = "Biological processes")
  
  meanPoint <- range(pretty(processesArr))[2] * 4/7
  
  abline(v = -log10(0.05), lty = 2, col = blue, lwd = 2)
  text(x = meanPoint, y = spacesPosition, labels = processesNames2)
  
}

## 36) howManyLoops():

# 

howManyLoops <- function(nfeat, wantedLoops){
  prob <- seq(1 - 0.001, 0.001, by = -0.001)
  iterationsVector <- c()
  for(i in 1:length(prob)){
    remaining <- c()
    remember <- nfeat
    k <- 1
    while(remember > 1){
      remember <- floor(remember * prob[i])
      remaining <- append(x = remaining, values = remember)
    }
    iterationsVector <- append(x = iterationsVector, values = length(remaining))
  }
  
  return(prob[iterationsVector == wantedLoops])
}

## 37) recursiveRandomForest3():

# 

recursiveRandomForest3 <- function(dataset,
                                   groups,
                                   treeNumber = 1000,
                                   retentionProportion = 0.8,
                                   chosenIterationsNumber = 5,
                                   withNormalization = TRUE,
                                   topNumber = 25,
                                   iterationsBest = 100,
                                   searchForBest = TRUE,
                                   howManyFeats = 8,
                                   minimization = TRUE){
  
  if(!("package:randomForest" %in% base::search())){
    library(randomForest)
  }
  
  importanceMatrix <- matrix(data = 0,
                             nrow = ncol(dataset),
                             ncol = 4)
  rownames(importanceMatrix) <- colnames(dataset)
  colnames(importanceMatrix) <- c("feature.index",
                                  "importance",
                                  "normalized.importance",
                                  "standard.cumulative")
  importanceMatrix[,"feature.index"] <- 1:nrow(importanceMatrix)
  
  iterations <- 1
  
  while(iterations <= chosenIterationsNumber){
    
    print(paste("Iteration: ", as.character(iterations), sep = ""))
    tempDataset <- dataset
    
    retentionLoop <- 1
    
    while(!(is.null(ncol(tempDataset)))){
      
      tempRF <- randomForest(x = tempDataset,
                             y = groups,
                             importance = TRUE,
                             ntree = treeNumber)
      
      tempImportance <- tempRF$importance
      
      vecImp <- tempImportance[,"MeanDecreaseAccuracy"]
      
      boolImp <- vecImp > 0
      boolNoImp <- vecImp < 0
      
      if(withNormalization){
        sumImp <- sum(vecImp[boolImp])
        sumNoImp <- sum(vecImp[boolNoImp])
        vecImp[boolImp] <- vecImp[boolImp] / sumImp
        vecImp[boolNoImp] <- vecImp[boolNoImp] / sumNoImp
      }
      
      importanceMatrix[rownames(tempImportance),"importance"] <-
        importanceMatrix[rownames(tempImportance),"importance"] + vecImp
      
      newOrder <- order(vecImp, decreasing = TRUE)
      tempImportance <- tempImportance[newOrder,]
      # View(tempImportance)
      
      retention <- floor(nrow(tempImportance) * retentionProportion)
      print(paste("Retained features (iteration: ",
                  as.character(iterations),
                  "; retention loop: ",
                  as.character(retentionLoop),
                  "): ",
                  as.character(retention),
                  sep = ""))
      featsRetained <- rownames(tempImportance)[1:retention]
      
      tempDataset <- tempDataset[,featsRetained]
      
      retentionLoop <- retentionLoop + 1
    }
    
    iterations <- iterations + 1
  }
  
  orderedMatrix <- order(importanceMatrix[,"importance"], decreasing = TRUE)
  orderedMatrix <- importanceMatrix[orderedMatrix,]
  
  boolOrdGreat <- orderedMatrix[,"importance"] > 0
  boolOrdLess <- orderedMatrix[,"importance"] < 0
  
  sumGreat <- sum(orderedMatrix[boolOrdGreat,"importance"])
  sumLess <- sum(orderedMatrix[boolOrdLess,"importance"])
  
  orderedMatrix[boolOrdGreat,"normalized.importance"] <-
    orderedMatrix[boolOrdGreat,"importance"] / sumGreat
  orderedMatrix[boolOrdLess,"normalized.importance"] <-
    -(orderedMatrix[boolOrdLess,"importance"] / sumLess)
  
  orderedMatrix[,"standard.cumulative"] <- cumsum(orderedMatrix[,"normalized.importance"])
  
  print("Importance matrix calculation ended. Generating best random forest.")
  
  iterations <- 1
  
  oobFeats <- c()
  oobBest <- 1
  
  top <- rownames(orderedMatrix)[1:topNumber]
  
  bestDataset <- dataset[,top]
  
  while(iterations <= iterationsBest){
    
    tempDataset <- bestDataset
    if(searchForBest){
      print(paste("Number of iteration: ", as.character(iterations), sep = ""))
    }
    
    if(searchForBest){
      
      while(!(is.null(ncol(tempDataset)))){
        
        tempRF <- randomForest(x = tempDataset,
                               y = groups,
                               importance = TRUE,
                               proximity = TRUE,
                               ntree = treeNumber)
        
        conf <- tempRF$confusion[,-ncol(tempRF$confusion)]
        oob <- 1 - (sum(diag(conf)) / sum(sum(conf)))
        
        if(minimization){
          if(oob < oobBest){
            print(paste("New classifier found. OOB: ",
                        as.character(round(oob, 4)),
                        sep = ""))
            oobBest <- oob
            oobFeats <- colnames(tempDataset)
            bestRF <- tempRF
          }
          if(oob == oobBest & length(oobFeats) > length(colnames(tempDataset))){
            print(paste("Number of proteins reduced in new classifier. New number of proteins: ",
                        as.character(length(colnames(tempDataset))),
                        sep = ""))
            oobBest <- oob
            oobFeats <- colnames(tempDataset)
            bestRF <- tempRF
          }
        }else{
          if(oob < oobBest){
            print(paste("New classifier found. OOB: ",
                        as.character(round(oob, 4)),
                        sep = ""))
            oobBest <- oob
            oobFeats <- colnames(tempDataset)
            bestRF <- tempRF
          }
        }
        
        importanceRF <- tempRF$importance[,"MeanDecreaseAccuracy"]
        retainedFeats <- rownames(tempRF$importance[order(importanceRF, decreasing = TRUE),])
        retainedFeats <- retainedFeats[-length(retainedFeats)]
        tempDataset <- tempDataset[,retainedFeats]
        
      }
      
    }else{
      
      while(!(is.null(ncol(tempDataset)))){
        
        tempRF <- randomForest(x = tempDataset,
                               y = groups,
                               importance = TRUE,
                               proximity = TRUE,
                               ntree = treeNumber)
        
        if(ncol(tempDataset) == howManyFeats){
          print(paste("Number of iteration (selected number of features: ",
                      as.character(howManyFeats),
                      "): ",
                      as.character(iterations),
                      sep = ""))
          conf <- tempRF$confusion[,-ncol(tempRF$confusion)]
          oob <- 1 - (sum(diag(conf)) / sum(sum(conf)))
          if(oob < oobBest){
            print(paste("New classifier found. OOB: ",
                        as.character(round(oob, 4)),
                        sep = ""))
            oobBest <- oob
            oobFeats <- colnames(tempDataset)
            bestRF <- tempRF
          }
        }
        
        importanceRF <- tempRF$importance[,"MeanDecreaseAccuracy"]
        retainedFeats <- rownames(tempRF$importance[order(importanceRF, decreasing = TRUE),])
        retainedFeats <- retainedFeats[-length(retainedFeats)]
        tempDataset <- tempDataset[,retainedFeats]
        
      }
      
    }
    
    iterations <- iterations + 1
    if(oobBest == 0){
      
    }
    
  }
  
  print("Done.")
  
  resultsList <- list()
  resultsList[[1]] <- orderedMatrix
  resultsList[[2]] <- oobFeats
  resultsList[[3]] <- bestRF
  
  names(resultsList) <- c("importance.matrix", "chosen.features", "best.RF")
  
  return(resultsList)
  
}

## 11) MDS2():

# It plots a multidimensional scaling (MDS) plot from the results of a random
# forest classification if "proximity = TRUE" was specified when creating the
# random forest.

MDS2 <- function(rfObj,
                 classCol,
                 areChosenColors = FALSE,
                 chosenColors,
                 main2 = "Multidimensional scaling (MDS) plot",
                 legPos = "bottomleft",
                 pointCex = 1.2,
                 legSize = 1,
                 cexMain = 1,
                 cexLabels = 1,
                 cexAxes = 1){
  
  if(!("package:randomForest" %in% base::search())){
    library(randomForest)
  }
  if(!("package:RColorBrewer" %in% base::search())){
    library(RColorBrewer)
  }
  
  prox <- cmdscale(1 - rfObj$proximity, eig = TRUE)
  uc <- unique(as.character(classCol))
  if(areChosenColors){
    myPal <- chosenColors
  }else{
    myPal <- brewer.pal(n = length(uc), name = "Set1")[1:length(uc)]
  }
  
  plot(x = prox$points[,1],
       y = prox$points[,2],
       xlab = "Dim. 1",
       ylab = "Dim. 2",
       xlim = range(pretty(c(min(prox$points[,1]), max(prox$points[,1])))),
       ylim = range(pretty(c(min(prox$points[,2]), max(prox$points[,2])))),
       pch = 21,
       col = "black",
       main = main2,
       cex = pointCex,
       cex.main = cexMain,
       cex.lab = cexLabels,
       cex.axis = cexAxes)
  
  for(i in 1:length(uc)){
    bool <- classCol == uc[i]
    points(x = prox$points[,1][bool],
           y = prox$points[,2][bool],
           pch = 21,
           col = "black",
           bg = myPal[i],
           cex = pointCex)
  }
  
  legend(x = legPos,
         legend = uc,
         pch = 21,
         col = "black",
         pt.bg = myPal,
         cex = legSize)
  
}

## 18) MDSplotter3D():

# 

MDSplotter3D <- function(randomF,
                         groupsVector,
                         main3D = "3D multidimensional scaling (MDS) plot",
                         type3D = "s",
                         size3D = 0.6,
                         is.biplot = FALSE,
                         legPos = "topright",
                         legSize = 1.5,
                         legInset = c(0.05),
                         is.colors = FALSE,
                         chosenColors){
  if(!("package:randomForest" %in% base::search())){
    library(randomForest)
  }
  if(!("package:RColorBrewer" %in% base::search())){
    library(RColorBrewer)
  }
  if(!("package:stringi" %in% base::search())){
    library(stringi)
  }
  if(!("package:rgl" %in% base::search())){
    library(rgl)
  }
  
  uniqueGroups <- unique(as.character(groupsVector))
  if(is.colors){
    tempPalette <- chosenColors
  }else{
    tempPalette <- brewer.pal(n = length(uniqueGroups),
                              name = "Set1")[c(1:length(uniqueGroups))]
  }
  
  sampleColor <- stri_replace_all_regex(as.character(groupsVector),
                                        pattern = uniqueGroups,
                                        replacement = tempPalette,
                                        vectorize = FALSE)
  if(is.biplot){
    dimensions <- MDSplot(rf = randomF, palette = tempPalette, k = 3, fac = groupsVector)$points
  }else{
    dimensions <- cmdscale(1 - randomF$proximity, eig = TRUE, k = 3)$points
  }
  
  xDif <- pretty(c(min(dimensions[,1]), max(dimensions[,1])))[2] -
    pretty(c(min(dimensions[,1]), max(dimensions[,1])))[1]
  xRangeMin <- min(dimensions[,1]) - xDif
  xRangeMax <- max(dimensions[,1]) + xDif
  
  yDif <- pretty(c(min(dimensions[,2]), max(dimensions[,2])))[2] -
    pretty(c(min(dimensions[,2]), max(dimensions[,2])))[1]
  yRangeMin <- min(dimensions[,2]) - yDif
  yRangeMax <- max(dimensions[,2]) + yDif
  
  zDif <- pretty(c(min(dimensions[,3]), max(dimensions[,3])))[2] -
    pretty(c(min(dimensions[,3]), max(dimensions[,3])))[1]
  zRangeMin <- min(dimensions[,3]) - zDif
  zRangeMax <- max(dimensions[,3]) + zDif
  
  plot3d(x = dimensions[,1],
         y = dimensions[,2],
         z = dimensions[,3],
         main = main3D,
         xlab = "Dimension 1",
         ylab = "Dimension 2",
         zlab = "Dimension 3",
         type = type3D,
         size = size3D,
         col = sampleColor,
         xlim = range(pretty(c(xRangeMin, xRangeMax))),
         ylim = range(pretty(c(yRangeMin, yRangeMax))),
         zlim = range(pretty(c(zRangeMin, zRangeMax))))
  
  legend3d(x = legPos, legend = uniqueGroups, col = tempPalette, pch = 19,
           cex = legSize, inset = legInset)
}