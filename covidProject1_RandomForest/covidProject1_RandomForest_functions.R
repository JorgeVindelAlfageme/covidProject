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

## howManyLoops():

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

## recursiveRandomForest3():

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

## MDS2():

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
       xlab = "Multidimensional scale (MDS) dimension 1",
       ylab = "Multidimensional scale (MDS) dimension 2",
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

## MDSplotter3D():

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

## heat():

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

## multiBox2():

multiBox2 <- function(protValues,
                      protName,
                      groups,
                      pvals,
                      onlySig = TRUE,
                      chosenPVals = FALSE,
                      whatPVals,
                      textCex = 1.5,
                      main2,
                      isNorm = TRUE){
  
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
  
  if(isNorm){
    boxplotYLab <- paste(protName, " normalized levels", sep = "")
  }else{
    boxplotYLab <- paste(protName, " raw levels", sep = "")
  }
  
  boxplot(protDataF[,1] ~ protDataF[,2],
          data = protDataF,
          pch = 21,
          bg = "grey",
          xaxt = "n",
          # xlab = "", # 18/11/2022
          xlab = "Groups",
          ylab = boxplotYLab,
          main = main2,
          ylim = range(pretty(c(min(protValues), finalPretty))))
  
  # axis(side = 1, at = 1:length(unique(groups)), labels = unique(groups), las = 2) # 18/11/2022 las
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

## photoBox2():

photoBox2 <- function(is.already = TRUE,
                      nameDir,
                      where,
                      dat,
                      chosenQuality = 100,
                      chosenWidth = 600,
                      chosenHeight = 480,
                      pvalues,
                      groups2,
                      is.ANOVA = TRUE,
                      main3,
                      isNorm2 = TRUE){
  
  if(!(is.already)){
    where <- paste(where, nameDir, "/", sep = "")
    dir.create(path = where)
  }
  
  for(i in 1:nrow(dat)){
    fileName <- paste(rownames(dat)[i], ".jpeg", sep = "")
    # fileName <- paste(as.character(i), ".jpeg", sep = "") # 18/11/2022
    filePath <- paste(where, fileName, sep = "")
    jpeg(file = filePath,
         quality = chosenQuality,
         width = chosenWidth,
         height = chosenHeight)
    if(is.ANOVA){
      multiBox2(protValues = dat[rownames(dat)[i],],
                protName = rownames(dat)[i],
                groups = groups2,
                pvals = pvalues[i,],
                onlySig = TRUE,
                chosenPVals = FALSE,
                main2 = main3[i],
                isNorm = isNorm2)
    }else{
      multiBox2(protValues = dat[rownames(dat)[i],],
                protName = rownames(dat)[i],
                groups = groups2,
                pvals = pvalues[i],
                onlySig = TRUE,
                chosenPVals = FALSE,
                main2 = main3[i],
                isNorm = isNorm2)
    }
    
    dev.off()
  }
  
}
