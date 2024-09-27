setwd(dir = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/")

dir.create(path = "rda")

table1 <- read.csv2(file = "20230803_Norm_Age_Area_Proteins_Tometaboloanalyst.csv")
table2 <- read.csv2(file = "20230803_Norm_Age_Area_Proteins_Tometaboloanalyst_GroupsbyAge.csv")
table3 <- read.csv2(file = "20230803_Norm_Area_Proteins_Tometaboloanalyst.csv")

library(xlsx)

uniprotDescription <- read.xlsx2(file = "idmapping_2023_08_03.xlsx",
                                 sheetIndex = 1)
rownames(uniprotDescription) <- uniprotDescription[,"Entry"]
uniprotDescription <- uniprotDescription[,c(4,5)]
colnames(uniprotDescription) <- c("function", "gene")

save(x = uniprotDescription, file = "rda/uniprotDescription.rda")

patientsKey <- cbind.data.frame(table1[, c("Sample", "Group", "Age")])
rownames(patientsKey) <- patientsKey[, "Sample"]
rownames(table2) <- table2[,1]
patientsKey[, 4] <- table2[rownames(patientsKey), 2]
colnames(patientsKey)[ncol(patientsKey)] <- "Subgroup"

covidData <- table1[, 4:ncol(table1)]
rownames(covidData) <- table1[,1]
colnames(covidData) <- unlist(lapply(X = strsplit(x = colnames(covidData),
                                                  split = ".",
                                                  fixed = TRUE),
                                     FUN = function(x) x[length(x)]))
covidData <- as.matrix(covidData)

log2CovidData <- log2(covidData)
scaledCovidData <- scale(x = log2CovidData)

covidPCA <- prcomp(x = scaledCovidData, scale. = FALSE)

covidColors <- c("darkgrey", "gold",
                 RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(3, 1)])
names(covidColors) <- unique(patientsKey[,"Group"])
allSamplesColors <- covidColors[patientsKey[,"Group"]]

firstPC <- 1
secondPC <- 2

firstPCName <- paste("PC", as.character(firstPC), sep = "")
secondPCName <- paste("PC", as.character(secondPC), sep = "")

firstLab <- paste(firstPCName, " (",
                  as.character(round(100 * summary(covidPCA)$importance[2, firstPC],1)),
                  " %)",
                  sep = "")
secondLab <- paste(secondPCName, " (",
                   as.character(round(100 * summary(covidPCA)$importance[2, secondPC],1)),
                   " %)",
                   sep = "")

plot(x = covidPCA$x[,firstPCName],
     y = covidPCA$x[,secondPCName],
     xlim = range(pretty(covidPCA$x[,firstPCName])),
     ylim = range(pretty(covidPCA$x[,secondPCName])),
     xaxt = "n",
     yaxt = "n",
     xlab = firstLab,
     ylab = secondLab,
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "PCA (PRM covid experiment)",
     cex.lab = 1.55)
axis(side = 1,
     at = pretty(covidPCA$x[,firstPCName]),
     labels = pretty(covidPCA$x[,firstPCName]),
     las = 1,
     cex.axis = 1.5)
axis(side = 2,
     at = pretty(covidPCA$x[,secondPCName]),
     labels = pretty(covidPCA$x[,secondPCName]),
     las = 1,
     cex.axis = 1.5)

legend(x = "topleft",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.3,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

library(umap)

umapdefault <- umap.defaults

umapdefault$random_state <- 1
umapdefault$n_components <- 2

covidUMAP <- umap::umap(d = scaledCovidData[, chosenFortunaFeats],
                        config = umapdefault,
                        method = "naive",
                        preserve.seed = TRUE)

plot(x = covidUMAP$layout[,1],
     y = covidUMAP$layout[,2],
     xlim = range(pretty(covidUMAP$layout[,1])),
     ylim = range(pretty(covidUMAP$layout[,2])),
     xaxt = "n",
     yaxt = "n",
     xlab = "UMAP dimension 1",
     ylab = "UMAP dimension 2",
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "UMAP (PRM covid experiment)",
     cex.lab = 1.55)
axis(side = 1,
     at = pretty(covidUMAP$layout[,1]),
     labels = pretty(covidUMAP$layout[,1]),
     las = 1,
     cex.axis = 1.5)
axis(side = 2,
     at = pretty(covidUMAP$layout[,2]),
     labels = pretty(covidUMAP$layout[,2]),
     las = 1,
     cex.axis = 1.5)

legend(x = "topright",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.3,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

library(dbscan)

chosenK <- 2 * ncol(covidUMAP$layout)
kdist <- dbscan::kNNdist(x = covidUMAP$layout, k = chosenK)
kdist <- kdist[order(kdist, decreasing = TRUE)]
names(kdist) <- 1:length(kdist)

newElbowMethod2 <- function(numericalVector,
                            wannaPlot = TRUE){
  
  numericalVector <- numericalVector[order(numericalVector, decreasing = TRUE)]
  
  dy <- numericalVector[length(numericalVector)] - numericalVector[1]
  dx <- length(numericalVector) - 1
  originalSlope <- dy/dx 
  
  newElbowRes <- c()
  
  for(i in 1:(length(numericalVector))){
    intercept <- numericalVector[i] - (originalSlope * i)
    newElbowRes <- append(x = newElbowRes, values = intercept)
  }
  
  finalIndex <- which.min(newElbowRes)
  finalValue <- numericalVector[finalIndex]
  
  if(wannaPlot){
    
    repBlack <- rep("black", length(numericalVector))
    
    plot(x = 1:length(numericalVector),
         y = numericalVector,
         xlim = range(pretty(1:length(numericalVector))),
         ylim = range(pretty(numericalVector)),
         xlab = "x",
         ylab = "y",
         xaxt = "n",
         yaxt = "n",
         cex.lab = 2,
         main = "Elbow method analysis (elbow point indicated)",
         cex.main = 2,
         pch = 21,
         col = "black",
         bg = repBlack,
         cex = 0.7)
    axis(side = 1,
         las = 1,
         cex.axis = 1.75,
         at = pretty(1:length(numericalVector)),
         labels = pretty(1:length(numericalVector)))
    axis(side = 2,
         las = 1,
         cex.axis = 1.75,
         at = pretty(numericalVector),
         labels = pretty(numericalVector))
    
    points(x = c(length(numericalVector), 1),
           y = c(numericalVector[length(numericalVector)], numericalVector[1]),
           type = "l",
           col = "forestgreen",
           lty = 1,
           lwd = 2)
    
    y2 <- newElbowRes[finalIndex] + originalSlope * length(numericalVector)
    newdy <- numericalVector[length(numericalVector)] - y2
    y1 <- numericalVector[1] - newdy
    
    points(x = rev(c(length(numericalVector), 1)),
           y = c(y1, y2),
           type = "l",
           col = "forestgreen",
           lty = 1,
           lwd = 2)
    points(x = finalIndex,
           y = numericalVector[finalIndex],
           cex = 0.7,
           col = "red",
           bg = "red",
           pch = 21)
    
    legend(x = "topright",
           legend = c(paste("Calculated elbow point index: ",
                            as.character(round(x = finalIndex, digits = 4)),
                            sep = ""),
                      paste("Calculated elbow point value: ",
                            as.character(round(x = finalValue, digits = 4)),
                            sep = "")),
           inset = c(0.02, 0.02),
           cex = 1.5)
    
  }
  
  return(list("interceptors" = newElbowRes,
              "finalIndex" = finalIndex,
              "finalValue" = finalValue))
  
}

elbowResults <- newElbowMethod2(numericalVector = kdist, wannaPlot = TRUE)
chosenEps <- elbowResults[["finalValue"]]

dbres <- fpc::dbscan(data = covidUMAP$layout,
                     eps = chosenEps,
                     MinPts = chosenK)

dbresCluster <- dbres$cluster

randomNColorSet <- function(ncolors){
  
  if(!("package:RColorBrewer" %in% search())){
    library(RColorBrewer)
  }
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- sample(col_vector, ncolors)
  return(col_vector)
  
}

if(min(dbresCluster) == 0){
  dbresCluster <- dbresCluster + 1
}

set.seed(9)
selectedColors <- randomNColorSet(ncolors = max(dbresCluster))
clustIntoColors <- sapply(X = dbresCluster, FUN = function(x) selectedColors[x])

plot(x = covidUMAP$layout[,1],
     y = covidUMAP$layout[,2],
     xlim = range(pretty(covidUMAP$layout[,1])),
     ylim = range(pretty(covidUMAP$layout[,2])),
     xaxt = "n",
     yaxt = "n",
     xlab = "UMAP dimension 1",
     ylab = "UMAP dimension 2",
     cex.lab = 1.5,
     main = "DBSCAN analysis (predicted groups)",
     cex.main = 1.75,
     pch = 21,
     col = "black",
     bg = clustIntoColors,
     cex = 1.75)
axis(side = 1,
     at = pretty(covidUMAP$layout[,1]),
     labels = pretty(covidUMAP$layout[,1]),
     las = 1,
     cex.axis = 1.5)
axis(side = 2,
     at = pretty(covidUMAP$layout[,2]),
     labels = pretty(covidUMAP$layout[,2]),
     las = 1,
     cex.axis = 1.5)

legend(x = "bottomright",
       legend = sprintf("DBSCAN cluster %s", 1:max(dbresCluster)),
       pch = 21,
       col = "black",
       pt.bg = selectedColors,
       inset = c(0.02, 0.05),
       cex = 1.3,
       pt.cex = 2)

library(limma)

covidData_noBatch <- removeBatchEffect(x = t(covidData_noBatch),
                                       batch = dbresCluster)
covidData_noBatch <- t(covidData_noBatch)

library(caret)

covidClasses <- as.factor(patientsKey[,"Group"])
partitionProp <- 0.7
set.seed(1)
trIndices <- as.numeric(createDataPartition(y = covidClasses,
                                            p = partitionProp,
                                            list = FALSE))
trData <- scaledCovidData[trIndices,]
testData <- scaledCovidData[-trIndices,]
trClasses <- covidClasses[trIndices]
testClasses <- covidClasses[-trIndices]

classColName <- "patient"
trDataFrame <- cbind.data.frame(trData,
                                "class" = trClasses)
colnames(trDataFrame)[ncol(trDataFrame)] <- classColName

library(Boruta)

borutaFormula <- as.formula(paste(classColName, " ~ .", sep = ""))
set.seed(1)
impFeats <- Boruta(borutaFormula, trDataFrame)

impFeats$finalDecision
impFeats$ImpHistory
impFeats$pValue
impFeats$maxRuns
impFeats$light
impFeats$mcAdj
impFeats$timeTaken
impFeats$roughfixed
impFeats$call
impFeats$impSource

borutaBool <- impFeats$finalDecision == "Confirmed" | impFeats$finalDecision == "Tentative"
borutaFeats <- names(impFeats$finalDecision)[borutaBool]

set.seed(1)
currentWrapper <- wrapperCholestasis2(startData = trDataFrame[,c(borutaFeats,
                                                                 classColName)],
                                      classColumn = classColName,
                                      defaultGrid = NULL,
                                      crossValidation = "LOOCV",
                                      bootSamples = 100,
                                      chosenAlgorithm = "svmRadial",
                                      RFchosenNTree = 1000)

predict(object = currentWrapper[["finalClassifier"]],
        newdata = testData[,borutaFeats])

confusionMatrix(data = predict(object = currentWrapper[["finalClassifier"]],
                               newdata = testData[,borutaFeats]),
                reference = testClasses)

nChosenFeats <- 10

set.seed(1)
covidFortuna <- fortuna(dataset = trData,
                        groups = trClasses, # as.factor(groups)
                        treeNumber = 1000,
                        retentionByProp = FALSE,
                        retentionProportion = 0.8,
                        retentionLoops = 30,
                        chosenIterationsNumber = 30,
                        # NAsubstitution = "none", # "none", "sampling", "zero"
                        isBinaryApproach = FALSE,
                        combBinaryClassMethod = "k", # "k", "elbow", "p"
                        chosenPValue = 0.05,
                        pAdjMethod = "fdr",
                        kFeats = nChosenFeats,
                        isElbowPlot = TRUE)

chosenFortunaFeats <- names(which(covidFortuna$importance.matrix$impArr[, "adj.pvalue"] <= 0.05))

newCovidPCA <- prcomp(x = scaledCovidData[, chosenFortunaFeats],
                      scale. = FALSE)

newFirstLab <- paste(firstPCName, " (",
                     as.character(round(100 * summary(newCovidPCA)$importance[2, firstPC],1)),
                     " %)",
                     sep = "")
newSecondLab <- paste(secondPCName, " (",
                      as.character(round(100 * summary(newCovidPCA)$importance[2, secondPC],1)),
                      " %)",
                      sep = "")

plot(x = newCovidPCA$x[,firstPCName],
     y = newCovidPCA$x[,secondPCName],
     xlim = range(pretty(newCovidPCA$x[,firstPCName])),
     ylim = range(pretty(newCovidPCA$x[,secondPCName])),
     xaxt = "n",
     yaxt = "n",
     xlab = newFirstLab,
     ylab = newSecondLab,
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "PCA (PRM covid experiment) (selected important features)",
     cex.lab = 1.5)
axis(side = 1,
     at = pretty(newCovidPCA$x[,firstPCName]),
     labels = pretty(newCovidPCA$x[,firstPCName]),
     las = 1,
     cex.axis = 1.75)
axis(side = 2,
     at = pretty(newCovidPCA$x[,secondPCName]),
     labels = pretty(newCovidPCA$x[,secondPCName]),
     las = 1,
     cex.axis = 1.75)

legend(x = "bottomright",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.4,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

allAlgorithms <- c("rf",
                   "lda",
                   "svmRadial",
                   "glmnet",
                   "knn",
                   "nb",
                   "nnet",
                   "rpart")
completeAlgorithmsName <- c("Random forest (RF)",
                            "Linear Discriminant Analysis (LDA)",
                            "Support Vector Machines (SVM) (Gaussian kernel)",
                            "Logistic regression",
                            "K-Nearest Neighbors (KNN)",
                            "Naive Bayes",
                            "Artificial Neural Network (ANN)",
                            "Decision tree")

performanceMeasurements <- c("accuracy", "kappa", "caret.pvalue")

CVsetResultsShell_rownames <- sprintf(rep(paste(allAlgorithms,
                                                ".combn.%s",
                                                sep = ""),
                                          ((2^length(chosenFortunaFeats)) - 1)),
                                      rep(1:((2^length(chosenFortunaFeats)) - 1),
                                          each = length(allAlgorithms)))

shellForTest <- matrix(data = NA,
                       nrow = length(CVsetResultsShell_rownames),
                       ncol = length(performanceMeasurements),
                       dimnames = list(CVsetResultsShell_rownames,
                                       performanceMeasurements))
CVsetResultsShell <- shellForTest[,c("accuracy", "kappa")]

combnKeys <- matrix(data = NA,
                    nrow = (2^length(chosenFortunaFeats)) - 1,
                    ncol = 2,
                    dimnames = list(sprintf("combn.%s",
                                            1:((2^length(chosenFortunaFeats)) - 1)),
                                    c("keys", "formula")))
combnKeys[, "keys"] <- rownames(combnKeys)

combnNumber <- 1
rowIndex <- 1

for(i in 1:length(chosenFortunaFeats)){
  
  icombn <- combn(x = chosenFortunaFeats, m = i)
  cat(paste("\nNumber of features selected: ", as.character(i), ".\n", sep = ""))
  
  for(j in 1:ncol(icombn)){
    
    jfeatCombn <- icombn[, j]
    
    jformula <- paste(classColName,
                      " ~ ",
                      paste(jfeatCombn,
                            collapse = " + "),
                      sep = "")
    
    cat(paste("Formula selected: '", jformula, "'.\n", sep = ""))
    
    jcombnName <- paste("combn.", as.character(combnNumber), sep = "")
    
    combnKeys[jcombnName, "formula"] <- jformula
    
    jformula <- as.formula(jformula)
    
    for(k in 1:length(allAlgorithms)){
      
      kalgorithm <- allAlgorithms[k]
      
      cat(paste("Algorithm selected: '", kalgorithm, "'.\n", sep = ""))
      
      if(kalgorithm != "glmnet" |
         (kalgorithm == "glmnet" & length(jfeatCombn) > 1)){
        
        if(kalgorithm == "knn" | kalgorithm == "rpart"){
          set.seed(1)
          ktrain <- try(train(jformula,
                              data = trDataFrame,
                              trControl = trainControl(method = "LOOCV",
                                                       classProbs = TRUE),
                              method = kalgorithm))
        }else{
          set.seed(1)
          ktrain <- try(train(jformula,
                              data = trDataFrame,
                              trControl = trainControl(method = "LOOCV",
                                                       classProbs = TRUE),
                              method = kalgorithm,
                              ntree = 500,
                              hidden = ceiling(0.67 * length(jfeatCombn))))
        }
        
        if(class(ktrain) != "try-error"){
          
          modelIndex <- which.max(ktrain$results$Accuracy)[1]
          CVsetResultsShell[rowIndex,"accuracy"] <- ktrain$results[modelIndex, "Accuracy"]
          CVsetResultsShell[rowIndex,"kappa"] <- ktrain$results[modelIndex, "Kappa"]
          
          kconfusionMatrix <- confusionMatrix(data = predict(object = ktrain,
                                                             newdata = testData),
                                              reference = testClasses)
          
          shellForTest[rowIndex, "accuracy"] <- kconfusionMatrix$overall["Accuracy"]
          shellForTest[rowIndex, "kappa"] <- kconfusionMatrix$overall["Kappa"]
          shellForTest[rowIndex, "caret.pvalue"] <- kconfusionMatrix$overall["AccuracyPValue"]
          
        }
        
      }
      
      rowIndex <- rowIndex + 1
      
      save(x = CVsetResultsShell,
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/CVsetResultsShell2.rda")
      save(x = shellForTest,
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/shellForTest2.rda")
      
    }
    
    combnNumber <- combnNumber + 1
    
  }
  
}

View(CVsetResultsShell)
View(shellForTest)

CVnoNA <- na.omit(CVsetResultsShell)
shellNoNA <- na.omit(shellForTest)

algorithmColors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")
names(algorithmColors) <- allAlgorithms

xConstantRange <- c(1, 200, 400, 600, 800, 1023)
yConstantRange <- pretty(c(CVnoNA[, "accuracy"], shellNoNA[, "accuracy"]))

for(i in 1:length(allAlgorithms[-which(allAlgorithms == "rpart")])){
  
  ialgorithm <- allAlgorithms[-which(allAlgorithms == "rpart")][i]
  
  ibool <- grepl(pattern = ialgorithm, x = rownames(CVnoNA))
  
  iCVnoNA <- CVnoNA[ibool, "accuracy"]
  ix <- 1:length(iCVnoNA)
  
  ibool <- grepl(pattern = ialgorithm, x = rownames(shellNoNA))
  ishellNoNA <- shellNoNA[ibool, "accuracy"]
  
  imain <- paste(completeAlgorithmsName[allAlgorithms == ialgorithm],
                 " cross-validation set",
                 " accuracy plot over feature selection",
                 sep = "")
  
  icol <- rep(algorithmColors[ialgorithm], length(iCVnoNA))
  
  plot(x = c(ix),
       y = c(iCVnoNA),
       xlim = c(xConstantRange[1], xConstantRange[length(xConstantRange)]),
       ylim = range(yConstantRange),
       xaxt = "n",
       yaxt = "n",
       main = imain,
       cex.main = 1.5,
       xlab = "Feature combination",
       ylab = "Accuracy (proportion of correctly classified samples)",
       cex.lab = 1.5,
       type = "l",
       lty = 1,
       lwd = 1.5,
       col = icol)
  axis(side = 1, at = xConstantRange, labels = xConstantRange, las = 1, cex.axis = 1.3)
  axis(side = 2, at = yConstantRange, labels = yConstantRange, las = 1, cex.axis = 1.3)
  
  icol <- rep(colorspace::darken(col = algorithmColors[ialgorithm], amount = 0.4),
              length(iCVnoNA))
  
  imain <- paste(completeAlgorithmsName[allAlgorithms == ialgorithm],
                 " test set",
                 " accuracy plot over feature selection",
                 sep = "")
  
  plot(x = c(ix),
       y = c(ishellNoNA),
       xlim = c(xConstantRange[1], xConstantRange[length(xConstantRange)]),
       ylim = range(yConstantRange),
       xaxt = "n",
       yaxt = "n",
       main = imain,
       cex.main = 1.5,
       xlab = "Feature combination",
       ylab = "Accuracy (proportion of correctly classified samples)",
       cex.lab = 1.5,
       type = "l",
       lty = 1,
       lwd = 1.5,
       col = icol)
  axis(side = 1, at = xConstantRange, labels = xConstantRange, las = 1, cex.axis = 1.3)
  axis(side = 2, at = yConstantRange, labels = yConstantRange, las = 1, cex.axis = 1.3)
  
}

maxTestAcc <- shellNoNA[max(shellNoNA[, "accuracy"]) == shellNoNA[, "accuracy"],]
maxTestAcc <- maxTestAcc[max(maxTestAcc[, "kappa"]) == maxTestAcc[, "kappa"],]
equivalentCV <- CVnoNA[rownames(maxTestAcc),]
algAndCombn <- rownames(equivalentCV)[which.max(equivalentCV[, "accuracy"])]
algAndCombn <- unlist(strsplit(x = algAndCombn, split = ".", fixed = TRUE))
selectedAlg <- algAndCombn[1]
selectedCombn <- paste(algAndCombn[2], ".", algAndCombn[3], sep = "")
selectedFormula <- combnKeys[selectedCombn, "formula"]
selectedFeatures <- unlist(strsplit(x = unlist(strsplit(x = selectedFormula,
                                                        split = " ~ ",
                                                        fixed = TRUE))[2],
                                    split = " + ",
                                    fixed = TRUE))

nBestAlgCombn <- 20

newRowNames <- rep(allAlgorithms,
                   each = nBestAlgCombn)
newRowNames <- paste(newRowNames, ".%s", sep = "")
newRowNames <- sprintf(newRowNames,
                       rep(1:20,
                           length(allAlgorithms)))

colorMatrixCV <- matrix(data = NA,
                        nrow = nBestAlgCombn * length(allAlgorithms),
                        ncol = length(chosenFortunaFeats),
                        dimnames = list(newRowNames, chosenFortunaFeats))
colorMatrixTest <- colorMatrixCV

colorMatrixIndex <- 1

for(i in 1:length(allAlgorithms)){
  
  ialgorithm <- allAlgorithms[i]
  
  iCVbool <- grepl(pattern = ialgorithm, x = rownames(CVnoNA), fixed = TRUE)
  ishellBool <- grepl(pattern = ialgorithm, x = rownames(shellNoNA), fixed = TRUE)
  
  iCV <- CVnoNA[iCVbool,]
  iCVrows <- rownames(iCV[order(iCV[, "accuracy"], decreasing = TRUE), ])[1:nBestAlgCombn]
  iCVcombns <- unlist(lapply(X = strsplit(x = iCVrows, split = ".", fixed = TRUE),
                             FUN = function(x) paste(x[2], ".", x[3], sep = "")))
  iCVfeats <- sapply(X = unlist(lapply(X = strsplit(x = combnKeys[iCVcombns, "formula"],
                                                    split = " ~ ",
                                                    fixed = TRUE),
                                       FUN = function(x) x[-1])),
                     FUN = function(x) strsplit(x = x,
                                                split = " + ",
                                                fixed = TRUE))
  
  ishell <- shellNoNA[ishellBool,]
  ishellRows <- rownames(ishell[order(ishell[, "accuracy"], decreasing = TRUE), ])[1:nBestAlgCombn]
  ishellCombns <- unlist(lapply(X = strsplit(x = ishellRows, split = ".", fixed = TRUE),
                                FUN = function(x) paste(x[2], ".", x[3], sep = "")))
  ishellFeats <- sapply(X = unlist(lapply(X = strsplit(x = combnKeys[ishellCombns, "formula"],
                                                       split = " ~ ",
                                                       fixed = TRUE),
                                          FUN = function(x) x[-1])),
                        FUN = function(x) strsplit(x = x,
                                                   split = " + ",
                                                   fixed = TRUE))
  
  for(j in 1:length(iCVfeats)){
    
    colorMatrixCV[colorMatrixIndex, ] <- as.numeric(colnames(colorMatrixCV) %in% iCVfeats[[j]])
    colorMatrixTest[colorMatrixIndex, ] <- as.numeric(colnames(colorMatrixTest) %in% ishellFeats[[j]])
    
    colorMatrixIndex <- colorMatrixIndex + 1
    
  }
  
}

library(plot.matrix)

par(mar = c(5.1, 4.1, 4.1, 2.1)) # Default values

plot(x = colorMatrixCV,
     col = c("royalblue3", "red"),
     # key = NULL,
     na.col = "black",
     main = paste("Cross-validation set best feature combination",
                  sep = ""),
     xlab = "",
     ylab = "",
     axis.row = list(cex.axis = 0.55, las = 1),
     axis.col = list(cex.axis = 1, las = 2),
     key = NULL)

plot(x = colorMatrixTest,
     col = c("royalblue3", "red"),
     # key = NULL,
     na.col = "black",
     main = paste("Test set best feature combination",
                  sep = ""),
     xlab = "",
     ylab = "",
     axis.row = list(cex.axis = 0.55, las = 1),
     axis.col = list(cex.axis = 1, las = 2),
     key = NULL)

summary(lm(formula = as.formula("y ~ x"),
           data = cbind.data.frame("x" = CVnoNA[, "kappa"],
                                   "y" = shellNoNA[, "kappa"])))

set.seed(1)
ldaFinalModel <- train(as.formula(selectedFormula),
                       data = trDataFrame,
                       method = "lda",
                       trControl = trainControl(method = "LOOCV",
                                                classProbs = TRUE))

confusionMatrix(data = predict(object = ldaFinalModel,
                               newdata = testData[, selectedFeatures]),
                reference = testClasses)

trDataPrediction <- predict(object = ldaFinalModel$finalModel,
                            newdata = trData[, selectedFeatures])

testDataPrediction <- predict(object = ldaFinalModel$finalModel,
                              newdata = testData[, selectedFeatures])

allLDACoordPred <- rbind(trDataPrediction$x, testDataPrediction$x)

firstLD <- 1
secondLD <- 2

svdReady <- ldaFinalModel$finalModel$svd^2

xLabel <- paste("LD",
                as.character(firstLD),
                " (",
                as.character(round(100 * (svdReady[firstLD] / sum(svdReady)), 2)),
                " %)",
                sep = "")
yLabel <- paste("LD",
                as.character(secondLD),
                " (",
                as.character(round(100 * (svdReady[secondLD] / sum(svdReady)), 2)),
                " %)",
                sep = "")

plot(x = allLDACoordPred[, firstLD],
     y = allLDACoordPred[, secondLD],
     xlim = range(pretty(allLDACoordPred[, firstLD])),
     ylim = range(pretty(allLDACoordPred[, secondLD])),
     xaxt = "n",
     yaxt = "n",
     col = "black",
     pch = 21,
     bg = covidColors[c(as.character(trClasses), as.character(testClasses))],
     cex = 2,
     main = "LDA (optimized feature selection model)",
     cex.main = 2,
     xlab = xLabel,
     ylab = yLabel,
     cex.lab = 1.5)
axis(side = 1,
     at = pretty(allLDACoordPred[, firstLD]),
     labels = pretty(allLDACoordPred[, firstLD]),
     las = 1,
     cex.axis = 1.5)
axis(side = 2,
     at = pretty(allLDACoordPred[, secondLD]),
     labels = pretty(allLDACoordPred[, secondLD]),
     las = 1,
     cex.axis = 1.5)

legend(x = "topright",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.3,
       pt.cex = 1.75,
       inset = c(0.02, 0.05))

trainingLDASVM <- cbind.data.frame(trDataPrediction$x,
                                   trClasses)
colnames(trainingLDASVM)[ncol(trainingLDASVM)] <- classColName

newCombinedModelFormula <- as.formula(paste(classColName,
                                            " ~ ",
                                            paste(colnames(trainingLDASVM)[-ncol(trainingLDASVM)],
                                                  collapse = " + ")))

LDASVMcombinedModel <- train(newCombinedModelFormula,
                             trainingLDASVM,
                             method = "svmRadial",
                             trControl = trainControl(method = "LOOCV",
                                                      classProbs = TRUE))

confusionMatrix(data = predict(object = LDASVMcombinedModel,
                               newdata = testDataPrediction$x),
                reference = testClasses)

################################################################################

save(x = scaledCovidData,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/scaledCovidData.rda")
save(x = patientsKey,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/patientsKey.rda")

set.seed(1)
completeDataFortuna <- fortuna(dataset = scaledCovidData,
                               groups = as.factor(patientsKey[,"Group"]), # as.factor(groups)
                               treeNumber = 1000,
                               retentionByProp = FALSE,
                               retentionProportion = 0.8,
                               retentionLoops = 30,
                               chosenIterationsNumber = 30,
                               # NAsubstitution = "none", # "none", "sampling", "zero"
                               isBinaryApproach = FALSE,
                               combBinaryClassMethod = "p", # "k", "elbow", "p"
                               chosenPValue = 0.05,
                               pAdjMethod = "fdr",
                               kFeats = 20,
                               isElbowPlot = TRUE)

View(completeDataFortuna$importance.matrix$impArr)

fullDataSelectedFeatures <- names(which(completeDataFortuna$importance.matrix$impArr[, "adjusted.pvalue"] <= 0.05))

###

firstPC <- 1
secondPC <- 2

fullDataSelectionPCA <- prcomp(x = scaledCovidData[, fullDataSelectedFeatures],
                               scale. = FALSE)

PCxLab <- paste("PC", as.character(firstPC), " (",
                as.character(round(100 * summary(fullDataSelectionPCA)$importance[2, firstPC],1)),
                " %)",
                sep = "")
PCyLab <- paste("PC", as.character(secondPC), " (",
                as.character(round(100 * summary(fullDataSelectionPCA)$importance[2, secondPC],1)),
                " %)",
                sep = "")

plot(x = fullDataSelectionPCA$x[,firstPC],
     y = fullDataSelectionPCA$x[,secondPC],
     xlim = range(pretty(fullDataSelectionPCA$x[,firstPC])),
     ylim = range(pretty(fullDataSelectionPCA$x[,secondPC])),
     xaxt = "n",
     yaxt = "n",
     xlab = PCxLab,
     ylab = PCyLab,
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "PCA (PRM covid experiment) (selected important features)",
     cex.lab = 1.5)
axis(side = 1,
     at = pretty(fullDataSelectionPCA$x[,firstPC]),
     labels = pretty(fullDataSelectionPCA$x[,firstPC]),
     las = 1,
     cex.axis = 1.75)
axis(side = 2,
     at = pretty(fullDataSelectionPCA$x[,secondPC]),
     labels = pretty(fullDataSelectionPCA$x[,secondPC]),
     las = 1,
     cex.axis = 1.75)

legend(x = "bottomright",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.4,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

allAlgorithms <- c("rf",
                   "lda",
                   "svmRadial",
                   "glmnet",
                   "knn",
                   "nb",
                   "nnet",
                   "rpart")
completeAlgorithmsName <- c("Random forest (RF)",
                            "Linear Discriminant Analysis (LDA)",
                            "Support Vector Machines (SVM) (Gaussian kernel)",
                            "Logistic regression",
                            "K-Nearest Neighbors (KNN)",
                            "Naive Bayes",
                            "Artificial Neural Network (ANN)",
                            "Decision tree")

performanceMeasurements <- c("accuracy", "kappa")

CVshell_rownames <- sprintf(rep(paste(allAlgorithms,
                                      ".combn.%s",
                                      sep = ""),
                                ((2^length(fullDataSelectedFeatures)) - 1)),
                            rep(1:((2^length(fullDataSelectedFeatures)) - 1),
                                each = length(allAlgorithms)))

CVshell <- matrix(data = NA,
                  nrow = length(CVshell_rownames),
                  ncol = length(performanceMeasurements),
                  dimnames = list(CVshell_rownames,
                                  performanceMeasurements))

fullSetCombnKeys <- matrix(data = NA,
                           nrow = (2^length(fullDataSelectedFeatures)) - 1,
                           ncol = 2,
                           dimnames = list(sprintf("combn.%s",
                                                   1:((2^length(fullDataSelectedFeatures)) - 1)),
                                           c("keys", "formula")))
fullSetCombnKeys[, "keys"] <- rownames(fullSetCombnKeys)

classColName <- "patient"
fullSetDataFrame <- cbind.data.frame(scaledCovidData,
                                     patientsKey[,"group"])
colnames(fullSetDataFrame)[ncol(fullSetDataFrame)] <- classColName

combnNumber <- 1
rowIndex <- 1

for(i in 1:length(fullDataSelectedFeatures)){
  
  icombn <- combn(x = fullDataSelectedFeatures, m = i)
  cat(paste("\nNumber of features selected: ", as.character(i), ".\n", sep = ""))
  
  for(j in 1:ncol(icombn)){
    
    jfeatCombn <- icombn[, j]
    
    jformula <- paste(classColName,
                      " ~ ",
                      paste(jfeatCombn,
                            collapse = " + "),
                      sep = "")
    
    cat(paste("Formula selected: '", jformula, "'.\n", sep = ""))
    
    jcombnName <- paste("combn.", as.character(combnNumber), sep = "")
    
    fullSetCombnKeys[jcombnName, "formula"] <- jformula
    
    jformula <- as.formula(jformula)
    
    for(k in 1:length(allAlgorithms)){
      
      kalgorithm <- allAlgorithms[k]
      
      cat(paste("Algorithm selected: '", kalgorithm, "'.\n", sep = ""))
      
      if(kalgorithm != "glmnet" |
         (kalgorithm == "glmnet" & length(jfeatCombn) > 1)){
        
        if(kalgorithm == "knn" | kalgorithm == "rpart"){
          set.seed(1)
          ktrain <- try(train(jformula,
                              data = fullSetDataFrame,
                              trControl = trainControl(method = "LOOCV",
                                                       classProbs = TRUE),
                              method = kalgorithm))
        }else{
          set.seed(1)
          ktrain <- try(train(jformula,
                              data = fullSetDataFrame,
                              trControl = trainControl(method = "LOOCV",
                                                       classProbs = TRUE),
                              method = kalgorithm,
                              ntree = 500,
                              hidden = ceiling(0.67 * length(jfeatCombn))))
        }
        
        if(class(ktrain) != "try-error"){
          
          modelIndex <- which.max(ktrain$results$Accuracy)[1]
          CVshell[rowIndex,"accuracy"] <- ktrain$results[modelIndex, "Accuracy"]
          CVshell[rowIndex,"kappa"] <- ktrain$results[modelIndex, "Kappa"]
          
        }
        
      }
      
      rowIndex <- rowIndex + 1
      
      save(x = CVshell,
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/CVshell.rda")
      
    }
    
    combnNumber <- combnNumber + 1
    
  }
  
}

View(CVshell)

maxCVAcc <- CVshell[max(CVshell[, "accuracy"]) == CVshell[, "accuracy"],]
maxCVAcc <- maxCVAcc[max(maxCVAcc[, "kappa"]) == maxCVAcc[, "kappa"],]
fullSequivalentCV <- CVshell[rownames(maxCVAcc),]
fullSAlgAndCombn <- rownames(equivalentCV)[which.max(fullSequivalentCV[, "accuracy"])]
fullSAlgAndCombn <- unlist(strsplit(x = fullSAlgAndCombn, split = ".", fixed = TRUE))
fullSAlg <- fullSAlgAndCombn[1]
fullSCombn <- paste(fullSAlgAndCombn[2], ".", fullSAlgAndCombn[3], sep = "")
fullSFormula <- combnKeys[fullSCombn, "formula"]
fullSFeatures <- unlist(strsplit(x = unlist(strsplit(x = fullSFormula,
                                                     split = " ~ ",
                                                     fixed = TRUE))[2],
                                 split = " + ",
                                 fixed = TRUE))

