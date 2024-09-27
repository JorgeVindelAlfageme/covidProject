setwd(dir = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/")

dir.create(path = "rda")

table1 <- read.csv2(file = "20230803_Norm_Age_Area_Proteins_Tometaboloanalyst.csv")
table2 <- read.csv2(file = "20230803_Norm_Age_Area_Proteins_Tometaboloanalyst_GroupsbyAge.csv")
table3 <- read.csv2(file = "20230803_Norm_Area_Proteins_Tometaboloanalyst.csv")

library(xlsx)

uniprotDescription <- read.xlsx2(file = "idmapping_2023_08_03.xlsx",
                                 sheetIndex = 1)
rownames(uniprotDescription) <- uniprotDescription[,"Entry"]
uniprotDescription <- uniprotDescription[, c(4, 5)]
colnames(uniprotDescription) <- c("function", "gene")

save(x = uniprotDescription, file = "rda/uniprotDescription.rda")

patientsKey <- cbind.data.frame(table1[, c("Sample", "Group", "Age")])
rownames(patientsKey) <- patientsKey[, "Sample"]
rownames(table2) <- table2[, 1]
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
names(covidColors) <- unique(patientsKey[, "Group"])
allSamplesColors <- covidColors[patientsKey[, "Group"]]

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

library(caret)

covidClasses <- as.factor(patientsKey[,"Group"])
partitionProp <- 0.7
set.seed(1)
rnorm(n = 1)
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
                        combBinaryClassMethod = "p", # "k", "elbow", "p"
                        chosenPValue = 0.05,
                        pAdjMethod = "fdr",
                        kFeats = 20,
                        isElbowPlot = TRUE)

fatimaDescription <- strsplit(x = colnames(table1)[-c(1:3)], split = ".", fixed = TRUE)
fatimaDescription <- data.frame("protein.code" = unlist(lapply(X = fatimaDescription,
                                                               FUN = function(x) x[2])),
                                "fatima.gene" = unlist(lapply(X = fatimaDescription,
                                                              FUN = function(x) x[3])))
rownames(fatimaDescription) <- fatimaDescription[,"protein.code"]
uniprotDescription <- cbind.data.frame(fatimaDescription,
                                       uniprotDescription[rownames(fatimaDescription),])
rownames(uniprotDescription) <- uniprotDescription[,"fatima.gene"]
save(x = uniprotDescription,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/uniprotDescription.rda")

library(xlsx)
write.xlsx(x = cbind.data.frame(uniprotDescription[rownames(covidFortuna$importance.matrix$impArr),],
                                covidFortuna$importance.matrix$impArr),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/allResults/fortunaTraining.xlsx",
           sheetName = "Sheet1",
           row.names = FALSE,
           col.names = TRUE,
           append = FALSE)

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

testChosenCombn <- rownames(shellNoNA)[shellNoNA[, "caret.pvalue"] <= 0.05]
kappaThres <- min(shellNoNA[testChosenCombn, "kappa"])
CVsigCombns <- rownames(CVnoNA)[CVnoNA[, "kappa"] >= kappaThres]
combnsIntersect <- intersect(x = testChosenCombn, y = CVsigCombns)

algAndCombn <- combnsIntersect[which.max(CVnoNA[combnsIntersect,])]
algAndCombn <- unlist(strsplit(x = algAndCombn, split = ".", fixed = TRUE))

selectedAlg <- algAndCombn[1]
selectedCombn <- paste(algAndCombn[2], ".", algAndCombn[3], sep = "")
selectedFormula <- combnKeys[selectedCombn, "formula"]
selectedFeatures <- unlist(strsplit(x = unlist(strsplit(x = selectedFormula,
                                                        split = " ~ ",
                                                        fixed = TRUE))[2],
                                    split = " + ",
                                    fixed = TRUE))

combnsForTestArr <- unlist(lapply(X = strsplit(x = testChosenCombn,
                                               split = ".",
                                               fixed = TRUE),
                                  FUN = function(x) paste(x[2], ".", x[3], sep = "")))
combnsForTestArr <- sapply(X = unlist(lapply(X = strsplit(x = combnKeys[combnsForTestArr, "formula"],
                                                          split = " ~ ",
                                                          fixed = TRUE),
                                             FUN = function(x) x[2])),
                           FUN = function(x) unlist(strsplit(x = x,
                                                             split = " + ",
                                                             fixed = TRUE)))

colorMatrixTest2 <- matrix(data = NA,
                           nrow = length(combnsForTestArr),
                           ncol = length(chosenFortunaFeats),
                           dimnames = list(unlist(lapply(X = strsplit(x = testChosenCombn,
                                                                      split = ".",
                                                                      fixed = TRUE),
                                                         FUN = function(x) paste(x[2], ".", x[3], sep = ""))),
                                           chosenFortunaFeats))

for(i in 1:nrow(colorMatrixTest2)){
  
  colorMatrixTest2[i,] <-
    as.numeric(colnames(colorMatrixTest2) %in% combnsForTestArr[[i]])
  
}

combnsForCVArr <- unlist(lapply(X = strsplit(x = CVsigCombns,
                                             split = ".",
                                             fixed = TRUE),
                                FUN = function(x) paste(x[2], ".", x[3], sep = "")))
combnsForCVArr <- sapply(X = unlist(lapply(X = strsplit(x = combnKeys[combnsForCVArr, "formula"],
                                                        split = " ~ ",
                                                        fixed = TRUE),
                                           FUN = function(x) x[2])),
                         FUN = function(x) unlist(strsplit(x = x,
                                                           split = " + ",
                                                           fixed = TRUE)))

colorMatrixCV2 <- matrix(data = NA,
                         nrow = length(combnsForCVArr),
                         ncol = length(chosenFortunaFeats),
                         dimnames = list(unlist(lapply(X = strsplit(x = CVsigCombns,
                                                                    split = ".",
                                                                    fixed = TRUE),
                                                       FUN = function(x) paste(x[2], ".", x[3], sep = ""))),
                                         chosenFortunaFeats))
for(i in 1:nrow(colorMatrixCV2)){
  
  colorMatrixCV2[i,] <-
    as.numeric(colnames(colorMatrixCV2) %in% combnsForCVArr[[i]])
  
}

library(plot.matrix)

par(mar = c(5.1, 4.1, 4.1, 2.1)) # Default values

plot(x = colorMatrixCV2,
     col = c("royalblue3", "red"),
     # key = NULL,
     na.col = "black",
     main = paste("Best cross-validation set significant feature combination",
                  sep = ""),
     xlab = "",
     ylab = "",
     axis.row = list(cex.axis = 0.55, las = 1),
     axis.col = list(cex.axis = 1, las = 2),
     key = NULL)

plot(x = colorMatrixTest2[1:200,],
     col = c("royalblue3", "red"),
     # key = NULL,
     na.col = "black",
     main = paste("Best test set significant feature combination",
                  sep = ""),
     xlab = "",
     ylab = "",
     axis.row = list(cex.axis = 0.55, las = 1),
     axis.col = list(cex.axis = 1, las = 2),
     key = NULL)

CVps2 <- apply(X = colorMatrixCV2,
               MARGIN = 2,
               FUN = function(x) binom.test(x = sum(x),
                                            n = length(x),
                                            p = sum(colorMatrixCV2) /
                                              prod(dim(colorMatrixCV2)),
                                            alternative = "greater")$p.value)
CVps2 <- CVps2[order(CVps2, decreasing = FALSE)]
adj_CVps2 <- p.adjust(p = CVps2, method = "BH")
colSums(colorMatrixCV2)[order(colSums(colorMatrixCV2), decreasing = TRUE)]

write.xlsx(x = cbind.data.frame(uniprotDescription[names(adj_CVps2),],
                                "times.chosen" = 
                                  colSums(colorMatrixCV2)[order(colSums(colorMatrixCV2),
                                                                decreasing = TRUE)],
                                "pvalue" = CVps2,
                                "adj.pvalue" = adj_CVps2),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/coloredMatrices/sig/coloredArrSigExcel.xlsx",
           row.names = FALSE,
           col.names = TRUE,
           sheetName = "CV",
           append = FALSE)

testps2 <- apply(X = colorMatrixTest2,
                 MARGIN = 2,
                 FUN = function(x) binom.test(x = sum(x),
                                              n = length(x),
                                              p = sum(colorMatrixTest2) /
                                                prod(dim(colorMatrixTest2)),
                                              alternative = "greater")$p.value)

testps2 <- testps2[order(testps2, decreasing = FALSE)]
adj_testps2 <- p.adjust(p = testps2, method = "BH")

write.xlsx(x = cbind.data.frame(uniprotDescription[names(adj_testps2),],
                                "times.chosen" = 
                                  colSums(colorMatrixTest2)[order(colSums(colorMatrixTest2),
                                                                  decreasing = TRUE)],
                                "pvalue" = testps2,
                                "adj.pvalue" = adj_testps2),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/coloredMatrices/sig/coloredArrSigExcel.xlsx",
           row.names = FALSE,
           col.names = TRUE,
           sheetName = "test",
           append = TRUE)

################################################################################

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

CVps <- apply(X = colorMatrixCV,
              MARGIN = 2,
              FUN = function(x) binom.test(x = sum(x),
                                           n = length(x),
                                           p = sum(colorMatrixCV) /
                                             prod(dim(colorMatrixCV)),
                                           alternative = "greater")$p.value)
CVps <- CVps[order(CVps, decreasing = FALSE)]
adj_CVps <- p.adjust(p = CVps, method = "BH")
colSums(colorMatrixCV)[order(colSums(colorMatrixCV), decreasing = TRUE)]

write.xlsx(x = cbind.data.frame(uniprotDescription[names(adj_CVps),],
                                "times.chosen" = 
                                  colSums(colorMatrixCV)[order(colSums(colorMatrixCV),
                                                               decreasing = TRUE)],
                                "pvalue" = CVps,
                                "adj.pvalue" = adj_CVps),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/coloredMatrices/thres/coloredArrThresExcel.xlsx",
           row.names = FALSE,
           col.names = TRUE,
           sheetName = "CV",
           append = FALSE)

testps <- apply(X = colorMatrixTest,
                MARGIN = 2,
                FUN = function(x) binom.test(x = sum(x),
                                             n = length(x),
                                             p = sum(colorMatrixTest) /
                                               prod(dim(colorMatrixTest)),
                                             alternative = "greater")$p.value)

testps <- testps[order(testps, decreasing = FALSE)]
adj_testps <- p.adjust(p = testps, method = "BH")
colSums(colorMatrixTest)[order(colSums(colorMatrixTest), decreasing = TRUE)]

write.xlsx(x = cbind.data.frame(uniprotDescription[names(adj_testps),],
                                "times.chosen" = 
                                  colSums(colorMatrixTest)[order(colSums(colorMatrixTest),
                                                                 decreasing = TRUE)],
                                "pvalue" = testps,
                                "adj.pvalue" = adj_testps),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/coloredMatrices/thres/coloredArrThresExcel.xlsx",
           row.names = FALSE,
           col.names = TRUE,
           sheetName = "test",
           append = TRUE)

summary(lm(formula = as.formula("y ~ x"),
           data = cbind.data.frame("x" = CVnoNA[, "kappa"],
                                   "y" = shellNoNA[, "kappa"])))

set.seed(1)
finalModel <- train(as.formula(selectedFormula),
                    data = trDataFrame,
                    method = selectedAlg,
                    trControl = trainControl(method = "LOOCV",
                                             classProbs = TRUE))

confusionMatrix(data = predict(object = finalModel,
                               newdata = testData[, selectedFeatures]),
                reference = testClasses)

set.seed(1)
ldaFinalModel <- train(as.formula(selectedFormula),
                       data = trDataFrame,
                       method = "lda",
                       trControl = trainControl(method = "LOOCV",
                                                classProbs = TRUE))

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
       inset = c(0.02, 0.02))

set.seed(1)
optimizedGLMNETmodel <- train(as.formula(selectedFormula),
                              data = trDataFrame,
                              method = "glmnet",
                              family = "multinomial",
                              type.multinomial = "grouped",
                              tuneGrid = expand.grid(alpha = seq(from = 0.9,
                                                                 to = 2,
                                                                 by = 0.01),
                                                     lambda = seq(from = 0.004,
                                                                  to = 0.01,
                                                                  by = 0.0001)),
                              trControl = trainControl(method = "LOOCV",
                                                       classProbs = TRUE))

save(x = optimizedGLMNETmodel,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/optimizedGLMNETmodel.rda")

confusionMatrix(data = predict(object = optimizedGLMNETmodel,
                               newdata = testData[, selectedFeatures]),
                reference = testClasses)

set.seed(1)
optimizedNNETmodel <- train(as.formula(selectedFormula),
                            data = trDataFrame,
                            method = "nnet",
                            tuneGrid = expand.grid(size = 1:10,
                                                   decay = seq(from = 0,
                                                               to = 0.1,
                                                               by = 0.0001)),
                            trControl = trainControl(method = "LOOCV",
                                                     classProbs = TRUE))

save(x = optimizedNNETmodel,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/optimizedNNETmodel.rda")

confusionMatrix(data = predict(object = optimizedNNETmodel,
                               newdata = testData[, selectedFeatures]),
                reference = testClasses)

install.packages(c("caTools", "LiblineaR", "stepPlr", "RWeka"))
library(caTools)

set.seed(1)
model <- train(x = trData[, selectedFeatures],
               y = trClasses,
               method = "regLogistic",
               trControl = trainControl(method = "LOOCV",
                                        classProbs = TRUE))
confusionMatrix(data = predict(object = model,
                               newdata = testData[, selectedFeatures]),
                reference = testClasses)

firstPC <- 1
secondPC <- 2

finalModelPCA <- prcomp(x = scaledCovidData[, selectedFeatures],
                        scale. = FALSE)

PCxLab <- paste("PC", as.character(firstPC), " (",
                as.character(round(100 * summary(finalModelPCA)$importance[2, firstPC],1)),
                " %)",
                sep = "")
PCyLab <- paste("PC", as.character(secondPC), " (",
                as.character(round(100 * summary(finalModelPCA)$importance[2, secondPC],1)),
                " %)",
                sep = "")

plot(x = finalModelPCA$x[,firstPC],
     y = finalModelPCA$x[,secondPC],
     xlim = range(pretty(finalModelPCA$x[,firstPC])),
     ylim = range(pretty(finalModelPCA$x[,secondPC])),
     xaxt = "n",
     yaxt = "n",
     xlab = PCxLab,
     ylab = PCyLab,
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "PCA (PRM covid experiment) (final model features)",
     cex.lab = 1.5)
axis(side = 1,
     at = pretty(finalModelPCA$x[,firstPC]),
     labels = pretty(finalModelPCA$x[,firstPC]),
     las = 1,
     cex.axis = 1.75)
axis(side = 2,
     at = pretty(finalModelPCA$x[,secondPC]),
     labels = pretty(finalModelPCA$x[,secondPC]),
     las = 1,
     cex.axis = 1.75)

legend(x = "topleft",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.35,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

library(MASS)

allLDA <- lda(x = scaledCovidData, grouping = covidClasses)
importantFeatsLDA <- lda(x = scaledCovidData[, chosenFortunaFeats],
                         grouping = covidClasses)
finalLDA <- lda(x = scaledCovidData[, selectedFeatures], grouping = covidClasses)

firstLD <- 1
secondLD <- 2

xLDA <- paste("LD", as.character(firstLD), " (",
              as.character(round(100 * ((finalLDA$svd^2)[firstLD] /
                                          sum(finalLDA$svd^2)), 2)),
              " %)",
              sep = "")
yLDA <- paste("LD", as.character(secondLD), " (",
              as.character(round(100 * ((finalLDA$svd^2)[secondLD] /
                                          sum(finalLDA$svd^2)), 2)),
              " %)",
              sep = "")

ldaPredict <- predict(object = finalLDA, newdata = scaledCovidData[, selectedFeatures])

plot(x = ldaPredict$x[, firstLD],
     y = ldaPredict$x[, secondLD],
     xlim = range(pretty(ldaPredict$x[, firstLD])),
     ylim = range(pretty(ldaPredict$x[, secondLD])),
     xaxt = "n",
     yaxt = "n",
     xlab = xLDA,
     ylab = yLDA,
     pch = 21,
     col = "black",
     bg = allSamplesColors,
     cex = 2,
     cex.main = 2,
     main = "LDA (PRM covid experiment) (final model features)",
     cex.lab = 1.5)
axis(side = 1,
     at = pretty(ldaPredict$x[, firstLD]),
     labels = pretty(ldaPredict$x[, firstLD]),
     las = 1,
     cex.axis = 1.75)
axis(side = 2,
     at = pretty(ldaPredict$x[, secondLD]),
     labels = pretty(ldaPredict$x[, secondLD]),
     las = 1,
     cex.axis = 1.75)

legend(x = "bottomright",
       legend = c("Non-hospitalized (NHOSP)",
                  "Hospitalized (HOSP)",
                  "Intensive Care Unit ('UCI')",
                  "Deceased (EXI)"),
       pch = 21,
       pt.bg = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
       cex = 1.2,
       pt.cex = 1.75,
       inset = c(0.02, 0.02))

rawLDA <- lda(x = trData[, selectedFeatures], grouping = trClasses)
confusionMatrix(data = predict(object = rawLDA,
                               newdata = testData[, selectedFeatures])$class,
                reference = testClasses)

set.seed(1)
lastResort <- train(as.formula("patient ~ LD1 + LD2 + LD3"),
                    data = cbind.data.frame(predict(object = rawLDA,
                                                    newdata = trData[, selectedFeatures])$x,
                                            "patient" = trClasses),
                    method = "svmRadial",
                    trControl = trainControl(method = "LOOCV", classProbs = TRUE))
confusionMatrix(data = predict(object = lastResort,
                               newdata = predict(object = rawLDA,
                                                 newdata = testData[, selectedFeatures])$x),
                reference = testClasses)


observedROC <- c()
predictedROC <- c()

testPredictionProbs <- predict(object = finalModel, newdata = testData, type = "prob")

for(i in 1:nlevels(trClasses)){
  
  iclass <- levels(trClasses)[i]
  iobserved <- as.numeric(as.character(testClasses) == iclass)
  observedROC <- append(x = observedROC, values = iobserved)
  
  ipredicted <- testPredictionProbs[, iclass]
  predictedROC <- append(x = predictedROC, values = ipredicted)
  
}

finalModelROC <- verification::roc.plot(x = observedROC,
                                        pred = predictedROC,
                                        xlab =  "1 - specificity",
                                        ylab = "Sensitivity",
                                        show.thres = FALSE)

# https://stats.stackexchange.com/questions/2151/how-to-plot-roc-curves-in-multiclass-classification

################################################################################