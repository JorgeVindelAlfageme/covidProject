setwd(dir = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/")
setwd(dir = "C:/Users/usuario/Downloads/fatPaper/partOfFornax/fatima/20230803_PRM_HUMS/")

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

write.xlsx(x = cbind.data.frame(uniprotDescription[rownames(completeDataFortuna$importance.matrix$impArr),],
                                completeDataFortuna$importance.matrix$impArr),
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/allResults/fortuna2.xlsx",
           sheetName = "Sheet1",
           row.names = FALSE,
           col.names = TRUE,
           append = FALSE)

fullDataSelectedFeatures <- names(which(completeDataFortuna$importance.matrix$impArr[, "adj.pvalue"] <= 0.05))

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

set.seed(1)
wrapperLinearRegression <- wrapperCholestasis2(startData = trDataFrame,
                                               classColumn = classColName,
                                               defaultGrid = NULL,
                                               crossValidation = "LOOCV",
                                               bootSamples = 100,
                                               chosenAlgorithm = "glmnet",
                                               RFchosenNTree = 1000)

confusionMatrix(data = predict(object = wrapperLinearRegression$finalClassifier,
                               newdata = testData),
                reference = testClasses)

set.seed(1)
wrapperLDA <- wrapperCholestasis2(startData = trDataFrame,
                                  classColumn = classColName,
                                  defaultGrid = NULL,
                                  crossValidation = "LOOCV",
                                  bootSamples = 100,
                                  chosenAlgorithm = "lda",
                                  RFchosenNTree = 1000)

confusionMatrix(data = predict(object = wrapperLDA$finalClassifier,
                               newdata = testData),
                reference = testClasses)


set.seed(1)
noTestwrapperLinearRegression <- wrapperCholestasis2(startData = cbind.data.frame(scaledCovidData,
                                                                                  "patient" = covidClasses),
                                                     classColumn = classColName,
                                                     defaultGrid = NULL,
                                                     crossValidation = "LOOCV",
                                                     bootSamples = 100,
                                                     chosenAlgorithm = "glmnet",
                                                     RFchosenNTree = 1000)

set.seed(1)
noTestwrapperLDA <- wrapperCholestasis2(startData = cbind.data.frame(scaledCovidData,
                                                                     "patient" = covidClasses),
                                        classColumn = classColName,
                                        defaultGrid = NULL,
                                        crossValidation = "LOOCV",
                                        bootSamples = 100,
                                        chosenAlgorithm = "lda",
                                        RFchosenNTree = 1000)

noTestLDA_allFeats <- lda(x = scaledCovidData, grouping = covidClasses)
predict(object = noTestLDA_allFeats,
        newdata = scaledCovidData)$x[, c("LD1", "LD2")]
cbind.data.frame(predict(object = noTestLDA_allFeats,
                         newdata = scaledCovidData)$x[, c("LD1", "LD2")],
                 "patient" = covidClasses)
set.seed(1)
LDA2D_SVM_LOOCV <- train(as.formula(paste(classColName, " ~ ", "LD1 + ", "LD2",
                                          sep = "")),
                         data = cbind.data.frame(predict(object = noTestLDA_allFeats,
                                                         newdata = scaledCovidData)$x[, c("LD1", "LD2")],
                                                 "patient" = covidClasses),
                         method = "svmRadial",
                         trControl = trainControl(method = "LOOCV",
                                                  classProbs = TRUE))
LDA2D_SVM_LOOCV <- train(as.formula(paste(classColName, " ~ ", "LD1 + ", "LD2",
                                          sep = "")),
                         data = cbind.data.frame(predict(object = noTestLDA_allFeats,
                                                         newdata = scaledCovidData)$x[, c("LD1", "LD2")],
                                                 "patient" = covidClasses),
                         method = "svmRadial",
                         trControl = trainControl(method = "LOOCV",
                                                  classProbs = TRUE),
                         tuneGrid = expand.grid(C = seq(from = 0.01, to = 1, by = 0.01),
                                                sigma = seq(from = 0.01, to = 2, by = 0.02)))

noTestLDA_wrapperFeats <- lda(x = scaledCovidData[, noTestwrapperLDA$bestFeats],
                              grouping = covidClasses)
predict(object = noTestLDA_wrapperFeats,
        newdata = scaledCovidData[, noTestwrapperLDA$bestFeats])$x[, c("LD1", "LD2")]
cbind.data.frame(predict(object = noTestLDA_wrapperFeats,
                         newdata = scaledCovidData[, noTestwrapperLDA$bestFeats])$x[, c("LD1", "LD2")],
                 "patient" = covidClasses)
set.seed(1)
LDA2D_SVM_LOOCV_wrapper <- train(as.formula(paste(classColName, " ~ ", "LD1 + ", "LD2",
                                                  sep = "")),
                                 data = cbind.data.frame(predict(object = noTestLDA_wrapperFeats,
                                                                 newdata = scaledCovidData[, noTestwrapperLDA$bestFeats])$x[, c("LD1", "LD2")],
                                                         "patient" = covidClasses),
                                 method = "svmRadial",
                                 trControl = trainControl(method = "LOOCV",
                                                          classProbs = TRUE))

set.seed(1)
LDA3D_SVM_LOOCV <- train(as.formula(paste(classColName, " ~ ", ".",
                                          sep = "")),
                         data = cbind.data.frame(predict(object = noTestLDA_allFeats,
                                                         newdata = scaledCovidData)$x,
                                                 "patient" = covidClasses),
                         method = "svmRadial",
                         trControl = trainControl(method = "LOOCV",
                                                  classProbs = TRUE),
                         tuneGrid = expand.grid(C = seq(from = 0.01, to = 1, by = 0.01),
                                                sigma = seq(from = 0.01, to = 2, by = 0.02)))

oldAlgorithms <- c("rf", "lda", "nb", "knn", "svmRadial", "glmnet")

allWrapperResults_LOOCV <- list()
allWrapperResults_boot <- list()

for(i in 1:length(oldAlgorithms)){
  
  ialg <- oldAlgorithms[i]
  
  set.seed(1)
  allWrapperResults_LOOCV[[ialg]] <- wrapperCholestasis2(startData = trDataFrame,
                                                         classColumn = classColName,
                                                         defaultGrid = NULL,
                                                         crossValidation = "LOOCV",
                                                         bootSamples = 100,
                                                         chosenAlgorithm = ialg,
                                                         RFchosenNTree = 1000)
  
  iconfLOOCV <- confusionMatrix(data = predict(object = allWrapperResults_LOOCV[[ialg]][["finalClassifier"]],
                                               newdata = testData[, allWrapperResults_LOOCV[[ialg]][["bestFeats"]]]),
                                reference = testClasses)
  
  allWrapperResults_LOOCV[[ialg]][["confArr"]] <- iconfLOOCV
  
  save(x = allWrapperResults_LOOCV,
       file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/allWrapperResults_LOOCV.rda")
  
  set.seed(1)
  allWrapperResults_boot[[ialg]] <- wrapperCholestasis2(startData = trDataFrame,
                                                        classColumn = classColName,
                                                        defaultGrid = NULL,
                                                        crossValidation = "boot",
                                                        bootSamples = 100,
                                                        chosenAlgorithm = ialg,
                                                        RFchosenNTree = 1000)
  
  iconfBoot <- confusionMatrix(data = predict(object = allWrapperResults_boot[[ialg]][["finalClassifier"]],
                                              newdata = testData[, allWrapperResults_boot[[ialg]][["bestFeats"]]]),
                               reference = testClasses)
  
  allWrapperResults_boot[[ialg]][["confArr"]] <- iconfBoot
  
  save(x = allWrapperResults_boot,
       file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/allWrapperResults_boot.rda")
  
}

featureSelectionShell <- matrix(data = NA,
                                nrow = length(allWrapperResults_LOOCV) +
                                  length(allWrapperResults_boot),
                                ncol = ncol(scaledCovidData))
colnames(featureSelectionShell) <- colnames(scaledCovidData)
krow <- 1

FSshellRownames <- c()

for(i in 1:length(allWrapperResults_LOOCV)){
  
  ialgLOOCV <- names(allWrapperResults_LOOCV)[i]
  
  ifeatsLOOCV <- lapply(X = allWrapperResults_LOOCV,
                        FUN = function(x) x[["bestFeats"]])[[ialgLOOCV]]
  
  featureSelectionShell[krow, ] <-
    as.numeric(colnames(featureSelectionShell) %in% ifeatsLOOCV)
  
  if(ialgLOOCV == "glmnet"){
    ialgLOOCV <- "lr"
  }
  if(ialgLOOCV == "svmRadial"){
    ialgLOOCV <- "svm"
  }
  
  FSshellRownames <- append(x = FSshellRownames,
                            values = paste(ialgLOOCV, ".", "LOOCV", sep = ""))
  
  krow <- krow + 1
  
  ###
  
  ialgBoot <- names(allWrapperResults_boot)[i]
  
  ifeatsBoot <- lapply(X = allWrapperResults_boot,
                       FUN = function(x) x[["bestFeats"]])[[ialgBoot]]
  
  featureSelectionShell[krow, ] <-
    as.numeric(colnames(featureSelectionShell) %in% ifeatsBoot)
  
  if(ialgBoot == "glmnet"){
    ialgBoot <- "lr"
  }
  if(ialgBoot == "svmRadial"){
    ialgBoot <- "svm"
  }
  
  FSshellRownames <- append(x = FSshellRownames,
                            values = paste(ialgBoot, ".", "Boot", sep = ""))
  
  krow <- krow + 1
  
}

rownames(featureSelectionShell) <- FSshellRownames

plot(x = featureSelectionShell,
     col = c("royalblue3", "red"),
     # key = NULL,
     na.col = "black",
     main = paste("Wrapper feature selection for classification algorithms",
                  sep = ""),
     xlab = "",
     ylab = "",
     axis.row = list(cex.axis = 0.7, las = 1),
     axis.col = list(cex.axis = 1, las = 2),
     key = NULL)

ps_wrapper <- apply(X = featureSelectionShell,
                    MARGIN = 2,
                    FUN = function(x) binom.test(x = sum(x),
                                                 n = length(x),
                                                 p = sum(featureSelectionShell) /
                                                   prod(dim(featureSelectionShell)),
                                                 alternative = "greater")$p.value)
ps_wrapper <- ps_wrapper[order(ps_wrapper, decreasing = FALSE)]
adj_ps_wrapper <- p.adjust(p = ps_wrapper, method = "BH")

wrapperChoiceResults <- matrix(data = c(ps_wrapper,
                                        adj_ps_wrapper,
                                        colSums(featureSelectionShell)[names(ps_wrapper)]),
                               ncol = 3)
rownames(wrapperChoiceResults) <- names(ps_wrapper)
colnames(wrapperChoiceResults) <- c("binomial.pvalue", "adj.pvalue", "frequency")

write.xlsx2(x = cbind.data.frame(uniprotDescription[rownames(wrapperChoiceResults),],
                                 wrapperChoiceResults),
            file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/wrappers/wrappersFeatures.xlsx",
            sheetName = "Sheet1",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)

allWrapperResultsMatrix <- matrix(data = NA,
                                  nrow = length(allWrapperResults_LOOCV) +
                                    length(allWrapperResults_boot),
                                  ncol = 5)
colnames(allWrapperResultsMatrix) <- c("acc.cv", "kappa.cv",
                                       "acc.test", "kappa.test", "acc.pvalue.test")
rowNamesWrapperResults <- c()
krow <- 1

for(i in 1:length(allWrapperResults_LOOCV)){
  
  ialgLOOCV <- names(allWrapperResults_LOOCV)[i]
  
  imaxAccLOOCV_CV_index <- which.max(allWrapperResults_LOOCV[[ialgLOOCV]]$finalClassifier$results[,"Accuracy"])
  imaxAccLOOCV_CV <- allWrapperResults_LOOCV[[ialgLOOCV]]$finalClassifier$results[imaxAccLOOCV_CV_index,
                                                                                  "Accuracy"]
  imaxKappaLOOCV_CV <- allWrapperResults_LOOCV[[ialgLOOCV]]$finalClassifier$results[imaxAccLOOCV_CV_index,
                                                                                    "Kappa"]
  iaccLOOCV_test <- allWrapperResults_LOOCV[[ialgLOOCV]]$confArr$overall["Accuracy"]
  ikappaLOOCV_test <- allWrapperResults_LOOCV[[ialgLOOCV]]$confArr$overall["Kappa"]
  ipvalueLOOCV_test <- allWrapperResults_LOOCV[[ialgLOOCV]]$confArr$overall["AccuracyPValue"]
  
  allWrapperResultsMatrix[krow, c("acc.cv", "kappa.cv",
                                  "acc.test", "kappa.test", "acc.pvalue.test")] <-
    c(imaxAccLOOCV_CV, imaxKappaLOOCV_CV,
      iaccLOOCV_test, ikappaLOOCV_test, ipvalueLOOCV_test)
  
  rowNamesWrapperResults <- append(x = rowNamesWrapperResults,
                                   values = paste(ialgLOOCV, ".", "LOOCV", sep = ""))
  
  krow <- krow + 1
  
  ialgBoot <- names(allWrapperResults_boot)[i]
  
  imaxAccBoot_CV_index <- which.max(allWrapperResults_boot[[ialgBoot]]$finalClassifier$results[,"Accuracy"])
  imaxAccBoot_CV <- allWrapperResults_boot[[ialgBoot]]$finalClassifier$results[imaxAccBoot_CV_index,
                                                                               "Accuracy"]
  imaxKappaBoot_CV <- allWrapperResults_boot[[ialgBoot]]$finalClassifier$results[imaxAccBoot_CV_index,
                                                                                 "Kappa"]
  iaccBoot_test <- allWrapperResults_boot[[ialgBoot]]$confArr$overall["Accuracy"]
  ikappaBoot_test <- allWrapperResults_boot[[ialgBoot]]$confArr$overall["Kappa"]
  ipvalueBoot_test <- allWrapperResults_boot[[ialgBoot]]$confArr$overall["AccuracyPValue"]
  
  allWrapperResultsMatrix[krow, c("acc.cv", "kappa.cv",
                                  "acc.test", "kappa.test", "acc.pvalue.test")] <-
    c(imaxAccBoot_CV, imaxKappaBoot_CV,
      iaccBoot_test, ikappaBoot_test, ipvalueBoot_test)
  
  rowNamesWrapperResults <- append(x = rowNamesWrapperResults,
                                   values = paste(ialgBoot, ".", "boot", sep = ""))
  
  krow <- krow + 1
  
}

rownames(allWrapperResultsMatrix) <- rowNamesWrapperResults

write.xlsx(x = allWrapperResultsMatrix,
           file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/wrappers/allWrapperResultsMatrix.xlsx",
           sheetName = "Sheet1",
           col.names = TRUE,
           row.names = TRUE,
           append = FALSE)

orderIWant <- c("NHOSP", "HOSP", "UCI", "EXI")

for(i in 1:ncol(scaledCovidData)){
  
  ifeat <- colnames(scaledCovidData)[i]
  
  idataFrame <- cbind.data.frame(scaledCovidData[, ifeat], covidClasses)
  colnames(idataFrame) <- c(ifeat, "class")
  
  idataFrame$class <- factor(idataFrame$class, levels = orderIWant)
  
  boxplot(idataFrame[,1] ~ idataFrame[,2],
          data = idataFrame,
          xaxt = "n",
          yaxt = "n",
          ylim = range(pretty(idataFrame[,1])),
          xlab = "",
          main = "",
          ylab = "",
          pch = 21,
          col = covidColors[orderIWant],
          bg = covidColors[orderIWant],
          outcex = 2,
          cex.lab = 1.8)
  axis(side = 1,
       at = 1:nlevels(idataFrame[,2]),
       labels = rep("", nlevels(idataFrame[,2])))
  axis(side = 2,
       las = 1,
       at = pretty(idataFrame[,1]),
       labels = pretty(idataFrame[,1]),
       cex.axis = 1.8)
  
  highest <- range(pretty(idataFrame[,1]))[2]
  lowest <- range(pretty(idataFrame[,1]))[1]
  difference <- highest - lowest
  labelsYPosition <- lowest - ((1/17) * difference)
  titleYPosition <- highest + ((1/10) * difference)
  xLabelPosition <- 0
  
  text(x = 1:nlevels(idataFrame[,2]),
       y = labelsYPosition,
       labels = orderIWant,
       adj = c(1, 1),
       xpd = NA,
       srt = 35,
       cex = 2)
  
  text(x = mean(1:nlevels(idataFrame[,2])),
       y = titleYPosition,
       labels = ifeat,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 4.5)
  text(x = xLabelPosition,
       y = mean(range(pretty(idataFrame[,1]))),
       labels = "Scaled protein abundance",
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 2)
  
}

covidCorArr <- corMatrix(dataMatrix = scaledCovidData,
                         colorRange = c("firebrick", "white", "forestgreen"),
                         corMethod = "cor",
                         adaptMargins = c(5.1, 4.1, 4.1, 3.1),
                         cexMain = 2,
                         mainTitle = "Correlation matrix between features",
                         fmtKey = "%.2f",
                         corTagsPos = 3,
                         corTagsOffset = 0.2,
                         corDigits = 2, 
                         corTagsCex = 0.8,
                         axisRowCex = 1,
                         axisColCex = 1,
                         argsLabelsCex = 0.5,
                         argsLabelsCorrection = 0.5)

allCombnByPairs <- t(combn(x = colnames(scaledCovidData), m = 2))

sigCombn <- c()
notSigCombn <- c()

for(i in 1:nrow(allCombnByPairs)){
  
  feat1 <- allCombnByPairs[i, 1]
  feat2 <- allCombnByPairs[i, 2]
  
  if(covidCorArr$pvalue[feat1, feat2] <= 0.05){
    sigCombn <- append(x = sigCombn, values = i)
  }else{
    notSigCombn <- append(x = notSigCombn, values = i)
  }
  
}

plotClusterComb(combResList = list("clusterSig" = allCombnByPairs[sigCombn,],
                                   "clusterNotSig" = allCombnByPairs[notSigCombn,]),
                dat = scaledCovidData,
                classVec = as.character(covidClasses),
                chosenCor = "pearson",
                namedColorVector = covidColors,
                cexPoint = 3,
                cexTitle = 2.75,
                cexAxesLabels = 1.8,
                cexAxesMarks = 2,
                legPositions = rep("topleft", nrow(allCombnByPairs)),
                isLeg = TRUE,
                legCex = 1.75,
                legInset = c(0.02, 0.02),
                setDir = "cor/cors",
                chosenWidth = 1200,
                chosenHeight = 800,
                legPointsCex = 2.2)

library(bnlearn)
library(Rgraphviz)

covidGS <- gs(x = as.data.frame(scaledCovidData))
dev.new()
options(repr.plot.width = 100, repr.plot.height = 100)
# changes text size. I believe its maximum value gets skewed by the longest
# feature name by its number of characters.
graphviz.plot(covidGS)

library(randomForest)

set.seed(1)
scaledCovidData_RF_all <- randomForest(x = scaledCovidData,
                                       y = covidClasses,
                                       proximity = TRUE,
                                       importance = TRUE,
                                       ntree = 1000)
scaledCovidData_RF_all$importance

RFClassEx(rfObject = scaledCovidData_RF_all,
          trueClassesOrder = c("NHOSP", "HOSP", "UCI", "EXI"),
          classesColors = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
          cexMain = 2,
          cexAxisNumbers = 1.1,
          cexPlotAxis = 1.2,
          tiltAngle = 40,
          cexLabels = 1.1,
          cexLeg = 1.6,
          legInset = c(0.02, 0.02),
          adjText = c(1, 1),
          yLabelsPosition = -0.001)

set.seed(1)
scaledCovidData_RF_Fortuna <- randomForest(x = scaledCovidData[, chosenFortunaFeats],
                                           y = covidClasses,
                                           proximity = TRUE,
                                           importance = TRUE,
                                           ntree = 1000)

RFClassEx(rfObject = scaledCovidData_RF_Fortuna,
          trueClassesOrder = c("NHOSP", "HOSP", "UCI", "EXI"),
          classesColors = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
          cexMain = 2,
          cexAxisNumbers = 1.1,
          cexPlotAxis = 1.1,
          tiltAngle = 35,
          cexLabels = 1.5,
          cexLeg = 1.75,
          legInset = c(0.02, 0.02),
          adjText = c(1, 1),
          yLabelsPosition = -0.005)

set.seed(1)
scaledCovidData_RF_selected <- randomForest(x = scaledCovidData[, selectedFeatures],
                                            y = covidClasses,
                                            proximity = TRUE,
                                            importance = TRUE,
                                            ntree = 1000)

RFClassEx(rfObject = scaledCovidData_RF_selected,
          trueClassesOrder = c("NHOSP", "HOSP", "UCI", "EXI"),
          classesColors = covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
          cexMain = 2,
          cexAxisNumbers = 1.2,
          cexPlotAxis = 1.1,
          tiltAngle = 30,
          cexLabels = 1.8,
          cexLeg = 1.5,
          legInset = c(0.02, 0.02),
          adjText = c(1, 1),
          yLabelsPosition = -0.01)

################################################################################

# Columna para la combinaci?n de prote?nas, de la precisi?n, del kappa, y para
# el test, incluir el p-valor de caret.

SFEcombns <- sum(ncol(scaledCovidData):1)

CVTableResults <- as.data.frame(matrix(data = NA,
                                       nrow = SFEcombns,
                                       ncol = length(c("combn",
                                                       "accuracy",
                                                       "kappa"))))
colnames(CVTableResults) <- c("combn", "accuracy", "kappa")

testTableResults <- as.data.frame(matrix(data = NA,
                                         nrow = SFEcombns,
                                         ncol = length(c("combn",
                                                         "accuracy",
                                                         "kappa",
                                                         "pvalue"))))
colnames(testTableResults) <- c("combn", "accuracy", "kappa", "pvalue")

CVTableAllSamplesResults <- CVTableResults

scaledCovidDataFrame <- cbind.data.frame(scaledCovidData, covidClasses)
colnames(scaledCovidDataFrame)[ncol(scaledCovidDataFrame)] <- classColName

consideredProteins <- colnames(scaledCovidData)

krow <- 1

pass <- TRUE

listIndex <- 1
roundsPreservedFeatures <- list()
maxAccuracyVec <- c()

while(length(consideredProteins) > 1){
  
  icombnAccuracy <- c()
  
  if(pass){
    icombn <- data.frame(consideredProteins)
  }else{
    icombn <- combn(x = consideredProteins,
                    m = (length(consideredProteins) - 1))
  }
  
  cat(paste("\nCurrent number of features: ",
            as.character(length(consideredProteins)),
            ".\n",
            sep = ""))
  
  for(j in 1:ncol(icombn)){
    
    jcombn <- icombn[, j]
    jformula <- paste(classColName,
                      " ~ ",
                      paste(jcombn, collapse = " + "),
                      sep = "")
    
    cat(paste("\nCurrently chosen formula: '", jformula, "'.\n", sep = ""))
    
    CVTableResults[krow, "combn"] <- jformula
    testTableResults[krow, "combn"] <- jformula
    
    jformula <- as.formula(jformula)
    
    set.seed(1)
    jlda <- train(jformula,
                  data = trDataFrame,
                  method = "lda",
                  trControl = trainControl(method = "LOOCV",
                                           classProbs = TRUE))
    
    jtrData <- trData[, jcombn, drop = FALSE]
    
    jxTrain <- predict(object = jlda$finalModel, newdata = jtrData)$x
    jtrainDataFrame <- cbind.data.frame(jxTrain, trClasses)
    colnames(jtrainDataFrame)[ncol(jtrainDataFrame)] <- classColName
    
    set.seed(1)
    jmodel <- train(as.formula(paste(classColName, " ~ ", ".", sep = "")),
                    data = jtrainDataFrame,
                    method = "rf",
                    trControl = trainControl(method = "LOOCV",
                                             classProbs = TRUE))
    
    jindexMaxAcc <- which.max(jmodel$results[, "Accuracy"])
    CVTableResults[krow, "accuracy"] <- jmodel$results[jindexMaxAcc, "Accuracy"]
    CVTableResults[krow, "kappa"] <- jmodel$results[jindexMaxAcc, "Kappa"]
    
    icombnAccuracy <- append(x = icombnAccuracy,
                             values = jmodel$results[jindexMaxAcc, "Accuracy"])
    
    jtestData <- testData[, jcombn, drop = FALSE]
    
    jxTest <- predict(object = jlda$finalModel, newdata = jtestData)$x
    
    jtestConfArr <- confusionMatrix(data = predict(object = jmodel,
                                                   newdata = jxTest),
                                    reference = testClasses)
    
    testTableResults[krow, "accuracy"] <- jtestConfArr$overall["Accuracy"]
    testTableResults[krow, "kappa"] <- jtestConfArr$overall["Kappa"]
    testTableResults[krow, "pvalue"] <- jtestConfArr$overall["AccuracyPValue"]
    
    krow <- krow + 1
    
  }
  
  if(pass){
    
    pass <- FALSE
    cat(paste("\n", sep = ""))
    
  }else{
    
    consideredProteins <- icombn[, which.max(icombnAccuracy)]
    
    cat(paste("\nPreserved features: ",
              paste(consideredProteins, collapse = ", "),
              ".\n",
              sep = ""))
    
  }
  
  maxAccuracyVec <- append(x = maxAccuracyVec,
                           values = icombnAccuracy[which.max(icombnAccuracy)])
  roundsPreservedFeatures[[listIndex]] <- consideredProteins
  
  listIndex <- listIndex + 1
  
}

SVMLDAsave <- list("CVTableResults" = CVTableResults,
                   "testTableResults" = testTableResults,
                   "maxAccuracyVec" = maxAccuracyVec,
                   "roundsPreservedFeatures" = roundsPreservedFeatures)

save(x = SVMLDAsave,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/SVMLDAsave.rda")

View(CVTableResults)
View(testTableResults)

SVMLDAbestTest <- which(SVMLDAsave$testTableResults$pvalue <= 0.05)
SVMLDAbestTest <-
  rownames(SVMLDAsave$CVTableResults[SVMLDAbestTest,])[which.max(SVMLDAsave$CVTableResults[SVMLDAbestTest,"accuracy"])]
SVMLDAbestTest <- as.numeric(SVMLDAbestTest)

SVMLDA_features <- unlist(strsplit(x = unlist(strsplit(x = CVTableResults[SVMLDAbestTest, "combn"],
                                                       split = " ~ ",
                                                       fixed = TRUE))[2],
                                   split = " + ",
                                   fixed = TRUE))

fortunaAndSVM_intersect <- intersect(x = chosenFortunaFeats, y = SVMLDA_features)

set.seed(1)
lda_SVMLDA <- train(as.formula(paste(classColName, " ~ .", sep = "")),
                    data = trDataFrame[, c(SVMLDA_features, classColName)],
                    method = "lda",
                    trControl = trainControl(method = "LOOCV",
                                             classProbs = TRUE))

trData_SVMLDA <- trData[, SVMLDA_features, drop = FALSE]

xTrain_SVMLDA <- predict(object = lda_SVMLDA$finalModel,
                         newdata = trData_SVMLDA)$x
trainDataFrame_SVMLDA <- cbind.data.frame(xTrain_SVMLDA, trClasses)
colnames(trainDataFrame_SVMLDA)[ncol(trainDataFrame_SVMLDA)] <- classColName

set.seed(1)
svm_SVMLDA <- train(as.formula(paste(classColName, " ~ ", ".", sep = "")),
                    data = trainDataFrame_SVMLDA,
                    method = "svmRadial",
                    trControl = trainControl(method = "LOOCV",
                                             classProbs = TRUE))

testData_SVMLDA <- testData[, SVMLDA_features, drop = FALSE]

xTest_SVMLDA <- predict(object = lda_SVMLDA$finalModel,
                        newdata = testData_SVMLDA)$x

confusionMatrix(data = predict(object = svm_SVMLDA,
                               newdata = xTest_SVMLDA),
                reference = testClasses)

tuneGridSVM <- expand.grid(C = seq(from = 0.01, to = 2, by = 0.02),
                           sigma = seq(from = 0.0001, to = 1, by = 0.01))

tuneGridResults <- expand.grid(C = seq(from = 0.01, to = 2, by = 0.02),
                               sigma = seq(from = 0.0001, to = 1, by = 0.01),
                               CVacc = NA,
                               CVkappa = NA,
                               testAcc = NA,
                               testKappa = NA,
                               pvalue = NA)

for(i in 1:nrow(tuneGridSVM)){
  
  set.seed(1)
  isvm <- train(as.formula(paste(classColName, " ~ ", ".", sep = "")),
                data = trainDataFrame_SVMLDA,
                method = "svmRadial",
                trControl = trainControl(method = "LOOCV",
                                         classProbs = TRUE),
                tuneGrid = expand.grid(C = tuneGridSVM[i, "C"],
                                       sigma = tuneGridSVM[i, "sigma"]))
  
  tuneGridResults[i, "CVacc"] <- isvm$results$Accuracy
  tuneGridResults[i, "CVkappa"] <- isvm$results$Kappa
  
  iconf <- confusionMatrix(data = predict(object = isvm,
                                          newdata = xTest_SVMLDA),
                           reference = testClasses)
  
  tuneGridResults[i, "testAcc"] <- iconf$overall["Accuracy"]
  tuneGridResults[i, "testKappa"] <- iconf$overall["Kappa"]
  tuneGridResults[i, "pvalue"] <- iconf$overall["AccuracyPValue"]
  
  cat(paste("C: ", as.character(tuneGridSVM[i, "C"]), "; ",
            "sigma: ", as.character(tuneGridSVM[i, "sigma"]), "; ",
            "CV accuracy: ", as.character(round(tuneGridResults[i, "CVacc"], 2)), "; ",
            "test accuracy: ", as.character(round(tuneGridResults[i, "testAcc"], 2)), ".\n",
            sep = ""))
  
}

save(x = tuneGridResults,
     file = "Z:/USUARIOS/JVINDEL/fatima/20230803_PRM_HUMS/rda/tuneGridResults.rda")

bestLDASVM <- tuneGridResults[tuneGridResults$pvalue == min(tuneGridResults$pvalue),]
bestLDASVM <- bestLDASVM[bestLDASVM$CVacc == max(bestLDASVM$CVacc),]
bestLDASVM <- bestLDASVM[bestLDASVM$CVkappa == max(bestLDASVM$CVkappa),]

set.seed(1)
bestLDASVM_svm <- train(as.formula(paste(classColName, " ~ ", ".", sep = "")),
                        data = trainDataFrame_SVMLDA,
                        method = "svmRadial",
                        trControl = trainControl(method = "LOOCV",
                                                 classProbs = TRUE),
                        tuneGrid = expand.grid(C = bestLDASVM[1, "C"],
                                               sigma = bestLDASVM[1, "sigma"]))

testData_SVMLDA <- testData[, SVMLDA_features, drop = FALSE]

xTest_SVMLDA <- predict(object = lda_SVMLDA$finalModel,
                        newdata = testData_SVMLDA)$x

confusionMatrix(data = predict(object = bestLDASVM_svm,
                               newdata = xTest_SVMLDA),
                reference = testClasses)

observedROC_LDASVM <- c()
predictedROC_LDASVM <- c()

testPredictionProbs_LDASVM <- predict(object = bestLDASVM_svm,
                                      newdata = xTest_SVMLDA,
                                      type = "prob")

for(i in 1:nlevels(trClasses)){
  
  iclass <- levels(trClasses)[i]
  iobserved <- as.numeric(as.character(testClasses) == iclass)
  observedROC_LDASVM <- append(x = observedROC_LDASVM, values = iobserved)
  
  ipredicted <- testPredictionProbs_LDASVM[, iclass]
  predictedROC_LDASVM <- append(x = predictedROC_LDASVM, values = ipredicted)
  
}

finalModelROC_LDASVM <- verification::roc.plot(x = observedROC_LDASVM,
                                               pred = predictedROC_LDASVM,
                                               xlab =  "1 - specificity",
                                               ylab = "Sensitivity",
                                               show.thres = FALSE)

# SVMLDA_features
# lda_SVMLDA
# xTrain_SVMLDA
# xTest_SVMLDA
# bestSigma <- 0.1601
# bestC <- 0.31
# bestLDASVM_svm

lda_SVMLDA$finalModel$scaling

a <- abs(lda_SVMLDA$finalModel$scaling[, 1])
a <- a[order(a, decreasing = TRUE)]
a <- (a / sum(a)) * 100

library(randomForest)

set.seed(1)
SVMLDA_RF <- randomForest(x = rbind(xTrain_SVMLDA, xTest_SVMLDA),
                          y = c(trClasses, testClasses),
                          proximity = TRUE,
                          importance = TRUE,
                          ntree = 1000)

randomForest(x = xTrain_SVMLDA,
             y = trClasses,
             proximity = TRUE,
             importance = TRUE,
             ntree = 1000)

randomForest(x = xTest_SVMLDA,
             y = testClasses,
             proximity = TRUE,
             importance = TRUE,
             ntree = 1000)

resultsRFEx_SVMLDA <- RFClassEx(rfObject = SVMLDA_RF,
                                trueClassesOrder =
                                  c("NHOSP", "HOSP", "UCI", "EXI"),
                                classesColors =
                                  covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
                                cexMain = 2,
                                cexAxisNumbers = 1.3,
                                cexPlotAxis = 1.4,
                                tiltAngle = 40,
                                cexLabels = 2,
                                cexLeg = 2,
                                legInset = c(0.02, 0.02),
                                adjText = c(1, 1),
                                yLabelsPosition = -0.01)

set.seed(1)
SVMLDA_features_RF <- randomForest(x = scaledCovidData[, SVMLDA_features],
                                   y = covidClasses,
                                   proximity = TRUE,
                                   importance = TRUE,
                                   ntree = 1000)

resultsRFEx_SVMLDA <- RFClassEx(rfObject = SVMLDA_features_RF,
                                trueClassesOrder =
                                  c("NHOSP", "HOSP", "UCI", "EXI"),
                                classesColors =
                                  covidColors[c("NHOSP", "HOSP", "UCI", "EXI")],
                                cexMain = 2,
                                cexAxisNumbers = 1.2,
                                cexPlotAxis = 1.2,
                                tiltAngle = 0,
                                cexLabels = 1.5,
                                cexLeg = 1.9,
                                legInset = c(0.02, 0.02),
                                adjText = c(1, 1),
                                yLabelsPosition = -0.0025)

library(plotrix)

chosenSpace <- 0.2
chosenBarWidth <- 1
xFactor1 <- 26.6
yFactor1 <- 0.1750001
xFactor2 <- 40
yFactor2 <- 0.1750001
myColorGradient <- c("royalblue3", "white", "orange")
colorGradientLength <- 10

constantPalette <- colorRampPalette(colors = myColorGradient)(colorGradientLength)

for(i in 1:ncol(lda_SVMLDA$finalModel$scaling)){
  
  ild <- lda_SVMLDA$finalModel$scaling[, i]
  ild <- ild[order(abs(ild), decreasing = TRUE)]
  
  ildSign <- ild
  ildSign[ildSign < 0] <- -1
  ildSign[ildSign >= 0] <- 1
  
  ild <- (abs(ild) / sum(abs(ild))) * ildSign
  
  variableIntervals <- unique(c(seq(from = min(ild), to = 0, length.out = 6),
                                seq(from = 0, to = max(ild), length.out = 6)))
  
  consecutiveMeans <- c()
  
  for(j in 1:(length(variableIntervals) - 1)){
    
    jmean <- (variableIntervals[j] + variableIntervals[j + 1]) / 2
    consecutiveMeans <- append(x = consecutiveMeans, values = jmean)
    
  }
  
  ichosenColors <- unlist(sapply(X = ild,
                                 FUN = function(x)
                                   constantPalette[which.min(abs(consecutiveMeans - x))]))
  
  par(mar = c(5.1, # the greater the values, the less space from below
              4.1, # the greater the values, the less space at left hand
              4.1, # the greater the values, the less space from above
              7.1 # the greater the values, the less space at right hand
  )
  ) # Default values
  
  barplot(height = ild,
          xaxt = "n",
          yaxt = "n",
          ylim = range(pretty(ild)),
          col = ichosenColors,
          space = chosenSpace,
          width = chosenBarWidth)
  
  axis(side = 2, at = pretty(ild), labels = rep("", length(pretty(ild))))
  
  xLabels <- rep(chosenSpace + (chosenBarWidth/2), length(ild)) +
    seq(from = 0,
        to = (chosenSpace + chosenBarWidth) * (length(ild) - 1),
        by = chosenSpace + chosenBarWidth)
  
  xLow <- xLabels[1]
  xHigh <- xLabels[length(xLabels)]
  yLow <- range(pretty(ild))[1]
  yHigh <- range(pretty(ild))[2]
  
  yInterval <- pretty(ild)[2] - pretty(ild)[1]
  
  xa1 <- chosenSpace + chosenBarWidth
  xa2 <- xHigh + (chosenBarWidth/2)
  
  prettyXPos1 <- xLow - (xFactor1 * (xa1 / xa2))
  prettyXPos2 <- xLow - (xFactor2 * (xa1 / xa2))
  
  prettyYPos1 <- yLow - (yFactor1 * (yInterval / (yHigh + abs(yLow))))
  prettyYPos2 <- yHigh + (yFactor2 * (yInterval / (yHigh + abs(yLow))))
  
  text(x = prettyXPos1,
       y = pretty(ild),
       labels = pretty(ild),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.1)
  text(x = xLabels,
       y = prettyYPos1,
       labels = names(ild),
       adj = c(0.7, 0),
       xpd = NA,
       srt = 35,
       cex = 1.1)
  
  text(x = prettyXPos2,
       y = mean(range(pretty(ild))),
       labels = "Variance proportion contributed by each feature",
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 1.4)
  text(x = mean(c(xLabels[length(xLabels)], xLabels[1])),
       y = prettyYPos2,
       labels = paste("Contributed variance proportion by feature (LD",
                      as.character(i),
                      ")",
                      sep = ""),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 2)
  
  xLegendl <- xHigh + (2 * (chosenBarWidth/2)) + (2 * chosenSpace)
  xLegendr <- xLegendl + (1.5 * chosenBarWidth)
  xLegendText <- (7.5 * ((2 * chosenBarWidth) / xLegendr)) + xLegendr
  
  color.legend(xl = xLegendl,
               yb = yLow,
               xr = xLegendr,
               yt = yHigh,
               legend = "",
               rect.col = constantPalette,
               gradient = "y")
  
  text(x = xLegendText,
       y = round(seq(from = yLow,
                     to = yHigh,
                     length.out = length(constantPalette) + 1),
                 2),
       labels = round(variableIntervals, 2),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.4)
  
}

allLDCombn <- combn(x = colnames(lda_SVMLDA$finalModel$scaling), m = 2)
isText <- FALSE

for(i in 1:ncol(allLDCombn)){
  
  icol <- allLDCombn[, i]
  firstCol <- lda_SVMLDA$finalModel$scaling[, icol[1]]
  
  firstCol <- firstCol[order(abs(firstCol), decreasing = TRUE)]
  
  firstSign <- firstCol
  firstSign[firstSign < 0] <- -1
  firstSign[firstSign >= 0] <- 1
  
  firstCol <- (abs(firstCol) / sum(abs(firstCol))) * firstSign
  
  secondCol <- lda_SVMLDA$finalModel$scaling[, icol[2]]
  
  secondCol <- secondCol[order(abs(secondCol), decreasing = TRUE)]
  
  secondSign <- secondCol
  secondSign[secondSign < 0] <- -1
  secondSign[secondSign >= 0] <- 1
  
  secondCol <- (abs(secondCol) / sum(abs(secondCol))) * firstSign
  
  if(max(abs(c(firstCol, secondCol))) > 0.27){
    
    squareLimits <- c(-0.5, 0.5)
    
  }else{
    
    squareLimits <- c(-0.25, 0.25)
    
  }
  
  
  if(min(squareLimits) == -0.25){
    
    prettyAxis <- c(-0.27, -0.2, -0.1, 0, 0.1, 0.2, 0.27)
    xPrettyLab <- -0.315
    yPrettyLab <- -0.31
    xPrettyAxis <- -0.285
    yPrettyAxis <- -0.288
    yPrettyTitle <- 0.3
    textDisplacement <- 0.025
    
  }
  if(min(squareLimits) == -0.5){
    
    prettyAxis <- c(-0.54, -0.4, -0.2, 0, 0.2, 0.4, 0.54)
    xPrettyLab <- -0.63
    yPrettyLab <- -0.62
    xPrettyAxis <- -0.57
    yPrettyAxis <- -0.576
    yPrettyTitle <- 0.6
    textDisplacement <- 0.05
    
  }
  
  icolData <- cbind(firstCol, secondCol[names(firstCol)])
  colnames(icolData) <- icol
  
  firstDimension <- as.numeric(gsub(pattern = "LD", replacement = "", x = icol[1]))
  secondDimension <- as.numeric(gsub(pattern = "LD", replacement = "", x = icol[2]))
  
  ixlab <- (((lda_SVMLDA$finalModel$svd^2) / sum(lda_SVMLDA$finalModel$svd^2)) * 100)
  iylab <- ixlab[secondDimension]
  iylab <- paste("Eigenvector feature multiplier absolute value proportion ",
                 "(LD",
                 as.character(secondDimension),
                 ") (",
                 as.character(round(iylab, 2)),
                 " %)",
                 sep = "")
  ixlab <- ixlab[firstDimension]
  ixlab <- paste("Eigenvector feature multiplier absolute value proportion ",
                 "(LD",
                 as.character(firstDimension),
                 ") (",
                 as.character(round(ixlab, 2)),
                 " %)",
                 sep = "")
  
  plot(x = 0,
       y = 0,
       xlim = squareLimits,
       ylim = squareLimits,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       col = "white")
  
  abline(v = pretty(squareLimits),
         h = pretty(squareLimits))
  abline(v = 0,
         h = 0,
         lwd = 3)
  
  # xAxis
  text(x = prettyAxis,
       y = xPrettyAxis,
       labels = prettyAxis,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.4)
  
  # yAxis
  text(x = yPrettyAxis,
       y = prettyAxis,
       labels = prettyAxis,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.4)
  
  # xlab
  text(x = 0,
       y = xPrettyLab,
       labels = ixlab,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.75)
  
  # ylab
  text(x = yPrettyLab,
       y = 0,
       labels = iylab,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 1.75)
  
  # main
  text(x = 0,
       y = yPrettyTitle,
       labels = paste("LDA biplot (", icol[1], " vs. ", icol[2], ")", sep = ""),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 3)
  
  for(j in 1:nrow(icolData)){
    
    arrows(x0 = 0,
           y0 = 0,
           x1 = icolData[j, 1],
           y1 = icolData[j, 2],
           length = 0.3,
           angle = 25,
           code = 2, # arrow sense
           col = "black",
           lty = 1,
           lwd = 3)
    
    if(icolData[j, 1] > 0){
      xSign <- 1
    }else{
      xSign <- -1
    }
    
    if(icolData[j, 2] > 0){
      ySign <- 1
    }else{
      ySign <- -1
    }
    
    
    if(isText){
      
      text(x = icolData[j, 1] + (textDisplacement * xSign),
           y = icolData[j, 2] + (textDisplacement * ySign),
           labels = rownames(icolData)[j],
           adj = c(0.5, 0.5),
           xpd = NA,
           srt = 0,
           cex = 1.5)
      
    }
    
  }
  
}

ldArr <- lda_SVMLDA$finalModel$scaling
ldArr[!is.na(ldArr)] <- NA

for(i in 1:ncol(lda_SVMLDA$finalModel$scaling)){
  
  ild <- lda_SVMLDA$finalModel$scaling[, i]
  ild <- ild[order(abs(ild), decreasing = TRUE)]
  
  ildSign <- ild
  ildSign[ildSign < 0] <- -1
  ildSign[ildSign >= 0] <- 1
  
  ild <- (abs(ild) / sum(abs(ild))) * ildSign
  
  ldArr[, i] <- ild[rownames(ldArr)]
  
}

lighteningFactor <- 0.6
lightenedCovidColors <- lighten(col = covidColors, amount = lighteningFactor)
names(lightenedCovidColors) <- names(covidColors)

xAll_SVMLDA <- rbind(xTrain_SVMLDA, xTest_SVMLDA)
allClasses_SVMLDA <- c(trClasses, testClasses)

allLDAVars <- (lda_SVMLDA$finalModel$svd^2 / sum(lda_SVMLDA$finalModel$svd^2)) * 100

xTrainDF <- cbind.data.frame(xTrain_SVMLDA, trClasses)
colnames(xTrainDF)[ncol(xTrainDF)] <- classColName

yPrettyFactor1 <- 33/400 # Numbers ('x' axis)
yPrettyFactor2 <- 6/40 # Label ('x' axis)

xPrettyFactor1 <- 30/400 # Numbers ('y' axis)
xPrettyFactor2 <- 475/4000 # Label ('y' axis)

yPrettyFactor3 <- 145/1500 # Title

for(i in 1:ncol(allLDCombn)){
  
  firstLD <- as.numeric(gsub(pattern = "LD", replacement = "", x = allLDCombn[1, i]))
  secondLD <- as.numeric(gsub(pattern = "LD", replacement = "", x = allLDCombn[2, i]))
  
  firstLDLabel <- paste("LD",
                        as.character(firstLD),
                        " (",
                        as.character(round(allLDAVars[firstLD], 2)),
                        " %)",
                        sep = "")
  secondLDLabel <- paste("LD",
                         as.character(secondLD),
                         " (",
                         as.character(round(allLDAVars[secondLD], 2)),
                         " %)",
                         sep = "")
  
  xLowest <- range(pretty(xAll_SVMLDA[, firstLD]))[1]
  xHighest <- range(pretty(xAll_SVMLDA[, firstLD]))[2]
  xInterval <- xLowest - xHighest
  
  yLowest <- range(pretty(xAll_SVMLDA[, secondLD]))[1]
  yHighest <- range(pretty(xAll_SVMLDA[, secondLD]))[2]
  yInterval <- yLowest - yHighest
  
  xAxisNumbers <- yLowest + (yInterval * yPrettyFactor1)
  xAxisLabel <- yLowest + (yInterval * yPrettyFactor2)
  
  yAxisNumbers <- xLowest + (xInterval * xPrettyFactor1)
  yAxisLabel <- xLowest + (xInterval * xPrettyFactor2)
  
  yMain <- yHighest + (-(yInterval) * yPrettyFactor3)
  
  plot(x = xAll_SVMLDA[, firstLD],
       y = xAll_SVMLDA[, secondLD],
       xlim = range(pretty(xAll_SVMLDA[, firstLD])),
       ylim = range(pretty(xAll_SVMLDA[, secondLD])),
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       pch = 21,
       col = "black",
       bg = covidColors[as.character(allClasses_SVMLDA)],
       cex = 2)
  
  axis(side = 1,
       at = pretty(xAll_SVMLDA[, firstLD]),
       labels = rep("", length(pretty(xAll_SVMLDA[, firstLD]))))
  axis(side = 2,
       at = pretty(xAll_SVMLDA[, secondLD]),
       labels = rep("", length(pretty(xAll_SVMLDA[, secondLD]))))
  
  # x
  text(x = pretty(xAll_SVMLDA[, firstLD]),
       y = xAxisNumbers,
       labels = pretty(xAll_SVMLDA[, firstLD]),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.7)
  text(x = mean(c(xLowest, xHighest)),
       y = xAxisLabel,
       labels = firstLDLabel,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 2)
  
  # y
  text(x = yAxisNumbers,
       y = pretty(xAll_SVMLDA[, secondLD]),
       labels = pretty(xAll_SVMLDA[, secondLD]),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 1.7)
  text(x = yAxisLabel,
       y = mean(c(yLowest, yHighest)),
       labels = secondLDLabel,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 2)
  
  # title
  text(x = mean(c(xLowest, xHighest)),
       y = yMain,
       labels = paste("LDA (LD",
                      as.character(firstLD),
                      " vs. ",
                      "LD",
                      as.character(secondLD),
                      ")",
                      sep = ""),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 2.5)
  
  ##############################################################################
  
  tempFormula <- as.formula(paste(classColName,
                                  " ~ ",
                                  paste(allLDCombn[, i],
                                        collapse = " + "),
                                  sep = ""))
  
  set.seed(1)
  tempSVM <- train(tempFormula,
                   data = xTrainDF,
                   method = "svmRadial",
                   trControl = trainControl(method = "LOOCV",
                                            classProbs = TRUE),
                   tuneGrid = data.frame(C = 0.31, sigma = 0.1601))
  
  xGrid1 <- xLowest - ((pretty(xAll_SVMLDA[, firstLD])[2] - pretty(xAll_SVMLDA[, firstLD])[1])/2)
  xGrid2 <- xHighest + ((pretty(xAll_SVMLDA[, firstLD])[2] - pretty(xAll_SVMLDA[, firstLD])[1])/2)
  yGrid1 <- yLowest - ((pretty(xAll_SVMLDA[, secondLD])[2] - pretty(xAll_SVMLDA[, secondLD])[1])/2)
  yGrid2 <- yHighest + ((pretty(xAll_SVMLDA[, secondLD])[2] - pretty(xAll_SVMLDA[, secondLD])[1])/2)
  xGrid <- seq(from = xGrid1, to = xGrid2, length.out = 707)
  yGrid <- seq(from = yGrid1, to = yGrid2, length.out = 707)
  Grid <- expand.grid(x = xGrid, y = yGrid)
  colnames(Grid) <- allLDCombn[, i]
  
  GridClasses <- predict(object = tempSVM, newdata = Grid)
  
  plot(x = Grid[, 1],
       y = Grid[, 2],
       xlim = range(pretty(xAll_SVMLDA[, firstLD])),
       ylim = range(pretty(xAll_SVMLDA[, secondLD])),
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       pch = 15,
       col = lightenedCovidColors[as.character(GridClasses)],
       cex = 0.2)
  
  points(x = xAll_SVMLDA[, firstLD],
         y = xAll_SVMLDA[, secondLD],
         pch = 21,
         col = "black",
         bg = covidColors[as.character(allClasses_SVMLDA)],
         cex = 2)
  
  axis(side = 1,
       at = pretty(xAll_SVMLDA[, firstLD]),
       labels = rep("", length(pretty(xAll_SVMLDA[, firstLD]))))
  axis(side = 2,
       at = pretty(xAll_SVMLDA[, secondLD]),
       labels = rep("", length(pretty(xAll_SVMLDA[, secondLD]))))
  
  # x
  text(x = pretty(xAll_SVMLDA[, firstLD]),
       y = xAxisNumbers,
       labels = pretty(xAll_SVMLDA[, firstLD]),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 1.7)
  text(x = mean(c(xLowest, xHighest)),
       y = xAxisLabel,
       labels = firstLDLabel,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 2)
  
  # y
  text(x = yAxisNumbers,
       y = pretty(xAll_SVMLDA[, secondLD]),
       labels = pretty(xAll_SVMLDA[, secondLD]),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 1.7)
  text(x = yAxisLabel,
       y = mean(c(yLowest, yHighest)),
       labels = secondLDLabel,
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 90,
       cex = 2)
  
  # title
  text(x = mean(c(xLowest, xHighest)),
       y = yMain,
       labels = paste("LDA (LD",
                      as.character(firstLD),
                      " vs. ",
                      "LD",
                      as.character(secondLD),
                      ")",
                      sep = ""),
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = 2.5)
  
}

set.seed(1)
l <- train(patient ~ .,
           data = trDataFrame[, c("A2MG", "RET4", "CYTC", "CO5", "CLUS",
                                  "FETUA", "FIBG", "FIBB", "GNPTG", "FINC",
                                  "K1C16", "ALBU", "K1C9", "IBP2", "patient")],
           method = "lda",
           trControl = trainControl(method = "LOOCV",
                                    classProbs = TRUE))

xtrainLDAdata <- predict(object = l$finalModel,
                         newdata = trData[, c("A2MG", "RET4", "CYTC", "CO5", "CLUS",
                                              "FETUA", "FIBG", "FIBB", "GNPTG", "FINC",
                                              "K1C16", "ALBU", "K1C9", "IBP2")])$x
xtestLDAdata <- predict(object = l$finalModel,
                        newdata = testData[, c("A2MG", "RET4", "CYTC", "CO5", "CLUS",
                                               "FETUA", "FIBG", "FIBB", "GNPTG", "FINC",
                                               "K1C16", "ALBU", "K1C9", "IBP2")])$x

xAll_SVMLDA <- rbind(xtrainLDAdata, xtestLDAdata)

plot3D(data3D = xAll_SVMLDA,
       classVector = as.character(c(trClasses, testClasses)),
       classColors = covidColors[order(as.character(unique(trClasses)))],
       main3D = "",
       dimension1Name = "LD1 (66.26 %)",
       dimension2Name = "LD2  (26.35 %)",
       dimension3Name = "LD3  (7.39 %)",
       cex3D = 1.5,
       cexLegend = 1.5,
       xAxisLimits = c(-5, 5),
       yAxisLimits = c(-4, 5),
       zAxisLimits = c(-3.5, 3),
       par3D = 2.5)

PCAGif()

################################################################################

trData

# ytrClasses <- matrix(data = NA,
#                      nrow = length(trClasses),
#                      ncol = nlevels(trClasses),
#                      dimnames = list(rownames(trData),
#                                      as.character(levels(trClasses))))
# for(i in 1:nrow(ytrClasses)){
#   
#   ytrClasses[i,] <- as.numeric(colnames(ytrClasses) %in% as.character(trClasses)[i])
#   
# }

# or

ytrClasses <- t(sapply(X = as.character(trClasses),
                       FUN = function(x) as.numeric(as.character(levels(trClasses)) %in% x)))
rownames(ytrClasses) <- rownames(trData)
colnames(ytrClasses) <- as.character(levels(trClasses))

set.seed(1)
sdaModel <- sda(trData,
                ytrClasses,
                lambda = 1e-6,
                stop = -1,
                maxIte = 25,
                trace = TRUE)

sdaModel$call
sdaModel$beta
sdaModel$theta
sdaModel$varNames
sdaModel$varIndex
sdaModel$origP
sdaModel$rss

sdaModel$fit$prior
sdaModel$fit$counts
sdaModel$fit$means
sdaModel$fit$scaling
sdaModel$fit$lev
sdaModel$fit$svd
sdaModel$fit$N
sdaModel$fit$call

sdaModel$classes
sdaModel$lambda
sdaModel$stop

sum(as.character(predict(object = sdaModel$fit,
                         newdata = trData[, sdaModel$varNames])$class) ==
      as.character(trClasses)) /
  length(trClasses)

set.seed(1)
smdaModel <- sparseLDA::smda(trData,
                             ytrClasses,
                             lambda = 1e-6,
                             stop = -5,
                             maxIte = 25,
                             tol = 1e-2)

smdaModel$call
smdaModel$beta
smdaModel$theta
smdaModel$Z
smdaModel$Rj
smdaModel$K
smdaModel$mu
smdaModel$Sigma
smdaModel$Sigma_inv
smdaModel$dpi
smdaModel$varNames
smdaModel$varIndex
smdaModel$origP
smdaModel$rss

smdaModel$fit$prior
smdaModel$fit$counts
smdaModel$fit$means
smdaModel$fit$scaling
smdaModel$fit$lev
smdaModel$fit$svd
smdaModel$fit$N
smdaModel$fit$call

smdaModel$classes
smdaModel$subClasses
smdaModel$lambda
smdaModel$stop

sum(as.character(predict(object = smdaModel$fit,
                         newdata = trData[, smdaModel$varNames])$class) ==
      as.character(trClasses)) /
  length(trClasses)

################################################################################

library(sparseLDA)

## load data
data(penicilliumYES)
X <- penicilliumYES$X
Y <- penicilliumYES$Y
colnames(Y) <- c("P. Melanoconidium",
                 "P. Polonicum",
                 "P. Venetum")
## test samples
Iout<-c(3,6,9,12)
Iout<-c(Iout,Iout+12,Iout+24)
## training data
Xtr<-X[-Iout,]
k<-3
n<-dim(Xtr)[1]
## Normalize data
Xc<-normalize(Xtr)
Xn<-Xc$Xc
p<-dim(Xn)[2]
## Perform SDA with one non-zero loading for each discriminative
## direction with Y as matrix input
out <- sda(Xn, Y,
           lambda = 1e-6,
           stop = -1,
           maxIte = 25,
           trace = TRUE)
## predict training samples
train <- predict(out, Xn)
## testing
Xtst<-X[Iout,]
Xtst<-normalizetest(Xtst,Xc)

test <- predict(out, Xtst)
print(test$class)
## Factor Y as input
Yvec <- factor(rep(colnames(Y), each = 8))
out2 <- sda(Xn, Yvec,
            lambda = 1e-6,
            stop = -1,
            maxIte = 25,
            trace = TRUE)

################################################################################

library(glmnet)

#perform k-fold cross-validation to find optimal lambda value
set.seed(1)
lassoModel <- cv.glmnet(x = trData,
                        y = as.numeric(trClasses) - 1,
                        alpha = 1,
                        nfolds = nrow(trData),
                        grouped = FALSE)

#find optimal lambda value that minimizes test MSE
bestLambda <- lassoModel$lambda.min

#produce plot of test MSE by lambda value
plot(lassoModel)

#find coefficients of best model
optimLasso <- glmnet(trData,
                     as.numeric(trClasses) - 1,
                     alpha = 1,
                     lambda = bestLambda)
lassoProts <- rownames(coef(optimLasso))[-1][which((abs(coef(optimLasso)[, 1]) > 0)[-1])]

set.seed(1)
trGLMNET <- train(as.formula(paste(classColName,
                                   " ~ ",
                                   paste(lassoProts, collapse = " + "))),
                  data = trDataFrame,
                  method = "glmnet",
                  trControl = trainControl(method = "LOOCV", classProbs = TRUE))
trGLMNET
confusionMatrix(data = predict(object = trGLMNET, newdata = testData[, lassoProts]),
                reference = testClasses)

set.seed(1)
trLDA <- train(as.formula(paste(classColName,
                                " ~ ",
                                paste(lassoProts, collapse = " + "))),
               data = trDataFrame,
               method = "lda",
               trControl = trainControl(method = "LOOCV", classProbs = TRUE))
trLDA
confusionMatrix(data = predict(object = trLDA, newdata = testData[, lassoProts]),
                reference = testClasses)

set.seed(1)
trSVM <- train(as.formula(paste(classColName,
                                " ~ ",
                                ".")),
               data = cbind.data.frame(predict(object = trLDA$finalModel,
                                               newdata = trData[, lassoProts])$x,
                                       "patient" = trClasses),
               method = "svmRadial",
               trControl = trainControl(method = "LOOCV", classProbs = TRUE))
trSVM
confusionMatrix(data = predict(object = trSVM,
                               newdata = predict(object = trLDA$finalModel,
                                                 newdata = testData[, lassoProts])$x),
                reference = testClasses)

################################################################################

orderIWant <- c("NHOSP", "HOSP", "UCI", "EXI")

wantColor <- TRUE

dev.new()

par(mfrow = c(2, 4))

dataToBoxs <- scaledCovidData[, c("A2MG", "RET4", "CYTC", "CO5", "CLUS",
                                  "FETUA", "FIBG", "FIBB", "GNPTG", "FINC",
                                  "K1C16", "ALBU", "K1C9", "IBP2")[1:8]]

boxsPointCex <- 2
yLabelsFactor <- 3/45
yAxisCex <- 2.5
labelsTilt <- 30
labelsCex <- 2.3
titleFactor <- 1/9
titleCex <- 4

for(i in 1:ncol(dataToBoxs)){
  
  idat <- cbind.data.frame(dataToBoxs[, i],
                           as.factor(covidClasses))
  rownames(idat) <- rownames(dataToBoxs)
  colnames(idat) <- c(colnames(dataToBoxs)[i], "class")
  idat$class <- factor(idat$class, levels = orderIWant)
  boxplotFormula <- paste(colnames(idat)[2], " ~ ", colnames(idat)[1], sep = "")
  boxplotFormula <- as.formula(boxplotFormula)
  
  if(wantColor){
    colorsToBoxs <- covidColors[orderIWant]
  }else{
    colorsToBoxs <- "grey"
  }
  
  boxplot(idat[,1] ~ idat[,2],
          data = idat,
          xaxt = "n",
          yaxt = "n",
          ylim = range(pretty(idat[,1])),
          xlab = "",
          main = "",
          ylab = "",
          # ylab = paste(colnames(idat)[1],
          #              " normalized abundance levels",
          #              sep = ""),
          pch = 21,
          col = colorsToBoxs,
          bg = colorsToBoxs,
          outcex = boxsPointCex)
  axis(side = 1,
       at = 1:nlevels(idat[, 2]),
       labels = rep("", nlevels(idat[, 2])))
  axis(side = 2,
       las = 1,
       at = pretty(idat[,1]),
       labels = pretty(idat[,1]),
       cex.axis = yAxisCex)
  
  highest <- range(pretty(idat[,1]))[2]
  lowest <- range(pretty(idat[,1]))[1]
  difference <- highest - lowest
  labelsYPosition <- lowest - (yLabelsFactor * difference)
  titleYPosition <- highest + (titleFactor * difference)
  
  text(x = 1:nlevels(idat[,2]),
       y = labelsYPosition,
       labels = orderIWant,
       adj = c(1, 1),
       xpd = NA,
       srt = labelsTilt,
       cex = labelsCex)
  
  text(x = mean(1:nlevels(idat[,2])),
       y = titleYPosition,
       labels = uniprotDescription[colnames(idat)[1], "alternate.gene"],
       adj = c(0.5, 0.5),
       xpd = NA,
       srt = 0,
       cex = titleCex)
  
}

l$finalModel$scaling

write.xlsx2(x = l$finalModel$scaling,
            file = "LDAeigenvectors.xlsx",
            sheetName = "eigen",
            row.names = TRUE,
            col.names = TRUE,
            append = FALSE)

for(i in 1:ncol(l$finalModel$scaling)){
  
  icol <- l$finalModel$scaling[, i]
  iorder <- order(abs(icol), decreasing = TRUE)
  icol <- icol[iorder]
  
  isense <- icol >= 0
  isense[isense == TRUE] <- "Positive"
  isense[isense != "Positive"] <- "Negative"
  
  idataFrame <- data.frame("original.values" = icol,
                           "sense" = isense,
                           "absolute" = abs(icol),
                           "weight.proportion" = abs(icol) / sum(abs(icol)))
  
  write.xlsx2(x = idataFrame,
              file = "LDAeigenvectors.xlsx",
              sheetName = colnames(l$finalModel$scaling)[i],
              row.names = TRUE,
              col.names = TRUE,
              append = TRUE)
  
}

################################################################################

combinedModelProteins <- c("A2MG", "RET4", "CYTC", "CO5", "CLUS",
                           "FETUA", "FIBG", "FIBB", "GNPTG", "FINC",
                           "K1C16", "ALBU", "K1C9", "IBP2")
set.seed(1)
finalSVMmodel <- train(patient ~ .,
                       data = cbind.data.frame(xtrainLDAdata,
                                               "patient" = trClasses),
                       method = "svmRadial",
                       trControl = trainControl(method = "LOOCV",
                                                classProbs = TRUE),
                       tuneGrid = data.frame("C" = 0.31,
                                             "sigma" = 0.1601))



confusionMatrix(data = finalSVMmodel$pred$pred, reference = trClasses)

CVpredictionProbs_LDASVM <- finalSVMmodel$pred[, levels(trClasses)]

confusionMatrix(data = predict(object = finalSVMmodel,
                               newdata = xtestLDAdata),
                reference = testClasses)

observedROC_LDASVM <- c()
predictedROC_LDASVM <- c()

testPredictionProbs_LDASVM <- predict(object = finalSVMmodel,
                                      newdata = xtestLDAdata,
                                      type = "prob")

for(i in 1:nlevels(trClasses)){
  
  iclass <- levels(trClasses)[i]
  iobserved <- as.numeric(as.character(testClasses) == iclass)
  observedROC_LDASVM <- append(x = observedROC_LDASVM, values = iobserved)
  
  ipredicted <- testPredictionProbs_LDASVM[, iclass]
  predictedROC_LDASVM <- append(x = predictedROC_LDASVM, values = ipredicted)
  
}

finalModelROC_LDASVM <- verification::roc.plot(x = observedROC_LDASVM,
                                               pred = predictedROC_LDASVM,
                                               xlab =  "1 - specificity",
                                               ylab = "Sensitivity",
                                               show.thres = FALSE)

xROC <- rev(matrix(finalModelROC_LDASVM[["plot.data"]],
                   ncol = ncol(finalModelROC_LDASVM[["plot.data"]]))[, 3])
yROC <- rev(matrix(finalModelROC_LDASVM[["plot.data"]],
                   ncol = ncol(finalModelROC_LDASVM[["plot.data"]]))[, 2])
newXROC <- c()
newYROC <- c()

for(j in 1:length(xROC)){
  
  jx1 <- xROC[j]
  jy1 <- yROC[j]
  
  if(j != 1){
    
    jx0 <- xROC[j - 1]
    jy0 <- yROC[j - 1]
    jxBool <- jx0 == jx1
    jyBool <- jy0 == jy1
    
    if(as.numeric(jxBool) + as.numeric(jyBool) != 1){
      
      newXROC <- append(x = newXROC, values = mean(c(xROC[j - 1], xROC[j])))
      newYROC <- append(x = newYROC, values = yROC[j - 1])
      
      newXROC <- append(x = newXROC, values = mean(c(xROC[j - 1], xROC[j])))
      newYROC <- append(x = newYROC, values = yROC[j])
      
    }
    
  }
  
  newXROC <- append(x = newXROC, values = xROC[j])
  newYROC <- append(x = newYROC, values = yROC[j])
  
}

plot(x = newXROC,
     y = newYROC,
     xlim = range(pretty(newXROC)),
     ylim = range(pretty(newYROC)),
     yaxt = "n",
     type = "l",
     col = "black",
     lty = 1,
     lwd = 3,
     main = "LDA and SVM combined model ROC curve",
     xlab = "1 - specificity",
     ylab = "Sensitivity",
     cex.axis = 1.5,
     cex.lab = 1.4,
     cex.main = 2)
axis(side = 2,
     las = 2,
     at = pretty(newYROC),
     labels = pretty(newYROC),
     cex.axis = 1.5)

abline(v = pretty(newXROC), h = pretty(newYROC), col = "grey", lty = 1, lwd = 1)
segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col = "darkgrey", lty = 1, lwd = 2)

points(x = newXROC,
       y = newYROC,
       type = "l",
       col = "black",
       lty = 1,
       lwd = 3)

ROClegend <- c(paste("AUC: ",
                     as.character(round(finalModelROC_LDASVM$roc.vol$Area, 4)),
                     sep = ""),
               paste("ROC p-value: ",
                     format(x = round(finalModelROC_LDASVM$roc.vol$p.value, 4),
                            scientific = FALSE),
                     sep = ""))

legend(x = "bottomright",
       inset = c(0.025, 0.04),
       legend = ROClegend,
       pch = 22,
       col = NA,
       pt.bg = NA,
       cex = 2,
       pt.cex = 2)
