library(Mfuzz)

# Load data

data("yeast")
View(yeast)
yeastArr <- yeast@assayData$exprs
View(yeastArr)

# Remove genes with missing values (threshold: 25 %)

yeast.r <- filter.NA(yeast, thres = 0.25)
View(yeast.r)

# Impute missing values

yeast.f <- fill.NA(yeast.r, mode = "mean") # available modes: 'mean', 'knn', 'wknn'.
# "wknn" mode does not completely impute NA values. Could then chained modes
# be used ("wknn" + "mean")?
View(yeast.f@assayData$exprs)

# Filter unchanged genes over time

tmp <- filter.std(yeast.f, min.std = 0)
# generates a plot.

# Standardize expression values

yeast.s <- standardise(yeast.f)

# Soft-cluster

cl <- mfuzz(yeast.s, c = 16, m = 1.25)
# 'c' represents the number of clusters that are wanted to be identified.
# 'm' is defined as the "fuzzification" parameter.
# An object that indicates the number of clusters and what feature belongs to
# what cluster.
# cl$cluster

View(yeast.s@assayData$exprs)
mfuzz.plot(eset = yeast.s,
           cl = cl,
           mfrow = c(4, 4),
           time.labels = seq(0, 160, 10))

m1 <- mestimate(yeast.s) # !!!
m1
c1 <- cselection(eset = yeast.s, m = m1, crange = 2:100)
m2 <- partcoef(eset = yeast.s)
View(m2)
m2[[1]]
m2[[2]]
m2[[3]]
longRange <- 2:100
d1 <- Dmin(eset = yeast.s, m = m1, crange = longRange) # !!!
plot(x = longRange,
     y = d1,
     pch = 16,
     cex = 0.8,
     col = "black",
     xlab = "Number of clusters",
     ylab = "Minimum distance 'D min' between cluster centroid",
     main = "Minimum distance between cluster centroid vs. Number of clusters",
     xlim = range(pretty(longRange)),
     ylim = range(pretty(d1)))
points(x = longRange,
       y = d1,
       col = "black",
       type = "l",
       cex = 1)

cl2 <- mfuzz(yeast.s, c = 25, m = m1)
mfuzz.plot(eset = yeast.s,
           cl = cl2,
           mfrow = c(5, 4),
           time.labels = seq(0, 160, 10))

# Cluster stability measure: if 'm' becomes greater, then the most stable
# clusters will not modify their content.

O <- overlap(cl)
# 'O' shows the overlapping between the clusters similarity of gene expression
# over time in a square matrix which number of rows and columns equals
# the number of clusters. The bigger the value of a comparison between two
# clusters named with their number, the more similar they are.

Ptmp <- overlap.plot(cl, over = O, thres = 0.05)

O2 <- overlap(cl2)
Ptmp2 <- overlap.plot(cl2, over = O2, thres = 0.05)

################################################################################

# Transform my data into the format I'm supposed to use.

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
aggFrame <- cbind.data.frame(t(normSARS), "classes" = as.factor(groups))
aggRes <- aggregate(aggFrame[,1:(ncol(aggFrame)-1)], list(aggFrame$classes), mean)
rownames(aggRes) <- aggRes[,1]
aggRes <- aggRes[,-1]
aggRes <- t(aggRes)
saveRows <- rownames(aggRes)
aggRes <- cbind(aggRes, "SARS.CoV2.t16" = aggRes[,2])
aggRes <- aggRes[,-c(2,6)]
# uniqueGroups <- unique(groups) # This is way better :/
# aggRes <- aggRes[,uniqueGroups]
save(x = aggRes,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/aggRes.rda")

# Convert mean columns matrix into a 'ExpressionSet' class object:

aggSet <- ExpressionSet(assayData = aggRes)
View(aggSet)
save(x = aggSet,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/aggSet.rda")

aggTmp <- filter.std(aggSet, min.std = 0)

originalProteins <- rownames(normSARS)
filteredProteins <- rownames(aggTmp@assayData$exprs)
filteredOut <- setdiff(x = originalProteins, y = filteredProteins)

standSARS <- standardise(aggTmp)
save(x = standSARS,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/standSARS.rda")

mSARS <- mestimate(standSARS)
save(x = mSARS,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mSARS.rda")
longRange <- 2:100
dSARS <- Dmin(eset = standSARS, m = mSARS, crange = longRange)
save(x = dSARS,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/dSARS.rda")
plot(x = longRange,
     y = dSARS,
     pch = 16,
     cex = 0.8,
     col = "black",
     xlab = "Number of clusters",
     ylab = "Minimum distance 'D min' between cluster centroid",
     main = "Minimum distance between cluster centroid vs. Number of clusters",
     xlim = range(pretty(longRange)),
     ylim = range(pretty(dSARS)))
points(x = longRange,
       y = dSARS,
       col = "black",
       type = "l",
       cex = 1)

cSARS <- 9
timeSARS <- c(0, 3, 6, 9, 16)
clSARS <- mfuzz(standSARS, c = cSARS, m = mSARS)
save(x = clSARS,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/clSARS.rda")
mfuzz.plot(eset = standSARS,
           cl = clSARS,
           mfrow = c(3, 3),
           time.labels = timeSARS)

library(Mfuzz)
cSARS2 <- 8
clSARS2 <- mfuzz(standSARS, c = cSARS2, m = mSARS)
mfuzz.plot2(eset = standSARS,
            cl = clSARS2,
            mfrow = c(3, 3),
            time.labels = timeSARS,
            xlab = "Time (h)")

mfuzz.plot2(eset = standSARS,
            cl = clSARS,
            mfrow = c(3, 3),
            time.labels = timeSARS,
            xlab = "Time (h)")

View(standSARS@assayData$exprs)

mfuzzClusters <- list()

SARSclusters <- clSARS$cluster
nClusters <- max(SARSclusters)
for(i in 1:nClusters){
  mfuzzClusters[[i]] <- names(SARSclusters)[SARSclusters == i]
}
names(mfuzzClusters) <- sprintf("Cluster%s", 1:nClusters)

save(x = mfuzzClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzClusters.rda")

pvalsMatrix <- adjANOVASARS[,3:12]

for(i in 1:length(mfuzzClusters)){
  directoryName <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = directoryName,
           where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/clusterBoxplot/",
           dat = normCopy[mfuzzClusters[[i]],],
           pvalues = pvalsMatrix[mfuzzClusters[[i]],],
           groups2 = groups,
           is.ANOVA = TRUE)
}

################################################################################

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

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/aggRes.rda")
nonHuman <- c("P0DTD1",
              "P0DTC3",
              "P0DTC7",
              "P0DTC8",
              "P0DTC5",
              "P0DTC2",
              "P0DTD2",
              "P0DTC9",
              "P00761",
              "P00883",
              "P02769")

aggRes2 <- aggRes[-which(rownames(adjANOVASARS) %in% nonHuman),]
onlyHumanMuzz <- mfuzzing(aggregatedResult = aggRes2)
save(x = onlyHumanMuzz,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/onlyHumanFuzz.rda")

cOnlyHuman <- 20
clOnlyHuman <- mfuzz(eset = onlyHumanMuzz[["standarized.data"]],
                     c = cOnlyHuman,
                     m = onlyHumanMuzz[["m.fuzzification"]])
mfuzz.plot2(eset = onlyHumanMuzz[["standarized.data"]],
            cl = clOnlyHuman,
            mfrow = c(3, 2),
            time.labels = c(0, 3, 6, 9, 16),
            xlab = "Time (h)")

mfuzzClusters <- list()
SARSclusters <- clOnlyHuman$cluster
nClusters <- max(SARSclusters)
for(i in 1:nClusters){
  mfuzzClusters[[i]] <- names(SARSclusters)[SARSclusters == i]
}
names(mfuzzClusters) <- sprintf("Cluster%s", 1:nClusters)
save(x = mfuzzClusters,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mfuzzClusters.rda")

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normCopy.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")
pvalsMatrix <- adjANOVASARS[,3:12]

for(i in 1:length(mfuzzClusters)){
  directoryName <- paste("Cluster", as.character(i), sep = "")
  photoBox(is.already = FALSE,
           nameDir = directoryName,
           where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_All/clustersBoxplots/",
           dat = normCopy[mfuzzClusters[[i]],],
           pvalues = pvalsMatrix[mfuzzClusters[[i]],],
           groups2 = groups,
           is.ANOVA = TRUE)
}

adjTable <- adjANOVASARS[!(rownames(adjANOVASARS) %in% nonHuman),]
adjTable <- cbind.data.frame(description[rownames(adjTable),], adjTable)
save(x = adjTable,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjTable.rda")

library(xlsx)

for(i in 1:length(mfuzzClusters)){
  clusterSheet <- paste("Cluster", as.character(i), sep = "")
  clusterFileConstant <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/RNorm_All/"
  clusterFileVariable <- paste("cluster", as.character(i), ".xlsx", sep = "")
  clusterFile <- paste(clusterFileConstant, clusterFileVariable, sep = "")
  
  clusterData <- adjTable[mfuzzClusters[[i]],]
  clusterData <- clusterData[order(clusterData[,"pvalue.Ftest"], decreasing = FALSE),]
  write.xlsx2(x = clusterData,
              file = clusterFile,
              row.names = FALSE,
              col.names = TRUE,
              sheetName = clusterSheet,
              append = FALSE)
}



################################################################################

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/aggRes.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/standSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/mSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/clSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")

################################################################################

BiocManager::install("GOSim")
library(GOSim)
library(UniprotR)

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/refinedMatrix.rda")

