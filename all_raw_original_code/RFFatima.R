# Old code lines:
# load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
# load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normSARS.rda")
# load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
# load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")

load(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/rda/normSARS.rda")
load(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/rda/adjANOVASARS.rda")
load(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/rda/SARSgroups.rda")

dataTable <- cbind.data.frame(description[rownames(adjANOVASARS),], adjANOVASARS)
SARSProteins <- c("P0DTD1", "P0DTC3", "P0DTC7", "P0DTC8", "P0DTC5", "P0DTC2",
                  "P0DTD2", "P0DTC9")
errorProteins <- c("P00761", "P00883", "P02769")
sigProteins <- rownames(adjANOVASARS)[adjANOVASARS[,"qvalue.Ftest"] <= 0.05]

# Old code line:
# plasmaPath <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/plasma.txt"

plasmaPath <- "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/7_Mfuzz/plasma.txt"
plasmaData <- read.delim(file = plasmaPath)
plasmaProteins <- plasmaData[, 1]
library(randomForest)

set.seed(1)

fatimaProteins <- c("P01023", "P05121", "Q13201", "Q08380", "P05067")
fatimaData <- normSARS[fatimaProteins,]
fatimaData <- t(fatimaData)
fatimaRF <- randomForest(x = fatimaData,
                         y = as.factor(groups),
                         ntree = 1000,
                         proximity = TRUE,
                         importance = TRUE)
fatimaRF
MDS2(rfObj = fatimaRF,
     classCol = groups,
     main2 = "Multidimensional (MDS) plot using selected plasma proteins")
MDSplotter3D(randomF = fatimaRF, groupsVector = groups)
PCAGif()

fatimaData <- t(fatimaData)
spareCols <- gsub(pattern = ".CoV2", replacement = "", fixed = TRUE, x = colnames(fatimaData))
originalCols <- colnames(fatimaData)
colnames(fatimaData) <- spareCols
heat(dataF = fatimaData, colD = FALSE, ycex = 1.1, xcex = 0.7)

pvalMat <- adjANOVASARS[,3:12]

photoBox(is.already = TRUE,
         where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/R_fatima/boxplots/",
         dat = fatimaData,
         pvalues = pvalMat[fatimaProteins,],
         groups2 = groups,
         is.ANOVA = TRUE)

write.xlsx2(x = dataTable[fatimaProteins,],
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/R_fatima/fatimaProteins.xlsx",
            sheetName = "Sheet1",
            col.names = TRUE,
            row.names = FALSE,
            append = FALSE)

################################################################################

sigData <- normSARS[sigProteins,]
sigData <- sigData[-which(rownames(sigData) %in% SARSProteins),]
sigData <- sigData[intersect(x = rownames(sigData), y = plasmaProteins),]
rrfSARS <- recursiveRandomForest3(dataset = t(sigData),
                                  groups = as.factor(groups),
                                  treeNumber = 1000,
                                  retentionProportion = howManyLoops(nfeat = nrow(sigData),
                                                                     wantedLoops = 30)[1],
                                  chosenIterationsNumber = 50,
                                  withNormalization = FALSE,
                                  topNumber = 50,
                                  iterationsBest = 1000,
                                  searchForBest = FALSE,
                                  howManyFeats = 15)

rrfSARS[["best.RF"]]
importanceTable <- cbind.data.frame(description[rownames(rrfSARS[["importance.matrix"]]),],
                                    rrfSARS[["importance.matrix"]])
write.xlsx2(x = importanceTable,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/R_15/importanceTable.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)
o <- dataTable[rrfSARS[["chosen.features"]],]
o <- o[order(o[,"pvalue.Ftest"], decreasing = FALSE),]
write.xlsx2(x = o,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/R_15/rrfProteins.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)

MDS2(rfObj = rrfSARS[["best.RF"]],
     classCol = groups,
     main2 = "Multidimensional (MDS) plot using best 15 classification candidates",
     legPos = "topright")
MDSplotter3D(randomF = rrfSARS[["best.RF"]], groupsVector = groups)
PCAGif(nameMe = "rrfSARS3D")

data15 <- normSARS[rrfSARS[["chosen.features"]],]
colnames(data15) <- spareCols
heat(dataF = data15, colD = FALSE, ycex = 1, xcex = 0.7)

photoBox(is.already = FALSE,
         nameDir = "boxplots",
         where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/R_15/",
         dat = normSARS[rrfSARS[["chosen.features"]],],
         pvalues = pvalMat[rrfSARS[["chosen.features"]],],
         groups2 = groups,
         is.ANOVA = TRUE)

################################################################################

abund <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/PD_mock.xlsx",
                    sheetName = "abundances")
abund <- abund[,-ncol(abund)]
colnames(abund) <- gsub(x = colnames(abund), pattern = "X", replacement = "", fixed = TRUE)
colnames(abund) <- gsub(x = colnames(abund), pattern = "..", replacement = ".", fixed = TRUE)
rownames(abund) <- abund[,1]
abund <- abund[,-1]
r <- rownames(abund)
abund <- apply(X = abund, MARGIN = 2, FUN = as.numeric)
rownames(abund) <- r
abund <- log2(abund)
save(x = abund,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/abund.rda")

pdttest <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/allTTestPD.xlsx",
                      sheetName = "ttest")
rownames(pdttest) <- pdttest[,1]
pdttest <- pdttest[,-1]
colnames(pdttest) <- colnames(adjANOVASARS[,3:12])
r <- rownames(pdttest)
pdttest <- apply(X = pdttest, MARGIN = 2, FUN = as.numeric)
rownames(pdttest) <- r
table2 <- cbind.data.frame(description[rownames(pdttest),], pdttest)

abund <- abund[intersect(rownames(pdttest), plasmaProteins),]

rrfPD <- recursiveRandomForest3(dataset = t(abund),
                                groups = as.factor(groups),
                                treeNumber = 1000,
                                retentionProportion = howManyLoops(nfeat = nrow(abund),
                                                                   wantedLoops = 30)[1],
                                chosenIterationsNumber = 50,
                                withNormalization = FALSE,
                                topNumber = 50,
                                iterationsBest = 1000,
                                searchForBest = FALSE,
                                howManyFeats = 15)

rrfPD[["best.RF"]]
importanceTable2 <- cbind.data.frame(description[rownames(rrfPD[["importance.matrix"]]),],
                                     rrfPD[["importance.matrix"]])
write.xlsx2(x = importanceTable2,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_15/importanceTable.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)
o <- table2[rrfPD[["chosen.features"]],]
write.xlsx2(x = o,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_15/rrfProteins.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)

MDS2(rfObj = rrfPD[["best.RF"]],
     classCol = groups,
     main2 = "Multidimensional (MDS) plot using best 15 classification candidates",
     legPos = "topright")
MDSplotter3D(randomF = rrfPD[["best.RF"]], groupsVector = groups)
PCAGif(nameMe = "rrfPD3D")

PDData15 <- abund[rrfPD[["chosen.features"]],]
colnames(PDData15) <- spareCols
heat(dataF = PDData15, colD = FALSE, ycex = 1, xcex = 0.7)

photoBox(is.already = FALSE,
         nameDir = "boxplots",
         where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_15/",
         dat = abund[rrfPD[["chosen.features"]],],
         pvalues = pdttest[rrfPD[["chosen.features"]],],
         groups2 = groups,
         is.ANOVA = TRUE)

################################################################################

# 30 prots.

# 1) Loading and preprocessing data.

library(xlsx)

# Old code line
# abund <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/PD_mock.xlsx",
#                     sheetName = "abundances")

abund <- read.xlsx2(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/7_Mfuzz/PD_mock.xlsx",
                    sheetName = "abundances")

abund <- abund[, -ncol(abund)]
colnames(abund) <- gsub(x = colnames(abund), pattern = "X", replacement = "", fixed = TRUE)
colnames(abund) <- gsub(x = colnames(abund), pattern = "..", replacement = ".", fixed = TRUE)
rownames(abund) <- abund[, 1]
abund <- abund[, -1]
r <- rownames(abund)
abund <- apply(X = abund, MARGIN = 2, FUN = as.numeric)
rownames(abund) <- r
abund <- log2(abund)

save(x = abund,
     file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/abund.rda")

pdttest <- read.xlsx2(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/7_Mfuzz/allTTestPD.xlsx",
                      sheetName = "ttest")

rownames(pdttest) <- pdttest[,1]
pdttest <- pdttest[,-1]
colnames(pdttest) <- colnames(adjANOVASARS[, 3:12])
r <- rownames(pdttest)
pdttest <- apply(X = pdttest, MARGIN = 2, FUN = as.numeric)
rownames(pdttest) <- r
table2 <- cbind.data.frame(description[rownames(pdttest),], pdttest)

abund <- abund[intersect(rownames(pdttest), plasmaProteins), ]

set.seed(1)
rrfPD30 <- recursiveRandomForest3(dataset = t(abund),
                                  groups = as.factor(groups),
                                  treeNumber = 1000,
                                  retentionProportion = howManyLoops(nfeat = nrow(abund),
                                                                     wantedLoops = 30)[1],
                                  chosenIterationsNumber = 50,
                                  withNormalization = FALSE,
                                  topNumber = 50,
                                  iterationsBest = 100,
                                  searchForBest = FALSE,
                                  howManyFeats = 30)

rrfPD30[["best.RF"]]
importanceTable2_30 <- cbind.data.frame(description[rownames(rrfPD30[["importance.matrix"]]),],
                                        rrfPD30[["importance.matrix"]])
write.xlsx2(x = importanceTable2_30,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_30/importanceTable.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)
o <- table2[rrfPD30[["chosen.features"]],]
write.xlsx2(x = o,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_30/rrfProteins.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)

MDS2(rfObj = rrfPD30[["best.RF"]],
     classCol = groups,
     main2 = "Multidimensional (MDS) plot using best 30 classification candidates",
     legPos = "bottom",
     areChosenColors = TRUE,
     chosenColors = RColorBrewer::brewer.pal(n = length(unique(groups)),
                                             name = "Set1"),
     pointCex = 1.75,
     legSize = 1.5,
     cexMain = 1.7,
     cexLabels = 1.6,
     cexAxes = 1.6)

MDSplotter3D(randomF = rrfPD30[["best.RF"]], groupsVector = groups)

PCAGif(nameMe = "rrfPD30_3D")

PDData30 <- abund[rrfPD30[["chosen.features"]],]
# p <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_30/rrfProteins.xlsx",
#                 sheetIndex = 1)[,1]
# PDData30 <- abund[p,]
# colnames(PDData30) <- gsub(pattern = ".CoV2", replacement = "", x = colnames(PDData30), fixed = TRUE)
colnames(PDData30) <- spareCols
heat(dataF = PDData30, colD = FALSE, ycex = 0.7, xcex = 0.7)
heat(dataF = PDData30, colD = TRUE, ycex = 0.7, xcex = 0.7)

photoBox2(is.already = FALSE,
          nameDir = "boxplots",
          where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/PD_30/",
          dat = abund[rrfPD30[["chosen.features"]],],
          pvalues = pdttest[rrfPD30[["chosen.features"]],],
          groups2 = groups,
          is.ANOVA = TRUE,
          main3 = reducedDescription(fullDescription = description)[rrfPD30[["chosen.features"]]])

################################################################################

library(xlsx)
# paxdbFile <- "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/paxdb/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.08.01-11.12.53.87.xlsx"
paxdbFile <- "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/8_plasmaRF/paxdb/data/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.08.01-11.12.53.87.xlsx"

paxdbExcel <- read.xlsx2(file = paxdbFile, sheetIndex = 1)
most300 <- paxdbExcel[1:300, "Entry"]

load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/abund.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/adjANOVASARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/normCopy.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/descriptionSARS.rda")
load(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/rda/SARSgroups.rda")

# pdttest <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/7_Mfuzz/allTTestPD.xlsx",
#                       sheetName = "ttest")
pdttest <- read.xlsx2(file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/fatima/7_Mfuzz/allTTestPD.xlsx",
                      sheetName = "ttest")

rownames(pdttest) <- pdttest[, 1]
pdttest <- pdttest[, -1]
# colnames(pdttest) <- colnames(adjANOVASARS[, 3:12])
user_defined_colnames <- c("adj.pval.SARS.CoV2.t16-MOCK.t0", "adj.pval.SARS.CoV2.t3-MOCK.t0",
                           "adj.pval.SARS.CoV2.t6-MOCK.t0", "adj.pval.SARS.CoV2.t9-MOCK.t0",
                           "adj.pval.SARS.CoV2.t3-SARS.CoV2.t16", "adj.pval.SARS.CoV2.t6-SARS.CoV2.t16",
                           "adj.pval.SARS.CoV2.t9-SARS.CoV2.t16", "adj.pval.SARS.CoV2.t6-SARS.CoV2.t3",
                           "adj.pval.SARS.CoV2.t9-SARS.CoV2.t3", "adj.pval.SARS.CoV2.t9-SARS.CoV2.t6")
colnames(pdttest) <- user_defined_colnames
proteome_discoverer_t_test_rownames <- rownames(pdttest)
pdttest <- apply(X = pdttest, MARGIN = 2, FUN = as.numeric)
rownames(pdttest) <- proteome_discoverer_t_test_rownames
table2 <- cbind.data.frame(description[rownames(pdttest),], pdttest)

abund <- abund[intersect(rownames(pdttest), most300),]

set.seed(1)
rrfPD30 <- recursiveRandomForest3(dataset = t(abund),
                                  groups = as.factor(groups),
                                  treeNumber = 1000,
                                  retentionProportion = howManyLoops(nfeat = nrow(abund),
                                                                     wantedLoops = 30)[1],
                                  chosenIterationsNumber = 50,
                                  withNormalization = FALSE,
                                  topNumber = 48,
                                  iterationsBest = 100,
                                  searchForBest = FALSE,
                                  howManyFeats = 30)

rrfPD30[["best.RF"]]
importanceTable2_30 <- cbind.data.frame(description[rownames(rrfPD30[["importance.matrix"]]),],
                                        rrfPD30[["importance.matrix"]])
write.xlsx2(x = importanceTable2_30,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/paxdb/importanceTable.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)
o <- table2[rrfPD30[["chosen.features"]],]
write.xlsx2(x = o,
            file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/paxdb/rrfProteins.xlsx",
            row.names = FALSE,
            col.names = TRUE,
            sheetName = "Sheet1",
            append = FALSE)

# MDS2(rfObj = rrfPD30[["best.RF"]],
#      classCol = groups,
#      main2 = "Multidimensional (MDS) plot using best 30 classification candidates",
#      legPos = "topright"
#      )

MDS2(rfObj = rrfPD30[["best.RF"]],
     classCol = groups,
     areChosenColors = TRUE,
     chosenColors = RColorBrewer::brewer.pal(n = length(unique(groups)),
                                             name = "Set1"),
     main2 = "Multidimensional (MDS) plot using best 30 classification candidates",
     legPos = "topleft",
     pointCex = 1.9,
     legSize = 1.5,
     cexMain = 1.75,
     cexLabels = 1.5,
     cexAxes = 1.5)

MDSplotter3D(randomF = rrfPD30[["best.RF"]], groupsVector = groups)
PCAGif(nameMe = "rrfPD30_3D")

PDData30 <- abund[rrfPD30[["chosen.features"]],]

spareCols <- colnames(PDData30)
spareCols <- gsub(pattern = ".CoV2", replacement = "", fixed = TRUE, x = spareCols)

colnames(PDData30) <- spareCols

# q <- read.xlsx2(file = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/paxdb/rrfProteins.xlsx",
#                 sheetIndex = 1)[,1]
# PDData30 <- abund[q,]
# colnames(PDData30) <- gsub(pattern = ".CoV2", replacement = "", x = colnames(PDData30), fixed = TRUE)
heat(dataF = PDData30, colD = FALSE, ycex = 0.7, xcex = 0.7)
heat(dataF = PDData30, colD = TRUE, ycex = 0.7, xcex = 0.7)

photoBox2(is.already = FALSE,
          nameDir = "boxplots",
          where = "C:/Users/Centaurus2/Desktop/bioinf/datos/fatima/8_plasmaRF/paxdb/",
          dat = abund[rrfPD30[["chosen.features"]],],
          pvalues = pdttest[rrfPD30[["chosen.features"]],],
          groups2 = groups,
          is.ANOVA = TRUE,
          main3 = reducedDescription(fullDescription = description)[rrfPD30[["chosen.features"]]])

