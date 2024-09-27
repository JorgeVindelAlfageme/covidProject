# Using abundance data to build a machine learning classifier, only by
# considering those proteins with a relatively high concentration in human
# serum.

# 1) Loading packages and protein descriptions

library(xlsx)
library(limma)

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

# 4) Processing numerical data

abundInputPath <- "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/input/PD_mock.xlsx"
abund <- read.xlsx2(file = abundInputPath, sheetName = "abundances")
abund <- abund[, -ncol(abund)]

colnames(abund) <- gsub(x = colnames(abund),
                        pattern = "X",
                        replacement = "",
                        fixed = TRUE)
colnames(abund) <- gsub(x = colnames(abund),
                        pattern = "..",
                        replacement = ".",
                        fixed = TRUE)
rownames(abund) <- abund[, 1]
abund <- abund[, -1]
r <- rownames(abund)
abund <- apply(X = abund, MARGIN = 2, FUN = as.numeric)
rownames(abund) <- r
abund <- log2(abund)

abundOutputPath <- "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/output/abund.rda"

save(x = abund, file = abundOutputPath)

################################################################################

paxdbFile <- "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/input/uniprot-download_2022.08.01.xlsx"
# Path to Paxdb file with data about serum protein concentration
paxdbExcel <- read.xlsx2(file = paxdbFile, sheetIndex = 1)
most300 <- paxdbExcel[1:300, "Entry"]

pdttest <- read.xlsx2(
  file = "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/input/allTTestPD.xlsx",
  sheetName = "ttest"
)

rownames(pdttest) <- pdttest[, 1]
pdttest <- pdttest[, -1]

user_defined_colnames <- c(
  "adj.pval.SARS.CoV2.t16-MOCK.t0", "adj.pval.SARS.CoV2.t3-MOCK.t0",
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

rrfPD30 <- recursiveRandomForest3(
  dataset = t(abund),
  groups = as.factor(groups),
  treeNumber = 1000,
  retentionProportion = howManyLoops(
    nfeat = nrow(abund),
    wantedLoops = 30)[1],
  chosenIterationsNumber = 50,
  withNormalization = FALSE,
  topNumber = 48,
  iterationsBest = 100,
  searchForBest = FALSE,
  howManyFeats = 30
  )

rrfPD30[["chosen.features"]]

importanceTable2_30 <- cbind.data.frame(
  description[rownames(rrfPD30[["importance.matrix"]]),],
  rrfPD30[["importance.matrix"]]
  )

write.xlsx2(
  x = importanceTable2_30,
  file = "C:/Users/usuario/Downloads/fatPaper/partOfFornax/paperCovid/allNeededData/covidProject1_3/output/importanceTable.xlsx",
  row.names = FALSE,
  col.names = TRUE,
  sheetName = "Sheet1",
  append = FALSE
  )

tableWithChosenRRFFeats <- table2[rrfPD30[["chosen.features"]],]

write.xlsx2(
  x = tableWithChosenRRFFeats,
  file = "C:/Users/jvalf/Desktop/Archivos/trabajo2/partOfFornax/paperCovid/allNeededData/covidProject1_3/output/rrfProteins.xlsx",
  row.names = FALSE,
  col.names = TRUE,
  sheetName = "Sheet1",
  append = FALSE
  )

MDS2(
  rfObj = rrfPD30[["best.rf"]],
  classCol = groups,
  areChosenColors = TRUE,
  chosenColors = RColorBrewer::brewer.pal(
    n = length(unique(groups)),
    name = "Set1"
    ),
  main2 = "Multidimensional (MDS) plot using\nbest 30 classification candidates",
  legPos = "topright",
  pointCex = 3,
  legSize = 1.5,
  cexMain = 1.75,
  cexLabels = 1.6,
  cexAxes = 1.8
  )

MDSplotter3D(randomF = rrfPD30[["best.RF"]], groupsVector = groups)
PCAGif(nameMe = "rrfPD30_3D")

PDData30 <- abund[rrfPD30[["chosen.features"]],]

spareCols <- colnames(PDData30)
spareCols <- gsub(pattern = ".CoV2", replacement = "", fixed = TRUE, x = spareCols)

colnames(PDData30) <- spareCols

heat(dataF = PDData30, colD = FALSE, ycex = 0.9, xcex = 0.7)

photoBox2(is.already = FALSE,
          nameDir = "boxplots",
          where = "C:/Users/jvalf/Desktop/Archivos/trabajo2/partOfFornax/paperCovid/allNeededData/covidProject1_3/output/",
          dat = abund[rrfPD30[["chosen.features"]],],
          pvalues = pdttest[rrfPD30[["chosen.features"]],],
          groups2 = groups,
          is.ANOVA = TRUE,
          main3 = reducedDescription(fullDescription = description)[rrfPD30[["chosen.features"]]])