## Selective constraints on coding sequences of nervous system genes as a major
## determinant of long-term retention of duplicate genes in vertebrates
## Julien Roux, Jialin Liu and Marc Robinson-Rechavi
## August 2016

## This example script allows to obtain the zebrafish anatomical
## structures significantly enriched in expression of zebrafish
## 3R ohnologs (FDR < 10%), similarly to Table 1.

## Install BgeeDB via Bioconductor. For this, R >= 3.3 is required
source("http://bioconductor.org/biocLite.R")
biocLite("BgeeDB")
## Alternative way of installing BgeeDB:
## Download project at https://github.com/BgeeDB/BgeeDB_R/archive/master.zip
## Unzip it and launch R >= 3.0
## In R: install.packages("./BgeeDB_R-master", repos = NULL, type="source")
## (Installation of dependencies tidyR and data.table might be required)

## Load the package
library(BgeeDB)

## Loading expression calls for zebrafish based on in situ data
myTopAnatData <- loadTopAnatData(species=7955, datatype="in_situ")
## Look at the data
lapply(myTopAnatData, head)

## Prepare the lists of foreground and background gene lists
## Unzip the supplemntary dataset gene_lists.zip
myInterestingGenes <- read.table("./gene_lists/zebrafish_3R_ohnologs.txt", sep="\t", h=T)[,1]
## The background is composed of all genes with in situ expression 
geneList <- factor(as.integer(names(myTopAnatData$gene2anatomy) %in% myInterestingGenes))
names(geneList) <- names(myTopAnatData$gene2anatomy)
summary(geneList)

## Prepare the modified topGO object
myTopAnatObject <-  topAnat(myTopAnatData, geneList, nodeSize = 5)

## Launch the anatomical expression enrichment test
results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')

## Display significantly enriched anatomical structures at a 10% FDR threshold
enrichment3Rohnologs <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 0.1)
## The p-values are slightly different at the end of this table compared to Table 1 in the paper (e.g., the term "sensory ganglion" is enriched here, but not in Table 1). This is due to a slight modification on how the anatomical ontology is handled in BgeeDB compared to our (older) approach, influencing the decorrelation algorithm.
