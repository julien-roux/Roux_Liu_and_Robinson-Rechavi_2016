## Selective constraints on coding sequences of nervous system genes are a major determinant of duplicate gene retention in vertebrates
## Julien Roux, Jialin Liu and Marc Robinson-Rechavi
## August 2016
## This script aims at reproducing analyses of the paper and regenerate its figures and supplementary figures

## Install RColorBrewer package for a nice color palette
##install.packages("RColorBrewer")
library(RColorBrewer)
myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))


#### Read gene lists ####
## mouse orthologs of zebrafish 3R duplicates and singletons
dup3R_orth <- read.table("gene_lists/mouse_orthologs_3R_ohnologs.txt", sep="\t", h=T)
sing3R_orth <- read.table("gene_lists/mouse_orthologs_3R_singletons.txt", sep="\t", h=T)
names(dup3R_orth)[1] <- "Ensembl.Gene.ID"
names(sing3R_orth)[1] <- "Ensembl.Gene.ID"

## mouse small scale duplicates
SSD <- read.table("gene_lists/mouse_SSD.txt", sep="\t", h=T)
names(SSD) <- "Ensembl.Gene.ID"


#### sequence evolution files ####
## mouse-rat dN and dS values from Ensembl 79
dnds <- read.table("sequence_evolution/mouse_rat_dn_ds_e79.txt", sep="\t", h=T)
## keep only one-to-one homology type
dnds <- subset(dnds, Homology.Type == "ortholog_one2one")
## keep only genes that diverged on the Murinae branch
dnds <- subset(dnds, Ancestor == "Murinae")
## remove columns not used
dnds <- dnds[,c(1,5,6)]
## take log10 of dN and dS
dnds <- subset(dnds, dnds$dN>0)
dnds <- subset(dnds, dnds$dS>0)
dnds$dN.log10 <- log10(dnds$dN)
dnds$dS.log10 <- log10(dnds$dS)


## Akashi's test result
akashi <- read.table("sequence_evolution/mouse_akashi_test_scores.txt", sep="\t", h=T)[,c(2,4)]
names(akashi)[1]  <- "Ensembl.Gene.ID"
akashi$Psi.log10 <- log10(akashi$Psi)
## merge data frames
dndsAkashi <- merge(dnds, akashi, by="Ensembl.Gene.ID")


#### nervous system structures ####
NSbroad <- read.table("organ_systems/nervous_system_organs_broad_10090.txt", sep="\t", h=F)[,1]  
NSstrict <- read.table("organ_systems/nervous_system_organs_strict_10090.txt", sep="\t", h=F)[,1]


#### Microarray data ####
## GSE16496 experiment
GSE16496.annot <- read.table("expression_microarrays_Bgee/GSE16496_10090_matrix_annotation.tsv", sep="\t", h=T, quote="")
GSE16496.intensities <- read.table("expression_microarrays_Bgee/GSE16496_10090_matrix_intensities.tsv", sep="\t", h=T)
GSE16496.calls <- read.table("expression_microarrays_Bgee/GSE16496_10090_matrix_calls.tsv", sep="\t", h=T, quote="")
## Average probesets mapped to the same gene
GSE16496.aggregated <- aggregate(GSE16496.intensities[,-c(1,2)], list(GSE16496.intensities$Gene_ID), mean, na.rm=T)
names(GSE16496.aggregated)[1] <- "Ensembl.Gene.ID"
## Summarized expression values
meanNS <- apply(GSE16496.aggregated[,-1], 1, mean, na.rm=T)
maxNS <- apply(GSE16496.aggregated[,-1], 1, max, na.rm=T)
## Idem, but excluding  pituitary gland, retina and pineal body samples that do not behave like other NS samples
drop <- GSE16496.annot$Chip.ID[GSE16496.annot$Anatomical.entity.name %in% c("pituitary gland", "retina", "pineal body")]
meanNSstar <- apply(GSE16496.aggregated[ , !(names(GSE16496.aggregated) %in% drop)][,-1], 1, mean, na.rm=T)
maxNSstar <- apply(GSE16496.aggregated[ , !(names(GSE16496.aggregated) %in% drop)][,-1], 1, max, na.rm=T)
## Merge them in a data frame
NS_expr <- data.frame(GSE16496.aggregated[, 1], meanNS, maxNS, meanNSstar, maxNSstar)
names(NS_expr)[1] <- "Ensembl.Gene.ID"


## GSE3594 expression
GSE3594.annot <- read.table("expression_microarrays_Bgee/GSE3594_10090_matrix_annotation.tsv", sep="\t", h=T, quote="")
GSE3594.intensities <-  read.table("expression_microarrays_Bgee/GSE3594_10090_matrix_intensities.tsv", sep="\t", h=T)
GSE3594.calls <-  read.table("expression_microarrays_Bgee/GSE3594_10090_matrix_calls.tsv", sep="\t", h=T)
## Average probesets mapped to the same gene
GSE3594.aggregated <- aggregate(GSE3594.intensities[,-c(1,2)], list(GSE3594.intensities$Gene_ID), mean, na.rm=T)
names(GSE3594.aggregated)[1] <- "Ensembl.Gene.ID"
## Summarized expression values across all tissues
meanAllTissues <- apply(GSE3594.aggregated[,-1], 1, mean, na.rm=T)
maxAllTissues <- apply(GSE3594.aggregated[,-1], 1, max, na.rm=T)
medianAllTissues <- apply(GSE3594.aggregated[,-1], 1, median, na.rm=T)
meanNS <- apply(GSE3594.aggregated[,-1][GSE3594.annot$Anatomical.entity.ID %in% NSbroad], 1, mean, na.rm=T)
maxNS <- apply(GSE3594.aggregated[,-1][GSE3594.annot$Anatomical.entity.ID %in% NSbroad], 1, max, na.rm=T)
## Merge them in a data frame
all_expr <- data.frame(GSE3594.aggregated[,1], meanAllTissues, maxAllTissues, medianAllTissues, meanNS, maxNS)
names(all_expr)[1] <- "Ensembl.Gene.ID"


## GSE10246 expression
GSE10246.annot <- read.table("expression_microarrays_Bgee/GSE10246_10090_matrix_annotation.tsv", sep="\t", h=T, quote="")
GSE10246.intensities <- read.table("expression_microarrays_Bgee/GSE10246_10090_matrix_intensities.tsv", sep="\t", h=T)
GSE10246.calls <- read.table("expression_microarrays_Bgee/GSE10246_10090_matrix_calls.tsv", sep="\t", h=T, quote="")
## Average probesets mapped to the same gene
GSE10246.aggregated <- aggregate(GSE10246.intensities[,-c(1,2)], list(GSE10246.intensities$Gene_ID), mean, na.rm=T)
names(GSE10246.aggregated)[1] <- "Ensembl.Gene.ID"


#### translation rates ####
stemcells <- read.table("translation_rate/mouse_ESC.txt", sep="\t", h=F)
names(stemcells) <- c("UCSC.ID", "rate")
fibroblasts  <- read.table("translation_rate/mouse_embryonic_fibroblast.txt", sep="\t", h=F)
names(fibroblasts) <- c("UCSC.ID", "rate")
neutrophils <- read.table("translation_rate/mouse_neutrophils.txt", sep="\t", h=F)
names(neutrophils) <- c("UCSC.ID", "rate")

## ensembl gene ID to UCSC gene ID mapping
ensemblID_ucscID <- read.table("translation_rate/mouse_EnsemblID_2_UCSC_ID.txt", sep="\t", h=T)
ensemblID_ucscID <- subset(ensemblID_ucscID, UCSC.ID != "")
## keep one-to-one mappings
ensemblID <- subset(as.data.frame(table(ensemblID_ucscID$Ensembl.Gene.ID)), Freq == 1)[,1]
ucscID <- subset(as.data.frame(table(ensemblID_ucscID$UCSC.ID)), Freq == 1)[,1]
ensemblID_ucscID <- ensemblID_ucscID[ensemblID_ucscID$Ensembl.Gene.ID %in% ensemblID & ensemblID_ucscID$UCSC.ID %in% ucscID,]

## Add Ensembl ID
stemcells <- merge(stemcells, ensemblID_ucscID, by="UCSC.ID")
fibroblasts <- merge(fibroblasts, ensemblID_ucscID, by="UCSC.ID")
neutrophils <- merge(neutrophils, ensemblID_ucscID, by="UCSC.ID")


#### connectivity ####
connectivity <- read.table("connectivity/OGEE_mouse_connectivity.txt", sep="\t", h=T)[,c(2,3)]
names(connectivity)[1] <- "Ensembl.Gene.ID"
## log10 transformation
connectivity$connectivity.log10 <- log10(connectivity$connectivity)


#### protein complexes ####
uniprot_heteromultimer <- read.table("mouse_complexes/UniProt/mouse_heteromultimers.txt", sep="\t", h=F, quote="")[,1]
uniprot_heterodimer <- read.table("mouse_complexes/UniProt/mouse_heterodimers.txt", sep="\t", h=F, quote="")[,1]
uniprot_undefined <- read.table("mouse_complexes/UniProt/mouse_undefined_complexes.txt", sep="\t", h=F, quote="")[,1]
uniprot_homomer <- read.table("mouse_complexes/UniProt/mouse_homomultimers.txt", sep="\t", h=F, quote="")[,1]
uniprot_monomer <- read.table("mouse_complexes/UniProt/mouse_monomers.txt", sep="\t", h=F, quote="")[,1]

## keep genes belonging to different categories in the "highest" category 
uniprot_monomer <- uniprot_monomer[!uniprot_monomer %in% unique(unlist(list(uniprot_heteromultimer, uniprot_heterodimer, uniprot_undefined, uniprot_homomer)))]
uniprot_homomer <- uniprot_homomer[!uniprot_homomer %in% unique(unlist(list(uniprot_heteromultimer, uniprot_heterodimer, uniprot_undefined)))]
uniprot_undefined <- uniprot_undefined[!uniprot_undefined %in% unique(unlist(list(uniprot_heteromultimer, uniprot_heterodimer)))]
uniprot_heterodimer <- uniprot_heterodimer[!uniprot_heterodimer %in% uniprot_heteromultimer]

## ensembl gene ID to uniprot ID mapping
ensemblID_uniprotID <- read.table("mouse_complexes/UniProt/mouse_EnsemblID_2_UniProtGeneName.txt",sep="\t",h=T)
ensemblID_uniprotID <- subset(ensemblID_uniprotID, UniProt.Gene.Name != "")
## keep one-to-one mappings
ensemblID <- subset(as.data.frame(table(ensemblID_uniprotID$Ensembl.Gene.ID)), Freq == 1)[,1]
uniprotID <- subset(as.data.frame(table(ensemblID_uniprotID$UniProt.Gene.Name)), Freq == 1)[,1]
ensemblID_uniprotID <- ensemblID_uniprotID[ensemblID_uniprotID$Ensembl.Gene.ID %in% ensemblID & ensemblID_uniprotID$UniProt.Gene.Name %in% uniprotID,]

## define non-annotated genes
uniprot_non_annotated <- setdiff(unique(ensemblID_uniprotID$Ensembl.Gene.ID), unique(c(as.vector(uniprot_heteromultimer), as.vector(uniprot_heterodimer), as.vector(uniprot_undefined), as.vector(uniprot_homomer), as.vector(uniprot_monomer))))

## Put all categories in a list
uniprot <- list(as.vector(uniprot_heteromultimer), as.vector(uniprot_heterodimer), as.vector(uniprot_undefined), as.vector(uniprot_homomer), as.vector(uniprot_monomer), uniprot_non_annotated)


## protein complexes from Corum database
corum <- read.table("mouse_complexes/Corum/allComplexes.csv", sep="\t", h=T, quote="", fill=T)
corum <- subset(corum, organism == "Mouse")
entrez <- strsplit(as.character(corum[,4]), ",")
corum[,5] <- unlist(lapply(entrez, length))
corum_formatted <- data.frame(complexID=rep(corum[,1], sapply(entrez, length)), entrezID=unlist(entrez), peptideNumbers=rep(corum[,5], sapply(entrez, length)))
corum_formatted[,2] <- gsub("\\(|\\)","", corum_formatted[,2]) 

## ensembl gene ID to entrez gene ID mapping
entrezEnsemblID <- read.table("mouse_complexes/Corum/mouse_entrezID_ensemblID_one2one.txt", sep="\t", h=F)
names(entrezEnsemblID) <- c("entrezID","Ensembl.Gene.ID")
## merge complexes file with ID conversion file
corum_formatted <- merge(corum_formatted, entrezEnsemblID, by="entrezID")

## classify genes in different categories
corum_heteromultimer <- subset(corum_formatted, peptideNumbers >= 3)
corum_heteromultimer <- unique(corum_heteromultimer$Ensembl.Gene.ID)
corum_heterodimer <- subset(corum_formatted, peptideNumbers == 2)
corum_heterodimer <- unique(corum_heterodimer$Ensembl.Gene.ID)

## keep genes belonging to different categories in the "highest" category 
corum_heterodimer <- corum_heterodimer[!corum_heterodimer %in% corum_heteromultimer]

## define non-annotated genes
corum_non_annotated <- setdiff(unique(entrezEnsemblID$Ensembl.Gene.ID), unique(c(as.vector(corum_heteromultimer), as.vector(corum_heterodimer))))

## Put all categories in a list
corum <- list(as.vector(corum_heteromultimer), as.vector(corum_heterodimer), corum_non_annotated)


#### in situ nervous system expression ####
onlyNS <- read.table("organ_systems/genes_only_nervous_system_expression_broad_10090.txt", sep="\t", h=F)
names(onlyNS) <- "Ensembl.Gene.ID"
NS <- read.table("organ_systems/genes_nervous_system_and_other_expression_broad_10090.txt", sep="\t", h=F)
names(NS) <- "Ensembl.Gene.ID"
noNS <- read.table("organ_systems/genes_no_nervous_system_expression_broad_10090.txt", sep="\t", h=F)
names(noNS) <- "Ensembl.Gene.ID"
## get all nervous genes of mouse
allNS <- rbind(onlyNS, NS)
names(allNS) <- "Ensembl.Gene.ID"


#### bin analysis generic function ####
binAnalysis <- function (geneFeature, feature, n, duplicates, singletons) {
  ## Classify all genes into n bins
  geneFeature$bins <- cut(geneFeature[,feature], breaks=quantile(geneFeature[,feature], probs=seq(0, 1, 1/n),  na.rm=T))

  ## create data frame with proportion of orthologs of duplicate  
  bins <- unique(geneFeature$bins)
  bins <- as.data.frame(bins[!is.na(bins)])
  bins$prop <- rep(NA, nrow(bins))
  bins$median <- rep(NA, nrow(bins))
  for (i in 1:length(bins$bins)) {
    dup <- sum(geneFeature$Ensembl.Gene.ID[geneFeature$bins == bins$bins[i]] %in% duplicates)
    sing <- sum(geneFeature$Ensembl.Gene.ID[geneFeature$bins == bins$bins[i]] %in% singletons)
    bins$prop[i] <- dup/(dup+sing)
    bins$median[i] <- median(geneFeature[,feature][geneFeature$bins == bins$bins[i]], na.rm=T)
  }
  return(bins)
}

#### bin analysis of duplicate ortholog proportion according to different gene properties ####
## Separates genes into 2 different groups (e.g., NS and non-NS genes) 
binAnalysis2Groups <- function (geneFeature, group1, group2, feature, xlimits, ylimits, colors=c(myPalette[9], myPalette[3], myPalette[4]), legends=c("Nervous system genes","All genes","Non-nervous system genes")) {
  ## Classify all genes into 10 bins
  bins <- binAnalysis(geneFeature, feature, 10, dup3R_orth[,1], sing3R_orth[,1])

  ## plot
  xlabel <- gsub("\\.", " ", feature)
  plot(bins$median, bins$prop, ylim=ylimits, xlim=xlimits, col=colors[1], xlab=xlabel, ylab="Proportion of duplicate orthologs", pch=16)
  ## linear regression
  fit_prop <- lm(bins$prop ~ bins$median)
  abline(fit_prop, col=colors[1])
  
  ## genes of group 1 only
  bins <- binAnalysis(geneFeature[geneFeature$Ensembl.Gene.ID %in% group1$Ensembl.Gene.ID, ], feature, 10, dup3R_orth[,1], sing3R_orth[,1])
  ## add to plot
  points(bins$median, bins$prop, ylim=ylimits, col=colors[2], pch=15)
  ## linear regression
  fit_prop_group1 <- lm(bins$prop ~ bins$median)
  abline(fit_prop_group1, col=colors[2])
  
  ## genes of group 2 only
  bins <- binAnalysis(geneFeature[geneFeature$Ensembl.Gene.ID %in% group2$Ensembl.Gene.ID, ], feature, 10, dup3R_orth[,1], sing3R_orth[,1])
  ## add to plot
  points(bins$median, bins$prop, ylim=ylimits, col=colors[3], pch=17)
  ## linear regression
  fit_prop_group2 <- lm(bins$prop ~ bins$median)
  abline(fit_prop_group2, col=colors[3])

  ## legends
  legend("topleft", legends, col=colors[c(2,1,3)], pch=c(15,16,17), bty="n")
  legend("topright", c(
                       paste0("Beta=", signif(summary(fit_prop_group1)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop_group1)$coefficients[2,4], 2)),
                       paste0("Beta=", signif(summary(fit_prop)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop)$coefficients[2,4], 2)),
                       paste0("Beta=", signif(summary(fit_prop_group2)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop_group2)$coefficients[2,4], 2))
                       ), lty=1, col=colors[c(2,1,3)], bty="n")
}

## studying the effect of one gene property on WGDs retention and then correcting by another gene property
binAnalysis20Percent <- function(geneFeature, feature1, feature2, duplicates, singletons, xlimits, ylimits, ylabel="Proportion of duplicate orthologs", colors=c("green", "red", "blue"), logx=F){
  if(logx==TRUE){
    untf <- TRUE
    log <- "x"
  } else {
    untf <- FALSE
    log <- ""
  }
  xlabel <- gsub("\\.", " ", feature1)
  
  bins <- binAnalysis(geneFeature, feature1, 10, duplicates, singletons)
  plot(bins$median, bins$prop, ylim=ylimits, xlim=xlimits, col=colors[1], xlab=xlabel, ylab=ylabel, pch=16, log=log)

  ## linear regression 
  fit_prop <- lm(bins$prop ~ bins$median)
  abline(fit_prop, col=colors[1], untf=untf)
  
  ## 20% highest feature2 genes
  bins <- binAnalysis(geneFeature[geneFeature[,feature2] >= quantile(geneFeature[,feature2], probs=0.8, na.rm=T), ], feature1, 10, duplicates, singletons)
  ## add to plot
  points(bins$median, bins$prop, ylim=ylimits, col=colors[2], pch=15)
  ## linear regression
  fit_prop_high <- lm(bins$prop ~ bins$median)
  abline(fit_prop_high, col=colors[2], untf=untf)

  ## 20% lowest feature2 genes
  bins <- binAnalysis(geneFeature[geneFeature[,feature2] <= quantile(geneFeature[,feature2], probs=0.2, na.rm=T), ], feature1, 10, duplicates, singletons)
  ## add to plot
  points(bins$median, bins$prop, ylim=ylimits, col=colors[3], pch=17)
  ## linear regression
  fit_prop_low <- lm(bins$prop ~ bins$median)
  abline(fit_prop_low, col=colors[3], untf=untf)

  ## legends
  legend("topleft", c(paste0("20% highest ", feature2), paste0("all ", feature2), paste0("20% lowest ", feature2)), col=colors[c(2,1,3)], pch=c(15,16,17), bty="n")  
  ## Add the regression line, and the p-value
  legend("topright", c(
                       paste0("Beta=", signif(summary(fit_prop_high)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop_high)$coefficients[2,4], 2)),
                       paste0("Beta=", signif(summary(fit_prop)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop)$coefficients[2,4], 2)),
                       paste0("Beta=", signif(summary(fit_prop_low)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop_low)$coefficients[2,4], 2))
                       ), lty=1, col=colors[c(2,1,3)], bty="n")
}

#### correlation ploting ####
plotCorrelation <- function(x, y, xlabel, ylabel, color, title, xlimit, ylimit){
  plot(x, y, xlab=xlabel, ylab=ylabel, pch=16, col=color, main=title, xlim=xlimit, ylim=ylimit)
  rho <- cor.test(x, y, method="spearman", exact=F)
  ## Add rho and p-value
  legend("topright", paste0("Rho=", signif(rho$estimate, 2), " / p=", signif(rho$p.value, 2)), col="black", bty="n")
  ## Add loess lines
  loess_fit <- loess(y ~ x)
  xl <- seq(xlimit[1], xlimit[2], by=(xlimit[2]-xlimit[1])/100)
  lines(xl, predict(loess_fit, xl), col="red", lwd=2)
}

#####***** Figure 1 *****#####
## Fig 1A
## Add duplication info to data frame
dndsAkashi$WGD <- NA
dndsAkashi[dndsAkashi$Ensembl.Gene.ID %in% dup3R_orth$Ensembl.Gene.ID,]$WGD <- "dup3R_orth"
dndsAkashi[dndsAkashi$Ensembl.Gene.ID %in% sing3R_orth$Ensembl.Gene.ID,]$WGD <- "sing3R_orth"
dndsAkashi$SSD <- NA
dndsAkashi[dndsAkashi$Ensembl.Gene.ID %in% SSD$Ensembl.Gene.ID,]$SSD <- "SSD"

## reconciliate calls across probesets mapped to same gene (a bit long)
GSE3594.calls.aggregated <- aggregate(GSE3594.calls[,-c(1,2)], list(GSE3594.calls$Gene_ID), function(x){ return(paste(unique(x), collapse="/")) })
names(GSE3594.calls.aggregated)[1] <- "Ensembl.Gene.ID"
## Merge with dN, dS and Psi and duplication status
dndsAkashiCalls <- merge(dndsAkashi, GSE3594.calls.aggregated, by="Ensembl.Gene.ID")

## Calculate % WGD
GSE3594.ratios.WGD <- apply(dndsAkashiCalls[,-c(1:9)], 2, function(x){ return(sum(dndsAkashiCalls$WGD == "dup3R_orth" & x == "present", na.rm=T)/sum((dndsAkashiCalls$WGD == "dup3R_orth" | dndsAkashiCalls$WGD == "sing3R_orth") & x == "present", na.rm=T)) })

## Plot with a dot per replicate
GSE3594.ratios.WGD.merged <- cbind(GSE3594.annot[, 2:4], GSE3594.ratios.WGD)
GSE3594.ratios.WGD.mean <- aggregate(GSE3594.ratios.WGD, list(GSE3594.annot$Anatomical.entity.ID), mean)
GSE3594.ratios.WGD.mean <- merge(GSE3594.ratios.WGD.mean, unique(GSE3594.annot[,3:4]), by.x=1, by.y=1, all.x=T)
GSE3594.ratios.WGD.merged <- merge(GSE3594.ratios.WGD.merged, GSE3594.ratios.WGD.mean[,1:2], by.x=2, by.y=1)
names(GSE3594.ratios.WGD.merged)[5] <- "mean.ratio"
GSE3594.ratios.WGD.merged <- GSE3594.ratios.WGD.merged[order(GSE3594.ratios.WGD.merged$mean.ratio), ]
GSE3594.ratios.WGD.merged$order <- as.numeric(reorder(factor(GSE3594.ratios.WGD.merged$Anatomical.entity.ID), GSE3594.ratios.WGD.merged$mean.ratio))

par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE3594.ratios.WGD.merged$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
plot(GSE3594.ratios.WGD.merged$order, GSE3594.ratios.WGD.merged$GSE3594.ratios.WGD, ylab="Proportion of duplicate orthologs", xlab="", xaxt="n", type="n", ylim=c(0.107,0.16))
axis(side = 1, at = unique(GSE3594.ratios.WGD.merged[,c(3,6)])[,2], labels=F)
abline(v=unique(GSE3594.ratios.WGD.merged[,c(3,6)])[,2], col="grey", lty=3)
points(GSE3594.ratios.WGD.merged$order, GSE3594.ratios.WGD.merged$"GSE3594.ratios.WGD", col=color, pch=16, cex=1.5)
text(unique(GSE3594.ratios.WGD.merged[,c(3,6)])[,2], 0.103, srt = 45, adj = 1,labels = unique(GSE3594.ratios.WGD.merged[,c(3,6)])[,1], xpd = TRUE)

## Fig 1B
## Calculate % SSD
GSE3594.ratios.SSD <- apply(dndsAkashiCalls[,-c(1:9)], 2, function(x){ return(sum(dndsAkashiCalls$SSD == "SSD" & x == "present", na.rm=T)/sum(x == "present", na.rm=T)) })

## Plot with a dot per replicate
GSE3594.ratios.SSD.merged <- cbind(GSE3594.annot[, 2:4], GSE3594.ratios.SSD)
GSE3594.ratios.SSD.mean <- aggregate(GSE3594.ratios.SSD, list(GSE3594.annot$Anatomical.entity.ID), mean)
GSE3594.ratios.SSD.mean <- merge(GSE3594.ratios.SSD.mean, unique(GSE3594.annot[,3:4]), by.x=1, by.y=1, all.x=T)
GSE3594.ratios.SSD.merged <- merge(GSE3594.ratios.SSD.merged, GSE3594.ratios.SSD.mean[,1:2], by.x=2, by.y=1)
names(GSE3594.ratios.SSD.merged)[5] <- "mean.ratio"
GSE3594.ratios.SSD.merged <- GSE3594.ratios.SSD.merged[order(GSE3594.ratios.SSD.merged$mean.ratio), ]
GSE3594.ratios.SSD.merged$order <- as.numeric(reorder(factor(GSE3594.ratios.SSD.merged$Anatomical.entity.ID), GSE3594.ratios.SSD.merged$mean.ratio))

par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE3594.ratios.SSD.merged$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
plot(GSE3594.ratios.SSD.merged$order, GSE3594.ratios.SSD.merged$GSE3594.ratios.SSD, ylab="Proportion of small-scale duplicates", xlab="", xaxt="n", type="n", ylim=c(0.011, 0.0208))
axis(side = 1, at = unique(GSE3594.ratios.SSD.merged[,c(3,6)])[,2], labels=F)
abline(v=unique(GSE3594.ratios.SSD.merged[,c(3,6)])[,2], col="grey", lty=3)
points(GSE3594.ratios.SSD.merged$order, GSE3594.ratios.SSD.merged$"GSE3594.ratios.SSD", col=color, pch=16, cex=1.5)
text(unique(GSE3594.ratios.SSD.merged[,c(3,6)])[,2], 0.0102, srt = 45, adj = 1,labels = unique(GSE3594.ratios.SSD.merged[,c(3,6)])[,1], xpd = TRUE)


#####***** Figure 2 *****#####
## merge dnds and duplicates/singletons lists
dup3R_orth_dnds <- merge(dup3R_orth, dnds, by="Ensembl.Gene.ID")
sing3R_orth_dnds <- merge(sing3R_orth, dnds, by="Ensembl.Gene.ID")
## merge with in situ expression
dup3R_orth_allNS_dnds <- merge(dup3R_orth_dnds, allNS, by="Ensembl.Gene.ID" )
dup3R_orth_noNS_dnds <- merge(dup3R_orth_dnds, noNS, by="Ensembl.Gene.ID" )
sing3R_orth_allNS_dnds <- merge(sing3R_orth_dnds, allNS, by="Ensembl.Gene.ID" )
sing3R_orth_noNS_dnds <- merge(sing3R_orth_dnds, noNS, by="Ensembl.Gene.ID" )

allBoxes <- list(log10(dup3R_orth_dnds$dN), log10(sing3R_orth_dnds$dN), log10(dup3R_orth_allNS_dnds$dN), log10(sing3R_orth_allNS_dnds$dN), log10(dup3R_orth_noNS_dnds$dN), log10(sing3R_orth_noNS_dnds$dN))
boxplot(allBoxes,
        notch=T, pch=16, outcex=0.5, boxwex=0.7,
        ylab=expression(paste(italic(d)[N], ~(log[10]))),
        at=c(1,2,4,5,7,8), col=myPalette[1:2],
        ylim=c(-5,1), names=rep(c("Dup.", "Sing."), 3)
        )
mtext(c("All", "Nervous system", "Non-nervous system"), side=1, line=3, at=c(1.5, 4.5, 7.5))
abline(v=3, lty=2, col="grey")
text(c(1,2,4,5,7,8), -4.8, labels=paste0("n=", unlist(lapply(allBoxes, length))))

#####***** Figure S15 *****#####
## Address potential power issues (there are more NS genes than non-NS genes)
## http://biorxiv.org/content/early/2016/09/01/072959#comment-2890987439
## -> Subsampling of NS genes to number of non-NS genes 10,000 times and record wilcoxon test p-values
subsampled <- rep(NA, times=10000)
for (i in 1:10000){
  ## print(i)
  subsampled[i] <- wilcox.test(sample(allBoxes[[3]], size=length(allBoxes[[5]])), sample(allBoxes[[4]], size=length(allBoxes[[6]])))$p.value
}
summary(subsampled)
hist(subsampled, breaks=100, main="", xlab="Subsampling p-value")
## Compare to non-NS p-value
summary(subsampled < 0.05)
summary(subsampled < wilcox.test(allBoxes[[5]], allBoxes[[6]])$p.value)
abline(v=0.05, col=myPalette[1], lty=2)
abline(v=wilcox.test(allBoxes[[5]], allBoxes[[6]])$p.value, col=myPalette[1], lty=1)


#####***** Figure 3 *****#####
par(mfrow=c(2,2))
## Fig 3A
binAnalysis2Groups(dndsAkashi, allNS, noNS, "dN.log10", xlimits=c(-2.8, -0.7), ylimits=c(0, 0.35))
## Fig 3B
binAnalysis2Groups(dndsAkashi, allNS, noNS, "Psi.log10", xlimits=c(-0.2, 0.32), ylimits=c(0,0.35))
## Fig 3C
binAnalysis20Percent(dndsAkashi[dndsAkashi$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "dN.log10", "Psi", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(-3, -0.8), ylimits=c(0,0.6), colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)
## Fig 3D
binAnalysis20Percent(dndsAkashi[dndsAkashi$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "Psi.log10", "dN", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(-0.2, 0.32), ylimits=c(0,0.6), colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)


#####***** Figure S2 *****#####
## Fig S2A
## reconciliate calls across probesets mapped to same gene (a bit long)
GSE10246.calls.aggregated <- aggregate(GSE10246.calls[,-c(1,2)], list(GSE10246.calls$Gene_ID), function(x){ return(paste(unique(x), collapse="/")) })
names(GSE10246.calls.aggregated)[1] <- "Ensembl.Gene.ID"
## Merge with dN, dS and Psi and duplication status
dndsAkashiCalls <- merge(dndsAkashi, GSE10246.calls.aggregated, by="Ensembl.Gene.ID")

## Calculate % WGD
GSE10246.ratios.WGD <- apply(dndsAkashiCalls[,-c(1:9)], 2, function(x){ return(sum(dndsAkashiCalls$WGD == "dup3R_orth" & x == "present", na.rm=T)/sum((dndsAkashiCalls$WGD == "dup3R_orth" | dndsAkashiCalls$WGD == "sing3R_orth") & x == "present", na.rm=T)) })

## Plot with a dot per replicate
GSE10246.ratios.WGD.merged <- cbind(GSE10246.annot[, 2:4], GSE10246.ratios.WGD)
GSE10246.ratios.WGD.mean <- aggregate(GSE10246.ratios.WGD, list(GSE10246.annot$Anatomical.entity.ID), mean)
GSE10246.ratios.WGD.mean <- merge(GSE10246.ratios.WGD.mean, unique(GSE10246.annot[,3:4]), by.x=1, by.y=1, all.x=T)
GSE10246.ratios.WGD.merged <- merge(GSE10246.ratios.WGD.merged, GSE10246.ratios.WGD.mean[,1:2], by.x=2, by.y=1)
names(GSE10246.ratios.WGD.merged)[5] <- "mean.ratio"
GSE10246.ratios.WGD.merged <- GSE10246.ratios.WGD.merged[order(GSE10246.ratios.WGD.merged$mean.ratio), ]
GSE10246.ratios.WGD.merged$order <- as.numeric(reorder(factor(GSE10246.ratios.WGD.merged$Anatomical.entity.ID), GSE10246.ratios.WGD.merged$mean.ratio))

par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE10246.ratios.WGD.merged$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
plot(GSE10246.ratios.WGD.merged$order, GSE10246.ratios.WGD.merged$GSE10246.ratios.WGD, ylab="Proportion of duplicate orthologs", xlab="", xaxt="n", type="n", ylim=c(0.095,0.14))
axis(side = 1, at = unique(GSE10246.ratios.WGD.merged[,c(3,6)])[,2], labels=F)
abline(v=unique(GSE10246.ratios.WGD.merged[,c(3,6)])[,2], col="grey", lty=3)
points(GSE10246.ratios.WGD.merged$order, GSE10246.ratios.WGD.merged$"GSE10246.ratios.WGD", col=color, pch=16, cex=1.5)
text(unique(GSE10246.ratios.WGD.merged[,c(3,6)])[,2], 0.09, srt = 45, adj = 1,labels = unique(GSE10246.ratios.WGD.merged[,c(3,6)])[,1], xpd = TRUE)


## Fig S2B
## GTEx human data
GTEx.ratio.WGD <- read.table("expression_RNA-seq_Bgee/GTEx_SRP012682_matrix_annotation+WGD_ratio.tsv", h=T, sep="\t", comment.char="", quote = "")
GTEx.ratio.WGD$order <- as.numeric(reorder(factor(GTEx.ratio.WGD$uberonId), GTEx.ratio.WGD$mean.ratio))

## Here there are a lot of points so we plot as a boxplot
par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(unique(GTEx.ratio.WGD$uberonId) %in% NSbroad, myPalette[3], myPalette[4])
plot(GTEx.ratio.WGD$order, GTEx.ratio.WGD$ratio.WGD, ylab="Proportion of duplicate orthologs", xlab="", xaxt="n", type="n", ylim=c(0.1,0.135))
abline(v=unique(GTEx.ratio.WGD$order), col="grey", lty=3)
boxplot(GTEx.ratio.WGD$ratio.WGD ~ GTEx.ratio.WGD$order, ylab="", xlab="", col=color, pch=16, cex=0.5, names=NA, add=T)
text(unique(GTEx.ratio.WGD$order), 0.097, srt = 45, adj = 1,labels = unique(GTEx.ratio.WGD$uberonName), xpd = TRUE)


## Fig S2C
## reconciliate calls across probesets mapped to same gene (a bit long)
GSE16496.calls.aggregated <- aggregate(GSE16496.calls[,-c(1,2)], list(GSE16496.calls$Gene_ID), function(x){ return(paste(unique(x), collapse="/")) })
names(GSE16496.calls.aggregated)[1] <- "Ensembl.Gene.ID"
## Merge with dN, dS and Psi and duplication status
dndsAkashiCalls <- merge(dndsAkashi, GSE16496.calls.aggregated, by="Ensembl.Gene.ID")

## Calculate % WGD
GSE16496.ratios.WGD <- apply(dndsAkashiCalls[,-c(1:9)], 2, function(x){ return(sum(dndsAkashiCalls$WGD == "dup3R_orth" & x == "present", na.rm=T)/sum((dndsAkashiCalls$WGD == "dup3R_orth" | dndsAkashiCalls$WGD == "sing3R_orth") & x == "present", na.rm=T)) })

## Plot with a dot per replicate
GSE16496.ratios.WGD.merged <- cbind(GSE16496.annot[, 2:4], GSE16496.ratios.WGD)
GSE16496.ratios.WGD.mean <- aggregate(GSE16496.ratios.WGD, list(GSE16496.annot$Anatomical.entity.ID), mean)
GSE16496.ratios.WGD.mean <- merge(GSE16496.ratios.WGD.mean, unique(GSE16496.annot[,3:4]), by.x=1, by.y=1, all.x=T)
GSE16496.ratios.WGD.merged <- merge(GSE16496.ratios.WGD.merged, GSE16496.ratios.WGD.mean[,1:2], by.x=2, by.y=1)
names(GSE16496.ratios.WGD.merged)[5] <- "mean.ratio"
GSE16496.ratios.WGD.merged <- GSE16496.ratios.WGD.merged[order(GSE16496.ratios.WGD.merged$mean.ratio), ]
GSE16496.ratios.WGD.merged$order <- as.numeric(reorder(factor(GSE16496.ratios.WGD.merged$Anatomical.entity.ID), GSE16496.ratios.WGD.merged$mean.ratio))

par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE16496.ratios.WGD.merged$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
plot(GSE16496.ratios.WGD.merged$order, GSE16496.ratios.WGD.merged$GSE16496.ratios.WGD, ylab="Proportion of duplicate orthologs", xlab="", xaxt="n", type="n", ylim=c(0.113,0.142))
axis(side = 1, at = unique(GSE16496.ratios.WGD.merged[,c(3,6)])[,2], labels=F)
abline(v=unique(GSE16496.ratios.WGD.merged[,c(3,6)])[,2], col="grey", lty=3)
points(GSE16496.ratios.WGD.merged$order, GSE16496.ratios.WGD.merged$"GSE16496.ratios.WGD", col=color, pch=16, cex=1.5)
text(unique(GSE16496.ratios.WGD.merged[,c(3,6)])[,2], 0.11, srt = 45, adj = 1,labels = unique(GSE16496.ratios.WGD.merged[,c(3,6)])[,1], xpd = TRUE)

#####***** Figure S3 *****#####
## Fig S3A
## Merge with dN, dS and Psi
dndsAkashiExpr <- merge(dndsAkashi, GSE3594.aggregated, by="Ensembl.Gene.ID")

## Calculate the correlation of each sample with dN
allCors <- apply(dndsAkashiExpr[,-c(1:9)], 2, function(x){ return(cor.test(dndsAkashiExpr$dN, x, method="spearman", exact=F)$estimate) })
GSE3594.cors <- cbind(GSE3594.annot[, 2:4], allCors)
GSE3594.cors <- GSE3594.cors[order(GSE3594.cors[,4], decreasing=T), ]
## order samples by mean value across replicates
GSE3594.cors$order <- length(unique(GSE3594.cors$Anatomical.entity.ID)) - as.numeric(reorder(factor(GSE3594.cors$Anatomical.entity.ID), GSE3594.cors$allCors, mean)) + 1

## Mean and Max expression
corMeanAllTissues <- cor.test(dndsAkashiExpr$dN, all_expr$meanAllTissues[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxAllTissues <- cor.test(dndsAkashiExpr$dN, all_expr$maxAllTissues[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMeanNS <- cor.test(dndsAkashiExpr$dN, all_expr$meanNS[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxNS <- cor.test(dndsAkashiExpr$dN, all_expr$maxNS[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate

## Merge with other correlation coefficients
GSE3594.cors$Anatomical.entity.name <- as.character(GSE3594.cors$Anatomical.entity.name)
GSE3594.cors <- rbind(GSE3594.cors, c(NA, NA, "mean expression", corMeanAllTissues, 36), c(NA, NA, "maximum expression", corMaxAllTissues, 37), c(NA, NA, "mean nervous system expression", corMeanNS, 38), c(NA, NA, "maximum nervous system expression", corMaxNS, 39))

## Plot with one dot per replicate
par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE3594.cors$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
color[(length(color)-3):(length(color)-2)] <- rep("black", times=2)
color[(length(color)-1):length(color)] <- rep(myPalette[3], times=2)
plot(GSE3594.cors$order, GSE3594.cors$allCors, ylim=c(-0.32,-0.05), ylab="Spearman's correlation coefficient", xlab="", xaxt="n", type="n")
axis(side = 1, at = unique(GSE3594.cors[,c(3,5)])[,2], labels=F)
abline(v=unique(GSE3594.cors[,c(3,5)])[,2], col="grey", lty=3)
points(GSE3594.cors$order, GSE3594.cors$allCors, col=color, pch=16, cex=1.5)
text(unique(GSE3594.cors[,c(3,5)])[,2], -0.345, srt = 45, adj = 1,labels = unique(GSE3594.cors[,c(3,5)])[,1], xpd = TRUE)

## Fig S3B
## Merge with dN, dS and Psi
dndsAkashiExpr <- merge(dndsAkashi, GSE16496.aggregated, by="Ensembl.Gene.ID")

## Calculate the correlation of each sample with dN
allCors <- apply(dndsAkashiExpr[,-c(1:9)], 2, function(x){ return(cor.test(dndsAkashiExpr$dN, x, method="spearman", exact=F)$estimate) })
GSE16496.cors <- cbind(GSE16496.annot[, 2:4], allCors)
GSE16496.cors <- GSE16496.cors[order(GSE16496.cors[,4], decreasing=T), ]
## order samples by mean value across replicates
GSE16496.cors$order <- length(unique(GSE16496.cors$Anatomical.entity.ID)) - as.numeric(reorder(factor(GSE16496.cors$Anatomical.entity.ID), GSE16496.cors$allCors, mean)) + 1

## Mean and Max expression
corMeanNS <- cor.test(dndsAkashiExpr$dN, NS_expr$meanNS[NS_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxNS <- cor.test(dndsAkashiExpr$dN, NS_expr$maxNS[NS_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMeanNSstar <- cor.test(dndsAkashiExpr$dN, NS_expr$meanNSstar[NS_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxNSstar <- cor.test(dndsAkashiExpr$dN, NS_expr$maxNSstar[NS_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate

## Merge with other correlation coefficients
GSE16496.cors$Anatomical.entity.name <- as.character(GSE16496.cors$Anatomical.entity.name)
GSE16496.cors <- rbind(GSE16496.cors, c(NA, NA, "mean nervous system expression", corMeanNS, 48), c(NA, NA, "maximum nervous system expression", corMaxNS, 49), c(NA, NA, "mean nervous system expression*", corMeanNSstar, 50), c(NA, NA, "maximum nervous system expression*", corMaxNSstar, 51))

## Plot with one dot per replicate
par(mar=c(15, 5, 2, 1) + 0.1) # increase margins c(bottom, left, top, right)'
## color: dark green for strict NS, light green for broad & !strict NS, violet for non-NS
color <- ifelse(GSE16496.cors$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
color[(length(color)-3):length(color)] <- rep(myPalette[3], times=2)
plot(GSE16496.cors$order, GSE16496.cors$allCors, ylim=c(-0.39,-0.31), ylab="Spearman's correlation coefficient", xlab="", xaxt="n", type="n")
axis(side = 1, at = unique(GSE16496.cors[,c(3,5)])[,2], labels=F)
abline(v=unique(GSE16496.cors[,c(3,5)])[,2], col="grey", lty=3)
points(GSE16496.cors$order, GSE16496.cors$allCors, col=color, pch=16, cex=1.5)
text(unique(GSE16496.cors[,c(3,5)])[,2], -0.396, srt = 45, adj = 1,labels = unique(GSE16496.cors[,c(3,5)])[,1], xpd = TRUE)


#####***** Figure S4 *****#####
## Merge with dN, dS and nervous system expression summaries 
dnds_NS_expr <- merge(dnds, NS_expr, by="Ensembl.Gene.ID")

## merge Akashi test scores and nervous system expression summaries 
akashi_NS_expr <-  merge(akashi, NS_expr, by="Ensembl.Gene.ID")

## merge connectivity with dN, dS / Akashi test scores / nervous system expression summaries
dnds_connectivity <- merge(dnds, connectivity, by="Ensembl.Gene.ID")
NS_expr_connectivity <- merge(NS_expr, connectivity, by="Ensembl.Gene.ID" )
akashi_connectivity <- merge(akashi, connectivity, by="Ensembl.Gene.ID" )

par(mfrow=c(3,3))
## Fig S4A: dN vs. mean NS expression
plotCorrelation(dnds_NS_expr$dN.log10, dnds_NS_expr$meanNS, "dN (log10)", "Mean nervous system expression (log2)", rgb(0,0,0,0.1), "", c(-4.5, 1.5), c(2.17, 16))
## Fig S4B: dN vs. max NS expression
plotCorrelation(dnds_NS_expr$dN.log10, dnds_NS_expr$maxNS, "dN (log10)", "Maximum nervous system expression (log2)", rgb(0,0,0,0.1), "", c(-4.5, 1.5), c(2.17, 16))
## Fig S4C: dN vs. dS
plotCorrelation(dnds$dN.log10, dnds$dS.log10, "dN (log10)", "dS (log10)", rgb(0,0,0,0.1), "", c(-4.5, 1.5), c(-2.5,2.5))
## Fig S4D: dS vs. max NS expression
plotCorrelation(dnds_NS_expr$dS.log10, dnds_NS_expr$maxNS, "dS (log10)", "Maximum nervous system expression (log2)", rgb(0,0,0,0.1), "", c(-2, 1), c(2.17, 16))
## Fig S4E: dN vs. Psi
plotCorrelation(dndsAkashi$dN.log10, dndsAkashi$Psi.log10, "dN (log10)", "Psi (log10)", rgb(0,0,0,0.1), "", c(-4.5, 1.5), c(-0.5,0.5))
## Fig S4F:  Psi vs. max NS expression
plotCorrelation(akashi_NS_expr$Psi.log10, akashi_NS_expr$maxNS, "Psi (log10)", "Maximum nervous system expression (log2)", rgb(0,0,0,0.1), "", c(-0.5, 0.7), c(2.17, 16))
## Fig S4G:  dN vs. connectivity
plotCorrelation(dnds_connectivity$dN.log10, dnds_connectivity$connectivity.log10, "dN (log10)", "Connectivity (log10)", rgb(0,0,0,0.1), "", c(-4.5, 1.5), c(0,3))
##  Fig S4H: Psi vs. connectivity
plotCorrelation(akashi_connectivity$Psi.log10, akashi_connectivity$connectivity.log10, "Psi (log10)", "Connectivity (log10)", rgb(0,0,0,0.1), "",c(-0.4, 0.5),c(0,3))
##  Fig S4I: max NS expression vs. connectivity
plotCorrelation(NS_expr_connectivity$maxNS, NS_expr_connectivity$connectivity.log10, "Maximum nervous system expression (log2)", "Connectivity (log10)", rgb(0,0,0,0.1), " ",c(2.17, 16),c(0,3))


#####***** Figure S5 *****#####
## Merge with dndsAkashi with NS_expr
dndsAkashiExpr <- merge(dndsAkashi, NS_expr, by="Ensembl.Gene.ID")
## filter many genes with low expression levels
dndsAkashiExpr <- subset(dndsAkashiExpr, meanNS >= 3)

## Fig S5A
binAnalysis2Groups(dndsAkashiExpr, allNS, noNS, "meanNS", xlimits=c(2, 14), ylimits=c(0, 0.6))
binAnalysis2Groups(dndsAkashiExpr, allNS, noNS, "maxNS", xlimits=c(2, 14), ylimits=c(0, 0.6))
## Only green lines used in Fig S5A

## Fig S5B and S5C
binAnalysis20Percent(dndsAkashiExpr[dndsAkashiExpr$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "maxNS", "dN", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(4,15), ylimits=c(0, 0.6), ylabel="Proportion of duplicate orthologs", colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)
binAnalysis20Percent(dndsAkashiExpr[dndsAkashiExpr$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "dN.log10", "maxNS", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(-3,-0.5), ylimits=c(0, 0.6), ylabel="Proportion of duplicate orthologs", colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)


#####***** Figure S6 *****#####
## Fig S6A
## Merge expression levels with dN, dS and Psi
dndsAkashiExpr <- merge(dndsAkashi, GSE3594.aggregated, by="Ensembl.Gene.ID")

## Calculate the correlation of each sample with dS
allCors <- apply(dndsAkashiExpr[,-c(1:9)], 2, function(x){ return(cor.test(dndsAkashiExpr$dS, x, method="spearman", exact=F)$estimate) })
GSE3594.cors <- cbind(GSE3594.annot[, 2:4], allCors)
GSE3594.cors <- GSE3594.cors[order(GSE3594.cors[,4], decreasing=T), ]
## order samples by mean value across replicates
GSE3594.cors$order <- length(unique(GSE3594.cors$Anatomical.entity.ID)) - as.numeric(reorder(factor(GSE3594.cors$Anatomical.entity.ID), GSE3594.cors$allCors, mean)) + 1

## Mean and Max expression
corMeanAllTissues <- cor.test(dndsAkashiExpr$dS, all_expr$meanAllTissues[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxAllTissues <- cor.test(dndsAkashiExpr$dS, all_expr$maxAllTissues[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMeanNS <- cor.test(dndsAkashiExpr$dS, all_expr$meanNS[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate
corMaxNS <- cor.test(dndsAkashiExpr$dS, all_expr$maxNS[all_expr$Ensembl.Gene.ID %in% dndsAkashiExpr$Ensembl.Gene.ID], method="spearman", exact=F)$estimate

## Merge with other correlation coefficients
GSE3594.cors$Anatomical.entity.name <- as.character(GSE3594.cors$Anatomical.entity.name)
GSE3594.cors <- rbind(GSE3594.cors, c(NA, NA, "mean expression", corMeanAllTissues, 36), c(NA, NA, "maximum expression", corMaxAllTissues, 37), c(NA, NA, "mean nervous system expression", corMeanNS, 38), c(NA, NA, "maximum nervous system expression", corMaxNS, 39))

## Plot with one dot per replicate
par(mar=c(15, 5, 2, 1) + 0.1)
color <- ifelse(GSE3594.cors$Anatomical.entity.ID %in% NSbroad, myPalette[3], myPalette[4])
color[(length(color)-3):(length(color)-2)] <- rep("black", times=2)
color[(length(color)-1):length(color)] <- rep(myPalette[3], times=2)
plot(GSE3594.cors$order, GSE3594.cors$allCors, ylim=c(-0.13,0.001), ylab="Spearman's correlation coefficient", xlab="", xaxt="n", type="n")
axis(side = 1, at = unique(GSE3594.cors[,c(3,5)])[,2], labels=F)
abline(v=unique(GSE3594.cors[,c(3,5)])[,2], col="grey", lty=3)
points(GSE3594.cors$order, GSE3594.cors$allCors, col=color, pch=16, cex=1.5)
text(unique(GSE3594.cors[,c(3,5)])[,2], -0.14, srt = 45, adj = 1,labels = unique(GSE3594.cors[,c(3,5)])[,1], xpd = TRUE)

## Fig S6B
binAnalysis2Groups(dndsAkashi, allNS, noNS, "dS.log10", xlimits=c(-1.0, -0.4), ylimits=c(0, 0.35))


#####***** figure S7 *****#####
## merge translation rate, Psi and ensemblUcscID
stemcellsAkashi <- merge(stemcells,     akashi, by="Ensembl.Gene.ID")
fibroblastsAkashi <- merge(fibroblasts, akashi, by="Ensembl.Gene.ID")
neutrophilsAkashi <- merge(neutrophils, akashi, by="Ensembl.Gene.ID")

par(mfrow=c(1,3))
## psi vs. translation rate neutrophils
plotCorrelation(log10(neutrophilsAkashi$Psi), log10(neutrophilsAkashi$rate), "Psi (log10)", "Translation rate (log10)", rgb(0,0,0,0.1), "Neutrophils", c(-0.5, 0.7), c(0.79, 0.85))
## psi vs. translation rate Embryonic fibroblasts
plotCorrelation(log10(fibroblastsAkashi$Psi), log10(fibroblastsAkashi$rate), "Psi (log10)", "Translation rate (log10)", rgb(0,0,0,0.1), "Embryonic fibroblasts", c(-0.5, 0.7), c(0.85, 0.91))
## psi vs. translation rate Embryonic stem cells
plotCorrelation(log10(stemcellsAkashi$Psi), log10(stemcellsAkashi$rate), "Psi (log10)", "Translation rate (log10)", rgb(0,0,0,0.1), "Embryonic stem cells", c(-0.5, 0.7), c(0.52, 0.62))


#####***** Figure S8 *****#####

## Fig S8A
## Duplicates and singletons in each categories
uniprotDupOrth <- lapply(uniprot, function(x) x[x %in% dup3R_orth[,1]])
uniprotSingOrth <- lapply(uniprot, function(x) x[x %in% sing3R_orth[,1]])
uniprotDupOrthNS <- lapply(uniprotDupOrth, function(x) length(x[x %in% allNS[,1]]))
uniprotSingOrthNS <- lapply(uniprotSingOrth, function(x) length(x[x %in% allNS[,1]]))
uniprotDupOrthNoNS <- lapply(uniprotDupOrth, function(x) length(x[x %in% noNS[,1]]))
uniprotSingOrthNoNS <- lapply(uniprotSingOrth, function(x) length(x[x %in% noNS[,1]]))
## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:6)) {
  propNS <- cbind(propNS, uniprotDupOrthNS[[n]] / (uniprotDupOrthNS[[n]] + uniprotSingOrthNS[[n]]))
  numberNS <- cbind(numberNS, (uniprotDupOrthNS[[n]] + uniprotSingOrthNS[[n]]))
  propNoNS <- cbind(propNoNS, uniprotDupOrthNoNS[[n]] / (uniprotDupOrthNoNS[[n]] + uniprotSingOrthNoNS[[n]]))
  numberNoNS <- cbind(numberNoNS, (uniprotDupOrthNoNS[[n]] + uniprotSingOrthNoNS[[n]]))
}
## barplot
par(mar=c(12,5,2,2))
mp <- barplot(c(propNS,propNoNS), las=2, main="", ylab="Proportion of duplicate orthologs", ylim=c(0,0.3), col=rep(c(myPalette[3], myPalette[4]), each=6))
xlabels <- c("Hetero-multimer genes", "Hetero-dimer genes","Uncharacterized complexes genes","Homo-multimer genes","Monomer genes","Non-annotated genes")
text(x=mp, y=-0.01, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.01, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[6])/2, (mp[7]+mp[12])/2), y=0.29, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))

## Fig S8B
## merge protein complexes with dnds
uniprotdnds <- lapply(uniprot, function(x) dnds[dnds[,1] %in% x,])

colors <- myPalette[c(5, 7, 8, 9, 10, 12)]
legends <- vector()
for (i in 1:length(uniprotdnds)) {
  bins <- binAnalysis(uniprotdnds[[i]], "dN.log10", 10, dup3R_orth[,1], sing3R_orth[,1])
  if (i == 1){
    plot(bins$median, bins$prop, ylim=c(0, 1), xlim=c(-2.9, -0.7), col=colors[i], main="", xlab="dN (log10)", ylab="Proportion of duplicate orthologs", pch=i)
  } else {
    points(bins$median, bins$prop, pch=i,col=colors[i])
  }
  fit_prop <- lm(bins$prop ~ bins$median)
  abline(fit_prop, col=colors[i])
  legends <- c(legends, paste0("Beta=", signif(summary(fit_prop)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop)$coefficients[2,4], 2)))
}
legend("topleft", xlabels, col=colors, pch=c(1:6), bty="n")
legend("topright", legends, lty=1, col=colors, bty="n")


#####***** Figure S9 *****#####
metabolic <- read.table("GO/mouse_metabolic", sep="\t",h=T)
names(metabolic)[1] <- "UniProt.Gene.Name"
metabolic <- merge(metabolic, ensemblID_uniprotID, by="UniProt.Gene.Name")$Ensembl.Gene.ID
metabolic <- list(as.vector(metabolic), setdiff(unique(ensemblID_uniprotID$Ensembl.Gene.ID), metabolic))

## Duplicates and singletons in each categories
metabolicDupOrth <- lapply(metabolic, function(x) x[x %in% dup3R_orth[,1]])
metabolicSingOrth <- lapply(metabolic, function(x) x[x %in% sing3R_orth[,1]])
metabolicDupOrthNS <- lapply(metabolicDupOrth, function(x) length(x[x %in% allNS[,1]]))
metabolicSingOrthNS <- lapply(metabolicSingOrth, function(x) length(x[x %in% allNS[,1]]))
metabolicDupOrthNoNS <- lapply(metabolicDupOrth, function(x) length(x[x %in% noNS[,1]]))
metabolicSingOrthNoNS <- lapply(metabolicSingOrth, function(x) length(x[x %in% noNS[,1]]))

## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:length(metabolic))) {
  propNS <- cbind(propNS, metabolicDupOrthNS[[n]] / (metabolicDupOrthNS[[n]] + metabolicSingOrthNS[[n]]))
  numberNS <- cbind(numberNS, (metabolicDupOrthNS[[n]] + metabolicSingOrthNS[[n]]))
  propNoNS <- cbind(propNoNS, metabolicDupOrthNoNS[[n]] / (metabolicDupOrthNoNS[[n]] + metabolicSingOrthNoNS[[n]]))
  numberNoNS <- cbind(numberNoNS, (metabolicDupOrthNoNS[[n]] + metabolicSingOrthNoNS[[n]]))
}

## barplot
par(mar=c(9,5,2,2))
mp <- barplot(c(propNS, propNoNS), las=2, main="", ylab="Proportion of duplicate orthologs", ylim=c(0,0.3), col=rep(c(myPalette[3], myPalette[4]), each=2))
xlabels <- c("Metabolic process genes", "Non-Metabolic process genes")
text(x=mp, y=-0.01, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.01, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[2])/2, (mp[3]+mp[4])/2), y=0.29, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))


#####***** Figure S10 *****#####
## Fig S10A
dup3R_orth_connectivity <- merge(dup3R_orth, connectivity, by="Ensembl.Gene.ID")
sing3R_orth_connectivity <- merge(sing3R_orth, connectivity, by="Ensembl.Gene.ID")
dup3R_orth_allNS_connectivity <- merge(dup3R_orth_connectivity, allNS, by="Ensembl.Gene.ID" )
dup3R_orth_noNS_connectivity <- merge(dup3R_orth_connectivity, noNS, by="Ensembl.Gene.ID" )
sing3R_orth_allNS_connectivity <- merge(sing3R_orth_connectivity, allNS, by="Ensembl.Gene.ID" )
sing3R_orth_noNS_connectivity <- merge(sing3R_orth_connectivity, noNS, by="Ensembl.Gene.ID" )

allBoxes <- list(dup3R_orth_connectivity$connectivity.log10, sing3R_orth_connectivity$connectivity.log10, dup3R_orth_allNS_connectivity$connectivity.log10, sing3R_orth_allNS_connectivity$connectivity.log10, dup3R_orth_noNS_connectivity$connectivity.log10, sing3R_orth_noNS_connectivity$connectivity.log10)
boxplot(allBoxes,
        notch=T, pch=16, outcex=0.5, boxwex=0.7,
        ylab=expression(paste(Connectivity~(log[10]))),
        at=c(1,2,4,5,7,8), col=myPalette[1:2],
        ylim=c(-1,4), names=rep(c("Dup.", "Sing."), 3)
        )
mtext(c("All", "Nervous system", "Non-nervous system"), side=1, line=3, at=c(1.5, 4.5, 7.5))
abline(v=3, lty=2, col="grey")
text(c(1,2,4,5,7,8), -0.8, labels=paste0("n=", unlist(lapply(allBoxes, length))))

## Fig S10B and S10C
binAnalysis20Percent(dnds_connectivity[dnds_connectivity$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "dN.log10", "connectivity", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(-2.7, -0.8), ylimits=c(0,0.6), colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)
binAnalysis20Percent(dnds_connectivity[dnds_connectivity$Ensembl.Gene.ID %in% allNS$Ensembl.Gene.ID, ], "connectivity.log10", "dN", dup3R_orth[,1], sing3R_orth[,1], xlimits=c(0.2, 2.2), ylimits=c(0,0.6), colors=c(myPalette[3], myPalette[1], myPalette[2]), logx=F)

## Fig S10D
## merge protein complexes with connectivity
uniprotConnectivity <- lapply(uniprot, function(x) connectivity[connectivity$Ensembl.Gene.ID %in% x,])

xlabels <- c("Hetero-multimer genes", "Hetero-dimer genes","Uncharacterized complexes genes","Homo-multimer genes","Monomer genes","Non-annotated genes")
colors <- myPalette[c(5, 7, 8, 9, 10, 12)]
legends <- vector()
for (i in 1:length(uniprotConnectivity)) {
  bins <- try(binAnalysis(uniprotConnectivity[[i]], "connectivity.log10", 10, dup3R_orth[,1], sing3R_orth[,1]), silent=T)
  if (class(bins) == "try-error") {
    print("Not enough levels in studied feature to make 10 bins, trying with 9")
    bins <- try(binAnalysis(uniprotConnectivity[[i]], "connectivity.log10", 9, dup3R_orth[,1], sing3R_orth[,1]))
  }
  if (i == 1){
    plot(bins$median, bins$prop, ylim=c(0, 1), xlim=c(0.2, 2.2), col=colors[i], main="", xlab="Connecticity (log10)", ylab="Proportion of duplicate orthologs", pch=i)
  } else {
    points(bins$median, bins$prop, pch=i,col=colors[i])
  }
  fit_prop <- lm(bins$prop ~ bins$median)
  abline(fit_prop, col=colors[i])
  legends <- c(legends, paste0("Beta=", signif(summary(fit_prop)$coefficients[2,1], 2), " / p=", signif(summary(fit_prop)$coefficients[2,4], 2)))
}
legend("topleft", xlabels, col=colors, pch=c(1:6), bty="n")
legend("topright", legends, lty=1, col=colors, bty="n")


#####***** Figure S11 *****#####
dnds$SSD <- NA
dnds[dnds$Ensembl.Gene.ID %in% SSD$Ensembl.Gene.ID,]$SSD <- "SSD"

## Fig S11A
SSDOther_dn <- list(dnds$dN.log10[dnds$SSD %in% "SSD"],
                    dnds$dN.log10[is.na(dnds$SSD)],
                    dnds$dN.log10[dnds$SSD %in% "SSD" & dnds$Ensembl.Gene.ID %in% allNS[,1]],
                    dnds$dN.log10[is.na(dnds$SSD) & dnds$Ensembl.Gene.ID %in% allNS[,1]],
                    dnds$dN.log10[dnds$SSD %in% "SSD" & dnds$Ensembl.Gene.ID %in% noNS[,1]],
                    dnds$dN.log10[is.na(dnds$SSD) & dnds$Ensembl.Gene.ID %in% noNS[,1]])
                    
boxplot(SSDOther_dn,
        notch=T, pch=16, outcex=0.5, boxwex=0.7,
        ylab=expression(paste(italic(d)[N], ~(log[10]))),
        at=c(1,2,4,5,7,8),
        col=myPalette[1:2],
        ylim=c(-4.5,2),
        names=rep(c("Dup.", "Other."), 3)
        )
mtext(c("All", "Nervous system", "Non-nervous system"), side=1, line=3, at=c(1.5, 4.5, 7.5))
abline(v=3, lty=2, col="grey")
text(c(1,2,4,5,7,8), -4.3, labels=paste0("n=", unlist(lapply(SSDOther_dn, length))))

## Fig S11B and C
binAnalysis20Percent(dndsAkashi, "dN", "Psi", SSD[,1], dndsAkashi$Ensembl.Gene.ID[!dndsAkashi$Ensembl.Gene.ID %in% SSD$Ensembl.Gene.ID], xlimits=c(0.002,0.2), ylimits=c(0,0.15), ylabel="Proportion of small-scale duplicates", colors=c("black", myPalette[1], myPalette[2]), logx=T)
binAnalysis20Percent(dndsAkashi, "Psi.log10", "dN", SSD[,1], dndsAkashi$Ensembl.Gene.ID[!dndsAkashi$Ensembl.Gene.ID %in% SSD$Ensembl.Gene.ID], xlimits=c(-0.2,0.3), ylimits=c(0,0.15), ylabel="Proportion of small-scale duplicates", colors=c("black", myPalette[1], myPalette[2]), logx=F)


#####***** Figure S12 *****#####
## Fig S12A and B
dndsAkashiExpr <- merge(dndsAkashi, all_expr, by="Ensembl.Gene.ID")

binAnalysis20Percent(dndsAkashiExpr, "meanAllTissues", "dN", SSD[,1], dndsAkashiExpr$Ensembl.Gene.ID[!dndsAkashiExpr$Ensembl.Gene.ID %in% SSD[,1]], xlimits=c(2,11), ylimits=c(0, 0.3), ylabel="Proportion of duplicate orthologs", colors=c("black", myPalette[1], myPalette[2]), logx=F)
binAnalysis20Percent(dndsAkashiExpr, "maxAllTissues", "dN", SSD[,1], dndsAkashiExpr$Ensembl.Gene.ID[!dndsAkashiExpr$Ensembl.Gene.ID %in% SSD[,1]], xlimits=c(2,11), ylimits=c(0, 0.3), ylabel="Proportion of duplicate orthologs", colors=c("black", myPalette[1], myPalette[2]), logx=F)
## Only black lines of above figures shown in Figure S12A

## Fig S12C
binAnalysis20Percent(dndsAkashiExpr, "dN", "meanAllTissues", SSD[,1], dndsAkashiExpr$Ensembl.Gene.ID[!dndsAkashiExpr$Ensembl.Gene.ID %in% SSD[,1]], xlimits=c(0.001,0.12), ylimits=c(0, 0.3), ylabel="Proportion of duplicate orthologs", colors=c("black", myPalette[1], myPalette[2]), logx=T)


#####***** Figure S13 *****#####
## SSD genes in each categories
uniprotSSD <- lapply(uniprot, function(x) x[x %in% SSD[,1]])
uniprotSSDNS <- lapply(uniprotSSD, function(x) length(x[x %in% allNS[,1]]))
uniprotSSDNoNS <- lapply(uniprotSSD, function(x) length(x[x %in% noNS[,1]]))
uniprotNS <- lapply(uniprot, function(x) length(x[x %in% allNS[,1]]))
uniprotNoNS <- lapply(uniprot, function(x) length(x[x %in% noNS[,1]]))

## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:6)) {
  propNS <- cbind(propNS, uniprotSSDNS[[n]] / uniprotNS[[n]])
  numberNS <- cbind(numberNS, uniprotNS[[n]])
  propNoNS <- cbind(propNoNS, uniprotSSDNoNS[[n]] / uniprotNoNS[[n]])
  numberNoNS <- cbind(numberNoNS, uniprotNoNS[[n]])
}
## barplot
par(mar=c(12,5,2,2))
mp <- barplot(c(propNS,propNoNS), las=2, main="", ylab="Proportion of small scale duplicates", ylim=c(0,0.1), col=rep(c(myPalette[3], myPalette[4]), each=6))
xlabels <- c("Hetero-multimer genes", "Hetero-dimer genes","Uncharacterized complexes genes","Homo-multimer genes","Monomer genes","Non-annotated genes")
text(x=mp, y=-0.005, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.005, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[6])/2, (mp[7]+mp[12])/2), y=0.095, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))


#####***** Figure S14 *****#####
pValue <- read.table("FDR/all_pvalues.txt",sep="\t",h=F)
hist(pValue[,2], breaks=100, xlab="Linear regression p-values", xlim=c(0,1), main="Histogram of linear regression p-values")
abline(v=0.05, lty=2, col="red")

## calculate q-values and FDR
##source("http://bioconductor.org/biocLite.R")
##biocLite("qvalue")
library(qvalue)
fdr <- data.frame(p.adjust(pValue[,2], method="BH"))
qValue <- qvalue(pValue[,2])$qvalues


#####***** Figure S16 *****#####
## Fig S16A
## retrieve duplicates and singletons from each categories
corumDupOrth <- lapply(corum, function(x) x[x %in% dup3R_orth[,1]])
corumSingOrth <- lapply(corum, function(x) x[x %in% sing3R_orth[,1]])
corumDupOrthNS <- lapply(corumDupOrth, function(x) length(x[x %in% allNS[,1]]))
corumSingOrthNS <- lapply(corumSingOrth, function(x) length(x[x %in% allNS[,1]]))
corumDupOrthNoNS <- lapply(corumDupOrth, function(x) length(x[x %in% noNS[,1]]))
corumSingOrthNoNS <- lapply(corumSingOrth, function(x) length(x[x %in% noNS[,1]]))

## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:3)) {
  propNS <- cbind(propNS, corumDupOrthNS[[n]] / (corumDupOrthNS[[n]] + corumSingOrthNS[[n]]))
  numberNS <- cbind(numberNS, (corumDupOrthNS[[n]] + corumSingOrthNS[[n]]))
  propNoNS <- cbind(propNoNS, corumDupOrthNoNS[[n]] / (corumDupOrthNoNS[[n]] + corumSingOrthNoNS[[n]]))
  numberNoNS <- cbind(numberNoNS, (corumDupOrthNoNS[[n]] + corumSingOrthNoNS[[n]]))
}
## barplot
par(mar=c(9,5,2,2))
mp <- barplot(c(propNS,propNoNS), las=2, main="", ylab="Proportion of duplicate orthologs", ylim=c(0,0.3), col=rep(c(myPalette[3], myPalette[4]), each=3))
xlabels <- c(">= 3 genes", "2 genes","Non-annotated genes")
text(x=mp, y=-0.01, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.01, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[3])/2, (mp[4]+mp[6])/2), y=0.29, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))


## Fig S16B
## retrieve SSD 
corumSSD <- lapply(corum, function(x) x[x %in% SSD[,1]])
corumSSDNS <- lapply(corumSSD, function(x) length(x[x %in% allNS[,1]]))
corumSSDNoNS <- lapply(corumSSD, function(x) length(x[x %in% noNS[,1]]))
corumNS <- lapply(corum, function(x) length(x[x %in% allNS[,1]]))
corumNoNS <- lapply(corum, function(x) length(x[x %in% noNS[,1]]))

## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:3)) {
  propNS <- cbind(propNS, corumSSDNS[[n]] / corumNS[[n]])
  numberNS <- cbind(numberNS, corumNS[[n]])
  propNoNS <- cbind(propNoNS, corumSSDNoNS[[n]] / corumNoNS[[n]])
  numberNoNS <- cbind(numberNoNS, corumNoNS[[n]])
}
## barplot
par(mar=c(9,5,2,2))
mp <- barplot(c(propNS,propNoNS), las=2, main="", ylab="Proportion of small scale duplicates", ylim=c(0,0.04), col=rep(c(myPalette[3], myPalette[4]), each=3))
text(x=mp, y=-0.001, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.001, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[3])/2, (mp[4]+mp[6])/2), y=0.038, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))


#####***** Figure S17 *****#####
## Fig S17A
membrane <- read.table("GO/mouse_membrane", sep="\t",h=T)
names(membrane)[1] <- "UniProt.Gene.Name"
membrane <- merge(membrane, ensemblID_uniprotID, by="UniProt.Gene.Name")$Ensembl.Gene.ID
membrane <- list(as.vector(membrane), setdiff(unique(ensemblID_uniprotID$Ensembl.Gene.ID), membrane))

## Duplicates and singletons in each categories
membraneDupOrth <- lapply(membrane, function(x) x[x %in% dup3R_orth[,1]])
membraneSingOrth <- lapply(membrane, function(x) x[x %in% sing3R_orth[,1]])
membraneDupOrthNS <- lapply(membraneDupOrth, function(x) length(x[x %in% allNS[,1]]))
membraneSingOrthNS <- lapply(membraneSingOrth, function(x) length(x[x %in% allNS[,1]]))
membraneDupOrthNoNS <- lapply(membraneDupOrth, function(x) length(x[x %in% noNS[,1]]))
membraneSingOrthNoNS <- lapply(membraneSingOrth, function(x) length(x[x %in% noNS[,1]]))

## calculate dup orth proportion
propNS <-  vector()
propNoNS <- vector()
numberNS <- vector()
numberNoNS <- vector()
for(n in c(1:length(membrane))) {
  propNS <- cbind(propNS, membraneDupOrthNS[[n]] / (membraneDupOrthNS[[n]] + membraneSingOrthNS[[n]]))
  numberNS <- cbind(numberNS, (membraneDupOrthNS[[n]] + membraneSingOrthNS[[n]]))
  propNoNS <- cbind(propNoNS, membraneDupOrthNoNS[[n]] / (membraneDupOrthNoNS[[n]] + membraneSingOrthNoNS[[n]]))
  numberNoNS <- cbind(numberNoNS, (membraneDupOrthNoNS[[n]] + membraneSingOrthNoNS[[n]]))
}

## barplot
par(mar=c(9,5,2,2))
mp <- barplot(c(propNS, propNoNS), las=2, main="", ylab="Proportion of duplicate orthologs", ylim=c(0,0.3), col=rep(c(myPalette[3], myPalette[4]), each=2))
xlabels <- c("Membrane genes", "Non-membrane genes")
text(x=mp, y=-0.01, srt=45, cex.lab=0.8, adj=1, labels=rep(xlabels, 2), xpd=TRUE)
text(x=mp, y=0.01, cex=0.8, labels=paste0("n=", c(numberNS, numberNoNS)))
text(x=c((mp[1]+mp[2])/2, (mp[3]+mp[4])/2), y=0.29, cex=1, labels=c("Nervous system genes", "Non-nervous system genes"))

## Fig S17B
allMembrane <- data.frame(membrane[[1]])
names(allMembrane) <- "Ensembl.Gene.ID"
noMembrane <- data.frame(membrane[[2]])
names(noMembrane) <- "Ensembl.Gene.ID"
binAnalysis2Groups(dnds, allMembrane, noMembrane, "dN.log10", xlimits=c(-2.8, -0.7), ylimits=c(0, 0.6), colors=c(myPalette[9], myPalette[1], "black"), legends=c("Membrane genes","All genes","Non-membrane genes"))


