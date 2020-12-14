#Quality control using scater

#R initiation
suppressPackageStartupMessages({
  
  library(here)
  library(scater)
  library(scran)
  library(dplyr)
  library(SingleCellExperiment)
  library(AnnotationHub)
  library(DropletUtils)
})

#Data loading
#library(SingleCellExperiment)
lmmp_mf_sce <- readRDS(here("data", "sce_combined_TH_TL.rds")) 

#library(scater)
rownames(lmmp_mf_sce) <- uniquifyFeatureNames(
  rowData(lmmp_mf_sce)$gene_id, rowData(lmmp_mf_sce)$SYMBOL)

exprs(lmmp_mf_sce) <- log2(
  calculateCPM(lmmp_mf_sce, use.size.factors = FALSE) + 1)

rowRanges.out <- rowRanges(lmmp_mf_sce)
rowData(lmmp_mf_sce)$chr.loc <- rowRanges.out@seqnames
rowData(lmmp_mf_sce)$chr.loc <- as.character(rowData(lmmp_mf_sce)$chr.loc)

#Examine and exclude expty droplets with Dropletutils
#library(DropletUtils)
set.seed(100)
e.out <- emptyDrops(counts(lmmp_mf_sce))

summary(e.out$FDR <= 0.001)
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)

set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(lmmp_mf_sce), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value", main="", col="grey80") 

lmmp_mf_sce <- lmmp_mf_sce[,which(e.out$FDR <= 0.001)]

bcrank <- barcodeRanks(counts(lmmp_mf_sce))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2)
abline(h=metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
abline(h=metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "dodgerblue"), lty = 2, cex = 1.2)
dev.off()

#Quality Control
is.mito <- rowData(lmmp_mf_sce)$chr.loc == "chrM"
stats <- perCellQCMetrics(lmmp_mf_sce, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch = lmmp_mf_sce$sample_id)
colSums(as.matrix(qc))

colData(lmmp_mf_sce) <- cbind(colData(lmmp_mf_sce), stats)
lmmp_mf_sce$treatment <- factor(lmmp_mf_sce$sample_id)
lmmp_mf_sce$discard <- qc$discard


discard.mito <- isOutlier(lmmp_mf_sce$subsets_Mito_percent, type ="higher", log=FALSE, 
                          batch=lmmp_mf_sce$treatment)
attr(discard.mito, "thresholds")
colData(lmmp_mf_sce) <- cbind(colData(lmmp_mf_sce),discard.mito)

discard.features <- isOutlier(lmmp_mf_sce$detected, type = "lower", log = TRUE, 
                              batch = lmmp_mf_sce$treatment)
attr(discard.features, "thresholds")
colData(lmmp_mf_sce) <- cbind(colData(lmmp_mf_sce),discard.features)

discard.sum <- isOutlier(lmmp_mf_sce$sum, type = "lower", log = TRUE, 
                         batch = lmmp_mf_sce$treatment)
attr(discard.sum, "thresholds")


gridExtra::grid.arrange(
  plotColData(lmmp_mf_sce, y="sum", colour_by="discard") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(lmmp_mf_sce, x="treatment", y="detected", colour_by = "discard") 
  + scale_y_log10() + ggtitle("Detected features"),
  plotColData(lmmp_mf_sce, x="treatment", y="subsets_Mito_percent", colour_by = "discard") 
  + ggtitle("Mito percent"),
  nrow=2, ncol=2
)


unfiltered <- lmmp_mf_sce
lmmp_mf_sce <- lmmp_mf_sce[,!qc$discard]
sce.nonorm <- lmmp_mf_sce
#Normalization
#library(scran)
set.seed(1000)
clusters <- quickCluster(lmmp_mf_sce)
table(clusters)
lmmp_mf_sce <- computeSumFactors(lmmp_mf_sce, cluster = clusters)
lmmp_mf_sce <- logNormCounts(lmmp_mf_sce)
summary(sizeFactors(lmmp_mf_sce))


plot(librarySizeFactors(lmmp_mf_sce), sizeFactors(lmmp_mf_sce), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")
abline(a=0, b=1, col="red")
