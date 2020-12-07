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

#library(AnnotationHub)
ah <- AnnotationHub()
ens.mm.v97 <- ah[['AH73905']]
rowData(lmmp_mf_sce)$chr.loc <- mapIds(ens.mm.v97, keys=rownames(lmmp_mf_sce),
                                   keytype="GENEID", column="SEQNAME")

#Quality Control

#Exclude expty droplets with Dropletutils
#Drawing Kneeplot, unfiltered
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
pdf(file = here("plots", "Knee_plot_0.pdf"))
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2)
abline(h=metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
abline(h=metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "dodgerblue"), lty = 2, cex = 1.2)
dev.off()


unfiltered <- lmmp_mf_sce
is.mito <- rowData(lmmp_mf_sce)$chr.loc == "MT"
stats <- perCellQCMetrics(lmmp_mf_sce, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch = lmmp_mf_sce$treatment)
lmmp_mf_sce <- lmmp_mf_sce[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent", 
              colour_by="discard") + ggtitle("Mito percent"),
  ncol=2
)
colSums(as.matrix(qc))

#Normalization
#library(scran)
set.seed(1000)
clusters <- quickCluster(lmmp_mf_sce)
table(clusters)
lmmp_mf_sce <- computeSumFactors(lmmp_mf_sce, clusters=clusters)
lmmp_mf_sce <- logNormCounts(lmmp_mf_sce)
summary(sizeFactors(lmmp_mf_sce))

plot(librarySizeFactors(lmmp_mf_sce), sizeFactors(lmmp_mf_sce), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")

