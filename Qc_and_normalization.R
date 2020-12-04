#Quality control using scater

#R initiation
suppressPackageStartupMessages({
   library(monocle3)
   library(here)
   library(scater)
   library(scran)
   library(dplyr)
   library(SingleCellExperiment)
   library(AnnotationHub)
   library(DropletUtils)
   library(pheatmap)
   })

#Data loading
#library(SingleCellExperiment)
lmmp_mf <- readRDS(here("data", "sce_combined_TH_TL.rds")) 
lmmp_mf_sce <- as.SingleCellExperiment(lmmp_mf)

#library(scater)
rownames(lmmp_mf_sce) <- uniquifyFeatureNames(
  rowData(lmmp_mf_sce)$gene_id, rowData(lmmp_mf_sce)$SYMBOL)

exprs(lmmp_mf_sce) <- log2(
  calculateCPM(lmmp_mf_sce, use.size.factors = FALSE) + 1)

#library(AnnotationHub)
ah <- AnnotationHub()
ens.mm.v97 <- ah[['AH73905']]
rowData(lmmp_mf_sce)$SEQNAME <- mapIds(ens.mm.v97, keys=rowData(lmmp_mf_sce)$gene_id,
                                   keytype="GENEID", column="SEQNAME")

#Quality Control
unfiltered <- lmmp_mf_sce
is.mito <- rowData(lmmp_mf_sce)$SEQNAME == "MT"
stats <- perCellQCMetrics(lmmp_mf_sce, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent")
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