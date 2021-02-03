# Automatic cell type annotation with SingleR for reference

sce <- readRDS(here("data", "sce_raw.rds")) ##Notice object

#Append batch label to barcode to avoid duplicates
sce$barcode.1 <- paste(colnames(sce), sce$sample_id, sep = "_")
colnames(sce) <- sce$barcode.1

#Use SingleR with mouse ref genes
library(SingleR)
library(pheatmap)
ref <- celldex::MouseRNAseqData()
test <- sce

pred <- SingleR(test, ref, labels = ref$label.main)
table(pred$labels)
plotScoreHeatmap(pred, show.labels = TRUE,  annotation_col=
  data.frame(
  row.names=rownames(pred)))

tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(sce))
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))

#Use MetaNeighbor for automatic annotation
library(MetaNeighbor)
data(sce)
data(GOmouse)
AUROC_scores = MetaNeighbor(dat = sce,
                            experiment_labels = as.numeric(factor(sce$sample_id)),
                            celltype_labels = metadata(colData(sce))[["colLabels"]],
                            genesets = GOmouse,
                            bplot = TRUE)
#Marker Gene detection

###Combine every testing regime for comparison
combined <- multiMarkerStats(
  t=findMarkers(sce, direction="up"),
  wilcox=findMarkers(sce, test="wilcox", direction="up"),
  binom=findMarkers(sce, test="binom", direction="up")
)
colnames(combined[["1"]])
head(combined[["1"]][,1:20])


###Looking for any differences, any pairwise comparison
library(scran)
markers0 <- findMarkers(sce)
markers0
chosen16 <- "15"
interesting <- markers0[[chosen16]]
interesting[1:10,1:4]

cluster5.markers <- FindMarkers(sce, ident.1 = 5, ident.2 = c(7, 13), 
                                min.pct = 0.25)
head(cluster5.markers, n = 5)

best.set <- interesting[interesting$Top <= 10,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))

###Find cluster specific markers, upregulated, some PVAL.
markers.up1 <- findMarkers(
  sce, pval.type="some", direction="up")
chosen <- "22"
interesting.up1 <- markers.up1[[chosen]]
interesting.up1[1:10,1:3]


###Find cluster specific markers, upregulated, some PVAL, batch effects blocked
markers.up2 <- findMarkers(sce, pval.type= "some" , block=sce$sample_id,
                     direction="up")
chosen <- "20"
interesting.up2 <- markers.up2[[chosen]]
interesting.up2[1:10,1:3]

