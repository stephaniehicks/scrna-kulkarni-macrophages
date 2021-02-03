#Feature selection


sce <- readRDS(here("data", "sce_norm.rds"))
sce



# Removing lowly expressed genes, after QC
num_reads <- 1
num_cells <- 0.01*ncol(sce)
keep <- which(DelayedArray::rowSums(counts(sce) >= num_reads ) >= num_cells)
sce <- sce[keep,]
sce

#HVG selection 
library(scran)
dec.sce <- modelGeneVar(sce,
                        block=colData(sce)$sample_id)
chosen.hvgs <- getTopHVGs(dec.sce, prop = 0.1)
par(mfrow = c(1,2))
blocked.stats <- dec.sce$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main = i, pch = 16, cex = 0.5, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
  curfit <- metadata(current)
  curve(curfit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}

#To correct batch effects, use batchelor/fastMNN
library(batchelor)
sce.mnn <- fastMNN(sce, subset.row = chosen.hvgs, 
        batch = colData(sce)$sample_id, 
        BSPARAM = BiocSingular::RandomParam(deferred = TRUE))
altExp(sce, "mnn") <- sce.mnn
reducedDim(sce, "corrected") <- reducedDim(sce.mnn, "corrected")
str(reducedDim(sce.mnn))
dim(reducedDim(sce.mnn, "corrected"))
metadata(sce.mnn)$merge.info$lost.var
######TH               TL
######[1,] 0.02952356 0.04297248


#Graph based clustering, use MNN-corrected data
library(scran)
g.mnn <- buildSNNGraph(sce, use.dimred = "corrected")
clusters.mnn <- igraph::cluster_walktrap(g.mnn)$membership
table.mnn <- table(Cluster = clusters.mnn, Batch = sce.mnn$batch)
colLabels(sce) <- factor(clusters.mnn)

#TSNE plot by batch or cluster
set.seed(0010101010)
sce <- runTSNE(sce, subset_row = chosen.hvgs, dimred = "corrected")
TSNEplot0 <- plotTSNE(sce, colour_by = "sample_id")
TSNEplot0
library(scater)
TSNEplot1 <- plotReducedDim(sce, "TSNE", colour_by="label", text_by="label")
TSNEplot1

#UMAP plot by batch or cluster
sce <- runUMAP(sce, dimred = "corrected")
UMAPplot0 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "sample_id")
UMAPplot0
UMAPplot1 <- plotUMAP(sce, colour_by = "label", text_by="label")
UMAPplot1


saveRDS(sce, file = here("data", "sce_mnn_alt.rds"))
saveRDS(sce.mnn, file = here("data", "sce_mnn.rds"))

#Or proceed to dimensionality reduction w/o batch effect correction
library(Seurat)
set.seed(101010011)
rawsce <- runPCA(sce, subset_row = chosen.hvgs, 
       BSPARAM = BiocSingular::RandomParam())
rawsce <- runTSNE(rawsce, dimred = "PCA", perplexity = 40)
rawsce <- runUMAP(rawsce)

plotReducedDim(rawsce, dimred = "PCA", colour_by = "sample_id")
plotTSNE(rawsce, colour_by = "sample_id")
plotUMAP(rawsce, colour_by = "sample_id")

#clustering

#graph-based clustering
library(scran)
g <- buildSNNGraph(rawsce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

library(scater)
colLabels(rawsce) <- factor(clust)
plotReducedDim(rawsce, "PCA", colour_by="label", text_by = "label")
plotReducedDim(rawsce, "TSNE", colour_by="label", text_by = "label")
plotReducedDim(rawsce, "UMAP", colour_by="label", text_by = "label")

saveRDS(rawsce, file = here("data", "sce_raw.rds"))

#k-means clustering
#determining best K
library(cluster)
set.seed(110010101)
gaps <- clusGap(reducedDim(rawsce, "PCA"), kmeans, K.max=20)
best.k <- maxSE(gaps$Tab[,"gap"], gaps$Tab[,"SE.sim"])
best.k

plot(gaps$Tab[,"gap"], xlab="Number of clusters", ylab="Gap statistic")
abline(v=best.k, col="red")

### best.k
#####[1] 5

set.seed(100)
clust.kmeans <- kmeans(reducedDim(rawsce, "PCA"), centers=20)
table(clust.kmeans$cluster)


colLabels(rawsce) <- factor(clust.kmeans$cluster)
plotReducedDim(rawsce, "TSNE", colour_by="label", text_by = "label")


