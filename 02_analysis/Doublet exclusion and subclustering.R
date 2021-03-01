#Doublet detection with clusters
library(scDblFinder)
dbl.out <- findDoubletClusters(sce)
dbl.out
library(scater)
chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
chosen.doublet
#[1] "20"
markersd <- findMarkers(sce, direction="up")
dbl.markers <- markersd[[chosen.doublet]]
chosend <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(sce, order_columns_by="label", features=chosend, 
            center=TRUE, symmetric=TRUE, zlim=c(-5, 5)) 
            

#By simulation
library(BiocSingular)
set.seed(100)
dbl.dens <- computeDoubletDensity(sce, subset.row=chosen.hvgs, 
                                  d=ncol(reducedDim(sce)))
summary(dbl.dens)
sce$DoubletScore <- dbl.dens
plotTSNE(sce, colour_by="DoubletScore", text_by = "label")
plotColData(sce, x="label", y="DoubletScore", colour_by="label")

saveRDS(sce, file = here("data", "sce 0203.rds"))


#Subset macrophages (Cluster 7,10,14,19,22,23)
sce.mf <- sce[,sce$label %in% c(7, 10, 14, 19, 22, 23)]

dec.scemf <- modelGeneVar(sce.mf,
                          +                         block=colData(sce.mf)$sample_id)
chosen.hvgsmf <- getTopHVGs(dec.scemf, prop = 0.1)
par(mfrow = c(1,2))
blocked.statsmf <- dec.scemf$per.block
for (i in colnames(blocked.statsmf)) {
  +     current <- blocked.statsmf[[i]]
  +     plot(current$mean, current$total, main = i, pch = 16, cex = 0.5, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
  +     curfit <- metadata(current)
  +     curve(curfit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}

sce.mnnmf <- fastMNN(sce.mf, subset.row = chosen.hvgsmf, 
                     +                    batch = colData(sce.mf)$sample_id, 
                     +                    BSPARAM = BiocSingular::RandomParam(deferred = TRUE))
altExp(sce.mf, "mnn") <- sce.mnnmf
reducedDim(sce.mf, "corrected") <- reducedDim(sce.mnnmf, "corrected")
str(reducedDim(sce.mnnmf))
dim(reducedDim(sce.mnnmf, "corrected"))
metadata(sce.mnnmf)$merge.info$lost.var
#TH         TL
#[1,] 0.02000794 0.05620977
reducedDim(sce.mf, "corrected") <- reducedDim(sce.mnnmf, "corrected")

runTSNE(sce.mf, dimred= "corrected", subset_row = chosen.hvgsmf)
plotTSNE(sce.mf, colour_by = "sample_id")

g.mnnmf <- buildSNNGraph(sce.mf, k=20, use.dimred = "corrected")
clusters.mnnmf <- igraph::cluster_walktrap(g.mnnmf)$membership
table.mnnmf <- table(Cluster = clusters.mnnmf, Batch = sce.mnnmf$batch)
table.mnnmf
####Batch
#Cluster  TH  TL
#1 105  94
#2 239 197
#3 734 598
#4  15   6
#5   5   8

##When k=20
#Batch
#Cluster  TH  TL
#1 276 200
#2 742 625
#3  80  78

colLabels(sce.mf) <- factor(clusters.mnnmf)
table(colLabels(sce.mf))

#1    2    3    4    5 
#199  436 1332   21   13 

set.seed(101010011)
sce.mf <- runPCA(sce.mf, subset_row = chosen.hvgsmf, dimred ="corrected",
                   BSPARAM = BiocSingular::RandomParam())
sce.mf <- runTSNE(sce.mf, dimred = "PCA", perplexity = 40)
rawscemf <- runUMAP(rawscemf)
plotReducedDim(sce.mf, "PCA", colour_by="label", text_by = "label")
plotReducedDim(sce.mf, "TSNE", colour_by="label", text_by = "label")
TSNEplot1 <- plotReducedDim(sce.mf, "TSNE", colour_by="label", text_by="label")
TSNEplot1

markers.upmf <- findMarkers(sce.mf, block=sce.mf$sample_id,
                           direction="up")

chosen <- "2"
interesting.upmf <- markers.upmf[[chosen]]
interesting.upmf[1:10,1:3]



best.set <- interesting.upmf[interesting.upmf$Top <= 10,]
lfcs <- getMarkerEffects(best.set)
pheatmap(lfcs, breaks=seq(-5, 5, length.out=101))


plotExpression(sce.mf, features=c("C5ar1", "Cx3cr1", "Fcgr1", "Adgre1", "Cd63", 
                                  "Cd163", "Cd80", "Ly6c1", "Ccr2", "Timd4", "Cd4", "Mertk"), x="label", colour_by="label")

plotExpression(sce.mf, features=c("Cd14" ,"Cd24a", "Dpp4", "Runx3", "H2-Aa", "H2-Eb1", "H2-Ab1", "Ier3"), x="label", colour_by="label")
