# Quality control using scater

# R initiation
suppressPackageStartupMessages({
  library(here)
  library(scater)
  library(scran)
  library(dplyr)
  library(SingleCellExperiment)
  library(AnnotationHub)
  library(DropletUtils)
})

# Load data: SingleCellExperiment object
lmmp_mf_sce <- readRDS(here("data", "sce_combined_TH_TL.rds")) 

# make gene names unique
rownames(lmmp_mf_sce) <- scater::uniquifyFeatureNames(
  rowData(lmmp_mf_sce)$gene_id, rowData(lmmp_mf_sce)$SYMBOL)



##################################
#### Identify empty droplets #####
##################################

# Examine for empty droplets with DropletUtils
bcrank <- DropletUtils::barcodeRanks(counts(lmmp_mf_sce))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", 
     xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2)
abline(h=metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
abline(h=metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
legend("bottomleft", legend = c("Inflection", "Knee"), 
       col = c("darkgreen", "dodgerblue"), lty = 2, cex = 1.2)

# Exclude for empty droplets with DropletUtils:
#   emptyDrops() uses Monte Carlo simulations to compute p-values 
#   for the multinomial sampling transcripts from the ambient pool
set.seed(100)
e.out <- DropletUtils::emptyDrops(counts(lmmp_mf_sce))

summary(e.out$FDR <= 0.001)
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)
# > table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)
# Limited
# Sig     FALSE TRUE
# FALSE   603    0
# TRUE    437 5775

# The Limited field in the output indicates whether or not the computed  
# p-value for a particular barcode is bounded by the number of iterations. 
# If any non-significant barcodes are TRUE for Limited, we may need to 
# increase the number of iterations. A larger number of iterations will 
# result in a lower p-value for these barcodes, which may allow them to 
# be detected after correcting for multiple testing.
# Interpretation: looks good

# Here are we interested in getting the right background distribution.
# Look at p-values from emptyDrops() to look at assumption that barcodes 
# with low total UMI counts are empty droplets with `test.ambient=TRUE`. 
# Ideally, the distribution should be close to uniform.
# Large peaks near zero indicate that barcodes with total counts below lower 
# are not all ambient in origin. This can be resolved by decreasing lower 
# further to ensure that barcodes corresponding to droplets with very 
# small cells are not used to estimate the ambient profile.
set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(lmmp_mf_sce), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value", main="", col="grey80") 

table(e.out$FDR <= 0.001)

# Once we are satisified, we subset our SingleCellExperiment object 
# to retain only the detected cells. Also NAs are removed at this step. 
lmmp_mf_sce <- lmmp_mf_sce[,which(e.out$FDR <= 0.001)]


##################################
######## Quality Control #########
##################################

is.mito <- as.character(rowRanges(lmmp_mf_sce)@seqnames) == "chrM"
table(is.mito)
# > table(is.mito)
# is.mito
# FALSE  TRUE 
# 54310    37 

# addPerCellQC() adds qc metrics directly to colData() slot
lmmp_mf_sce <- scater::addPerCellQC(lmmp_mf_sce, subsets=list(Mito=is.mito))

# perCellQCMetrics() calculates qc metrics as stand alone table
stats <- scater::perCellQCMetrics(lmmp_mf_sce, subsets=list(Mito=which(is.mito)))

# quickPerCellQC() uses adaptive thresholds to identify outliers and 
#   does not work well here b/c it removes way too many from % mito. 
#   So we are going to investigate by hand
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch = lmmp_mf_sce$sample_id)
colSums(as.matrix(qc))
# > colSums(as.matrix(qc))
# low_lib_size            low_n_features high_subsets_Mito_percent 
# 32                       238                       494 
# discard 
# 744 


# Plot num of detected genes (x-axis) vs % mito (y-axis)
plotColData(lmmp_mf_sce, x="detected", y="subsets_Mito_percent", 
            colour_by = "sample_id") + 
  ggtitle("Detected features")
# Looks like TL sample has some higher % mito cells
# Subhash asked us to be liberal in what we accept using % mito so we 
#   will consider two thresholds (both very liberal)

###### Set qc threshold for % mito
# Below are two thresholds for % mito: one adaptive and one fixed (based on plot above)
colData(lmmp_mf_sce)$discard.mito.adapt <- 
  isOutlier(lmmp_mf_sce$subsets_Mito_percent, type ="higher", 
            log=FALSE, nmads = 10, batch=lmmp_mf_sce$sample_id)
attr(lmmp_mf_sce$discard.mito.adapt, "thresholds")

colData(lmmp_mf_sce)$discard.mito.fixed <- lmmp_mf_sce$subsets_Mito_percent > 30
table(lmmp_mf_sce$sample_id, lmmp_mf_sce$discard.mito.fixed)


###### Set qc threshold for num detected features
colData(lmmp_mf_sce)$discard.features <- 
  isOutlier(lmmp_mf_sce$detected, type = "lower", 
            log = TRUE, nmads = 4,
            batch = lmmp_mf_sce$sample_id)
attr(lmmp_mf_sce$discard.features, "thresholds")

###### Set qc threshold for 
colData(lmmp_mf_sce)$discard.sum <- 
  isOutlier(lmmp_mf_sce$sum, type = "lower", log = TRUE, nmads = 3,
                         batch = lmmp_mf_sce$sample_id)
attr(lmmp_mf_sce$discard.sum, "thresholds")

###### Plot results
gridExtra::grid.arrange(
  plotColData(lmmp_mf_sce, x="sample_id", y="sum", colour_by="discard.sum") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(lmmp_mf_sce, x="sample_id", y="detected", colour_by = "discard.features") 
  + scale_y_log10() + ggtitle("Detected features"),
  nrow=1, ncol=2
)

# plot by sample_id
plotColData(lmmp_mf_sce, x="detected", y="subsets_Mito_percent", 
            colour_by = "sample_id") + 
  ggtitle("Detected features")


gridExtra::grid.arrange(
  plotColData(lmmp_mf_sce, x="sample_id", y="subsets_Mito_percent", colour_by = "discard.mito.adapt") 
  + ggtitle("Mito percent"),
  plotColData(lmmp_mf_sce, x="sample_id", y="subsets_Mito_percent", colour_by = "discard.mito.fixed") 
  + ggtitle("Mito percent"),
  plotColData(lmmp_mf_sce, x="detected", y="subsets_Mito_percent", 
              colour_by = "discard.mito.adapt") + 
    ggtitle("Detected features"),
  plotColData(lmmp_mf_sce, x="detected", y="subsets_Mito_percent", 
              colour_by = "discard.mito.fixed") + 
    ggtitle("Detected features"),
  nrow=2, ncol=2
)

# Final qc discard column
lmmp_mf_sce$discard <- lmmp_mf_sce$discard.sum | 
                       lmmp_mf_sce$discard.features | 
                       lmmp_mf_sce$discard.mito.fixed

data.frame("total" = sum(lmmp_mf_sce$discard.sum), 
           "num_features" = sum(lmmp_mf_sce$discard.features),
           "perc_mito" = sum(lmmp_mf_sce$discard.mito.fixed), 
           "final_discard" = sum(lmmp_mf_sce$discard))

# This seems reasonable and very liberal in what we put through qc;
#   Might have to readjust later, but for now we assume it's OK

lmmp_mf_sce <- lmmp_mf_sce[,!lmmp_mf_sce$discard]


##################################
######### Normalization ##########
##################################

set.seed(1000)
clusters <- scran::quickCluster(lmmp_mf_sce)
table(clusters)
lmmp_mf_sce <- computeSumFactors(lmmp_mf_sce, cluster = clusters)
lmmp_mf_sce <- logNormCounts(lmmp_mf_sce)
summary(sizeFactors(lmmp_mf_sce))

data.frame(colData(lmmp_mf_sce), libsizfact = librarySizeFactors(lmmp_mf_sce)) %>% 
  ggplot(aes(x=libsizfact, y=sizeFactor, col = subsets_Mito_percent)) + 
  geom_point() + scale_x_log10() + scale_y_log10() + 
  xlab("Library size factors") + ylab("Deconvolution factors")

# Interpretation: We see that the deconvolution size factors exhibit 
#   cell type-specific deviations from the library size factors. 
#   This is consistent with the presence of composition biases that are 
#   introduced by strong differential expression between cell types. 
#   Use of the deconvolution size factors adjusts for these biases to improve 
#   normalization accuracy for downstream applications.

saveRDS(lmmp_mf_sce, file = here("data", "sce_norm.rds"))
