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
sce <- readRDS(here("data", "tc", "sce_p10.rds")) 

# make gene names unique
rownames(sce) <- scater::uniquifyFeatureNames(
  rowData(sce)$gene_id, rowData(sce)$SYMBOL)

# Check alevinQC results
#Check file integrity
baseDir <- here("01_quantification", "salmon_quants", "sample_p10_quant")
checkAlevinInputFiles(baseDir = baseDir)

outputDir <- here("Figures", "tc")
alevinQCReport(baseDir = baseDir, sampleId = "p10", 
               outputFile = "alevinReport_p10.pdf", 
               outputFormat = "pdf_document",
               outputDir = outputDir, forceOverwrite = TRUE)

##################################
#### Identify empty droplets #####
##################################

# Examine for empty droplets with DropletUtils
bcrank <- DropletUtils::barcodeRanks(counts(sce))
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
set.seed(10010)
e.out <- DropletUtils::emptyDrops(counts(sce))

summary(e.out$FDR <= 0.001)
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)



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
all.out <- emptyDrops(counts(sce), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value", main="", col="grey80") 

table(e.out$FDR <= 0.001)

# Once we are satisified, we subset our SingleCellExperiment object 
# to retain only the detected cells. Also NAs are removed at this step. 
sce <- sce[,which(e.out$FDR <= 0.001)]


##################################
######## Quality Control #########
##################################

is.mito <- as.character(rowRanges(sce)@seqnames) == "chrM"
table(is.mito)
# > table(is.mito)
# is.mito
# FALSE  TRUE 
# 54310    37 

# addPerCellQC() adds qc metrics directly to colData() slot
sce <- scater::addPerCellQC(sce, subsets=list(Mito=is.mito))

# perCellQCMetrics() calculates qc metrics as stand alone table
stats <- scater::perCellQCMetrics(sce, subsets=list(Mito=which(is.mito)))

# quickPerCellQC() uses adaptive thresholds to identify outliers and 
#   does not work well here b/c it removes way too many from % mito. 
#   So we are going to investigate by hand
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch = sce$sample_id)
colSums(as.matrix(qc))

# Plot num of detected genes (x-axis) vs % mito (y-axis)
plotColData(sce, x="detected", y="subsets_Mito_percent", 
            colour_by = "sample_id") + 
  ggtitle("Detected features")
# Looks like TL sample has some higher % mito cells
# Subhash asked us to be liberal in what we accept using % mito so we 
#   will consider two thresholds (both very liberal)

###### Set qc threshold for % mito
###### Adaptive threshold not expected to work well due to very high mito contents

colData(sce)$discard.mito.fixed <- sce$subsets_Mito_percent > 30
table(sce$sample_id, sce$discard.mito.fixed)

#w/ a 30% fixed threshold, too many cells are being thrown. SK suggests a 60% cutoff to preserve possible biologics. 
#      FALSE TRUE
#p10   318   18
#p11   502  125
#p20   313   35
#p21   669    4
#p60   541  314
#p61   251  114

colData(sce)$discard.mito.fixed <- sce$subsets_Mito_percent > 60
table(sce$sample_id, sce$discard.mito.fixed)

#Label those cells with 30%<n<60% mito contents in coldata
sce$highmito <- sce$subsets_Mito_percent > 30 & sce$subsets_Mito_percent <60

###### Set qc threshold for num detected features
colData(sce)$discard.features <- 
  isOutlier(sce$detected, type = "lower", 
            log = TRUE, nmads = 4,
            batch = sce$sample_id)
attr(sce$discard.features, "thresholds")

###### Set qc threshold for total UMI counts
colData(sce)$discard.sum <- 
  isOutlier(sce$sum, type = "lower", log = TRUE, nmads = 3,
            batch = sce$sample_id)
attr(sce$discard.sum, "thresholds")

###### Plot results
gridExtra::grid.arrange(
  plotColData(sce, x="sample_id", y="sum", colour_by="discard.sum") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="sample_id", y="detected", colour_by = "discard.features") 
  + scale_y_log10() + ggtitle("Detected features"),
  nrow=1, ncol=2
)

# plot by sample_id
plotColData(sce, x="detected", y="subsets_Mito_percent", 
            colour_by = "sample_id") + 
  ggtitle("Detected features")


gridExtra::grid.arrange(
  plotColData(sce, x="sample_id", y="subsets_Mito_percent", colour_by = "discard.mito.adapt") 
  + ggtitle("Mito percent"),
  plotColData(sce, x="sample_id", y="subsets_Mito_percent", colour_by = "discard.mito.fixed") 
  + ggtitle("Mito percent"),
  plotColData(sce, x="detected", y="subsets_Mito_percent", 
              colour_by = "discard.mito.adapt") + 
    ggtitle("Detected features"),
  plotColData(sce, x="detected", y="subsets_Mito_percent", 
              colour_by = "discard.mito.fixed") + 
    ggtitle("Detected features"),
  nrow=2, ncol=2
)

# Final qc discard column
sce$discard <- sce$discard.sum | 
  sce$discard.features | 
  sce$discard.mito.fixed

data.frame("total" = sum(sce$discard.sum), 
           "num_features" = sum(sce$discard.features),
           "perc_mito" = sum(sce$discard.mito.fixed), 
           "final_discard" = sum(sce$discard))

sce <- sce[,!sce$discard]

##################################
######### Normalization ##########
##################################

set.seed(1000)
clusters <- scran::quickCluster(sce)
table(clusters)
sce <- computeSumFactors(sce, cluster = clusters)
sce <- logNormCounts(sce)
summary(sizeFactors(sce))

data.frame(colData(sce), libsizfact = librarySizeFactors(sce)) %>% 
  ggplot(aes(x=libsizfact, y=sizeFactor, col = subsets_Mito_percent)) + 
  geom_point() + scale_x_log10() + scale_y_log10() + 
  xlab("Library size factors") + ylab("Deconvolution factors")

sce10 <- sce
#####################################
######### Integrating data ##########
#####################################

universe <- intersect(rownames(sce10), rownames(sce11), rownames(sce20), rownames(sce21),
                      rownames(sce60), rownames(sce61))
length(universe)




