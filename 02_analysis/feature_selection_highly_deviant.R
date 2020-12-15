# Feature selection using highly deviant genes (input must be un-normalized sce)


# R initiation
suppressPackageStartupMessages({
  library(here)
  library(scater)
  library(SingleCellExperiment)
  library(scry)
})

# load SingleCellExperiment object post QC
sce <- readRDS(here("data", "sce_norm.rds"))
sce

# Removing lowly expressed genes, after QC
num_reads <- 1
num_cells <- 0.01*ncol(sce)
keep <- which(DelayedArray::rowSums(counts(sce) >= num_reads ) >= num_cells)
sce <- sce[keep,]
sce

# Calculate the deviance residuals under a Poisson model 
sce <- scry::nullResiduals(sce, assay="counts", fam = "poisson", 
                           type="deviance", batch = "sample_id")

# Apply PCA to the results 
sce <- scater::runPCA(sce, ncomponents = 50,
                      ntop = 1000,
                      exprs_values = "poisson_deviance_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())
plotReducedDim(sce, dimred = "GLM-PCA", colour_by = "sample_id")
