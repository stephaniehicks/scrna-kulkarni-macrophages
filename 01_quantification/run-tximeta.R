# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo, Stephanie Hicks
# Date last modified: Nov 23, 2020
#
# Use tximeta to read in quant files (output from alevin) into SingleCellExperiment files

suppressPackageStartupMessages({
  library(here)
  library(tximeta)
  library(fishpond)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(org.Mm.eg.db) # org package for mouse
})

# Import with tximeta
# Note: alevin import currently only supports a single experiment/sample at a time; we have two samples here
sample_names = c("TH", "TL")
file_paths = here("01_quantification", "salmon_quants",  
                  paste0("sample_", sample_names, "_quant"), 
                  "alevin", "quants_mat.gz")
                 
for(i in seq_len(length(sample_names))){
    print(file_paths[i])
    se = tximeta(coldata = data.frame(names = sample_names[i],
                                      files = file_paths[i],
                                      stringsAsFactors = FALSE),
                 type = "alevin")
    
    # Check SummarizedExperiment object
    colData(se)
    assayNames(se)
    rowRanges(se)
    
    # Convert to a SingleCellExperiment (sce) object
    sce <- as(se, "SingleCellExperiment")
    
    # add gene IDs to sce object
    sce <- addIds(sce, "SYMBOL")
    
    # Save as a sce object
    saveRDS(sce, file = here("01_quantification", "salmon_quants", 
                             paste0("sce_", sample_names[i], ".rds")))
}

# Next, we load in the two `SingleCellExperiment` (`sce`) objects and combining them into one `sce` object. 
file_paths = here("01_quantification", "salmon_quants",  
                  paste0("sce_", sample_names[i], ".rds"))

sce_ls = vector("list", length(sample_names))
for(i in seq_len(sample_names)){
  sce_ls[[sample_names[i]]] <- readRDS(file = file_paths[i])  
}

# !! Still need to add code to combine list of objects together !!


# save clean data here
if(!file.exists(here("data"))){
  dir.create(here("data"))
}
saveRDS(sce_combined, file = here("data", "sce_macrophages.rds"))

