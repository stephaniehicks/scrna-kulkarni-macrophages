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

# create linkedTranscriptome for decoys pipeline
index_dir = here("salmon_index_files", "gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys")
fasta_path = here("salmon_index_files", "gencode.vM25.transcripts.mouse.fa.gz")
gtf_path = here("salmon_index_files", "gencode.vM25.annotation.gtf.gz")
json_file = here("salmon_index_files", paste0(basename(index_dir), ".json"))
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M25", genome="GRCm38", 
                fasta=fasta_path, gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

# Import with tximeta
# Note: alevin import currently only supports a single experiment/sample at a time; we have two samples here
sample_names = c("TH", "TL")
file_paths = here("01_quantification", "salmon_quants",  
                  paste0("sample_", sample_names, "_quant"), 
                  "alevin", "quants_mat.gz")
   
counts_mat = variance_mat = mean_mat = colnames_mat = coldata_mat <- NULL              
for(i in seq_len(length(sample_names))){
    print(file_paths[i])
    se = tximeta(coldata = data.frame(names = sample_names[i],
                                      files = file_paths[i],
                                      stringsAsFactors = FALSE),
                 type = "alevin", 
                 dropInfReps=TRUE)#, 
                 # alevinArgs=list(filterBarcodes=TRUE))

    colData(se)$sample_id <- rep(sample_names[i], ncol(se))

    # Save as a SummarizedExperiment object
    saveRDS(se, file = here("01_quantification", "salmon_quants", 
                             paste0("se_", sample_names[i], ".rds")))
    
    counts_mat <- cbind(counts_mat, assay(se, "counts"))
    variance_mat <- cbind(variance_mat, assay(se, "variance"))
    mean_mat <- cbind(mean_mat, assay(se, "mean"))
    colnames_mat <- c(colnames_mat, colnames(se))
    coldata_mat <- rbind(coldata_mat, colData(se))
}

# Create new SE
se_new = SummarizedExperiment(assays = SimpleList(counts = counts_mat, 
                                                  variance = variance_mat, 
                                                  mean = mean_mat),
                              colData = coldata_mat,
                              rowRanges = rowRanges(se),
                              metadata = metadata(se))

# Convert to a SingleCellExperiment (sce) object
sce <- as(se_new, "SingleCellExperiment")
    
# add gene IDs to sce object
sce <- addIds(sce, "SYMBOL")
mcols(sce)

# Check sce object
colData(sce)
assayNames(sce)
rowRanges(sce)

# save combined sce data
if(!file.exists(here("data"))){
  dir.create(here("data"))
}
saveRDS(sce, file = here("data", "sce_mouse_TH_TL.rds"))

