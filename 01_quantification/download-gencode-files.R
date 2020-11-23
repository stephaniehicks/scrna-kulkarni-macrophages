# download-gencode-files.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Nov 18, 2020
#
# Download GENCODE files and generate files for the pipelines


library(here)

# Create directories
if(!dir.exists(here("salmon_index_files"))){
  dir.create(here("salmon_index_files"))
}
if(!dir.exists(here("01_quantification", "salmon_quants"))){
  dir.create(here("01_quantification", "salmon_quants"))
}

# Download GENCODE Files
# download GENCODE primary assembly fasta file
if(!file.exists(here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
  download.file(tar_gz_file, 
                destfile = here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz"), 
                method = "wget")
}

# download GENCODE transcripts fasta file
if(!file.exists(here("salmon_index_files", "gencode.vM25.transcripts.fa.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz"
  download.file(tar_gz_file,
                destfile = here("salmon_index_files", "gencode.vM25.transcripts.fa.gz"),
                method = "wget")
}

# download GENCODE gtf file
if(!file.exists(here("salmon_index_files", "gencode.vM25.annotation.gtf.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
  download.file(tar_gz_file, 
                destfile = here("salmon_index_files", "gencode.vM25.annotation.gtf.gz"), 
                method = "wget")
}


# Make files for pipelines (code adapted from https://github.com/csoneson/rna_velocity_quant)
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(Biostrings)
  library(rtracklayer)
  library(GenomicFeatures)
  library(BSgenome)
})
source(here("01_quantification", "quantify-salmon-helpers.R")) 


#######################
# transcripts pipline #
#######################
## FASTA file
# Gtf path
gtf_file <- here("salmon_index_files", "gencode.vM25.annotation.gtf.gz")

# Read genomic (DNA) sequence from FASTA file
genome_fasta <- here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz") 
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1) # creates chr1, etc

# Extract transcript (tx) sequences (takes a few minutes)
tx <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "spliced")

# Write FASTA file
Biostrings::writeXStringSet(tx, file = here("salmon_index_files", "gencode.vM25.transcripts.mouse.fa.gz"), 
                            compress = T) # Compresses fa while saving, rather than compress later (takes a few min)

## tx2gene
# Read gtf
gtf <- rtracklayer::import(here("salmon_index_files",
                                "gencode.vM25.annotation.gtf.gz")) 
gtftx <- subset(gtf, type == "transcript")
gtfex <- subset(gtf, type == "exon")

df <- data.frame(gtftx, stringsAsFactors = FALSE) %>%
  dplyr::select(transcript_id, seqnames, start, end, strand, source, 
                gene_id, gene_type, gene_name, level, havana_gene, transcript_type,
                transcript_name, transcript_support_level, tag, havana_transcript) %>%
  dplyr::left_join(data.frame(gtfex, stringsAsFactors = FALSE) %>%
                     dplyr::group_by(transcript_id) %>% 
                     dplyr::summarize(transcript_length = sum(width)),
                   by = "transcript_id")

# Write table as txt, tsv, and rds
write.table(df %>% dplyr::select(transcript_id, gene_id), 
            file = here("salmon_index_files", "gencode.vM25.transcripts.tx2gene.mouse.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)
readr::write_tsv(x = df %>% dplyr::select(transcript_id, gene_id), 
            file = here("salmon_index_files", "gencode.vM25.transcripts.tx2gene.mouse.tsv"),
            col_names = FALSE)
saveRDS(df, file = here("salmon_index_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))

