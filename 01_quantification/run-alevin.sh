#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=5G

d=/fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages

# set project directory path (TH sample)
cd $d/fastq_files/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY

salmon alevin -l ISR \
    --index $d/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-nodecoys \
    -1 SK-TH_S1_L001_R1_001.fastq.gz SK-TH_S1_L001_R1_002.fastq.gz SK-TH_S1_L001_R1_003.fastq.gz SK-TH_S1_L001_R1_004.fastq.gz SK-TH_S1_L002_R1_001.fastq.gz SK-TH_S1_L002_R1_002.fastq.gz SK-TH_S1_L002_R1_003.fastq.gz SK-TH_S1_L002_R1_004.fastq.gz \
    -2 SK-TH_S1_L001_R2_001.fastq.gz SK-TH_S1_L001_R2_002.fastq.gz SK-TH_S1_L001_R2_003.fastq.gz SK-TH_S1_L001_R2_004.fastq.gz SK-TH_S1_L002_R2_001.fastq.gz SK-TH_S1_L002_R2_002.fastq.gz SK-TH_S1_L002_R2_003.fastq.gz SK-TH_S1_L002_R2_004.fastq.gz \
    --tgMap $d/salmon_index_files/gencode.vM25.transcripts.tx2gene.mouse.tsv \
    --chromium \
    --threads 6 \
    --output $d/01_quantification/salmon_quants/sample_TH_quant \
    --dumpFeatures --dumpBfh


# set project directory path (TL sample)
cd $d/fastq_files/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY

salmon alevin -l ISR \
    --index $d/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-nodecoys \
    -1 SK-TL_S1_L001_R1_001.fastq.gz SK-TL_S1_L001_R1_002.fastq.gz SK-TL_S1_L001_R1_003.fastq.gz SK-TL_S1_L002_R1_001.fastq.gz SK-TL_S1_L002_R1_002.fastq.gz SK-TL_S1_L002_R1_003.fastq.gz \
    -2 SK-TL_S1_L001_R2_001.fastq.gz SK-TL_S1_L001_R2_002.fastq.gz SK-TL_S1_L001_R2_003.fastq.gz SK-TL_S1_L002_R2_001.fastq.gz SK-TL_S1_L002_R2_002.fastq.gz SK-TL_S1_L002_R2_003.fastq.gz \
    --tgMap $d/salmon_index_files/gencode.vM25.transcripts.tx2gene.mouse.tsv \
    --chromium \
    --threads 10 \
    --output $d/01_quantification/salmon_quants/sample_TL_quant \
    --dumpFeatures --dumpBfh
