#salmon index at /media/2/Reference/mm/vM25/salmon_index

#P10
# set project directory path (P10 sample)
bash
cd /media/2/Mm_Timecourse_male_ileum_LMMP/p10
salmon alevin -l ISR \
    --index /media/2/Reference/mm/vM25/salmon_index \
    -1 p10_1_S11_L001_R1_001.fastq.gz  p10_1_S11_L002_R1_001.fastq.gz  p10_2_S20_L001_R1_001.fastq.gz  p10_2_S20_L002_R1_001.fastq.gz p10_3_S9_L001_R1_001.fastq.gz p10_3_S9_L002_R1_001.fastq.gz p10_4_S3_L001_R1_001.fastq.gz p10_4_S3_L002_R1_001.fastq.gz p10_a_S24_L001_R1_001.fastq.gz p10_a_S24_L002_R1_001.fastq.gz p10_b_S7_L001_R1_001.fastq.gz p10_b_S7_L002_R1_001.fastq.gz p10_c_S3_L001_R1_001.fastq.gz p10_c_S3_L002_R1_001.fastq.gz p10_d_S22_L001_R1_001.fastq.gz p10_d_S22_L002_R1_001.fastq.gz \
    -2 p10_1_S11_L001_R2_001.fastq.gz  p10_1_S11_L002_R2_001.fastq.gz  p10_2_S20_L001_R2_001.fastq.gz  p10_2_S20_L002_R2_001.fastq.gz p10_3_S9_L001_R2_001.fastq.gz p10_3_S9_L002_R2_001.fastq.gz p10_4_S3_L001_R2_001.fastq.gz p10_4_S3_L002_R2_001.fastq.gz p10_a_S24_L001_R2_001.fastq.gz p10_a_S24_L002_R2_001.fastq.gz p10_b_S7_L001_R2_001.fastq.gz p10_b_S7_L002_R2_001.fastq.gz p10_c_S3_L001_R2_001.fastq.gz p10_c_S3_L002_R2_001.fastq.gz p10_d_S22_L001_R2_001.fastq.gz p10_d_S22_L002_R2_001.fastq.gz \
    --tgMap /media/2/Reference/mm/vM25/txp2gene.tsv \
    --chromiumV3 \
    --threads 6\
    --output /home/yongbao/scrna-kulkarni-macrophages/01_quantification/salmon_quants/p10_quant \
    --dumpFeatures --dumpBfh \
    --numCellBootstraps 30

    
#P11

salmon alevin -l ISR \
    --index /media/2/Reference/mm/vM25/salmon_index \
    -1 p11_1_S7_L001_R1_001.fastq.gz p11_1_S7_L002_R1_001.fastq.gz p11_2_S22_L001_R1_001.fastq.gz p11_2_S22_L002_R1_001.fastq.gz p11_3_S23_L001_R1_001.fastq.gz p11_3_S23_L002_R1_001.fastq.gz p11_4_S15_L001_R1_001.fastq.gz p11_4_S15_L002_R1_001.fastq.gz p11_a_S23_L001_R1_001.fastq.gz p11_a_S23_L002_R1_001.fastq.gz p11_b_S14_L001_R1_001.fastq.gz p11_b_S14_L002_R1_001.fastq.gz p11_c_S21_L001_R1_001.fastq.gz p11_c_S21_L002_R1_001.fastq.gz p11_d_S13_L001_R1_001.fastq.gz p11_d_S13_L002_R1_001.fastq.gz \
    -2 p11_1_S7_L001_R2_001.fastq.gz p11_1_S7_L002_R2_001.fastq.gz p11_2_S22_L001_R2_001.fastq.gz p11_2_S22_L002_R2_001.fastq.gz p11_3_S23_L001_R2_001.fastq.gz p11_3_S23_L002_R2_001.fastq.gz p11_4_S15_L001_R2_001.fastq.gz p11_4_S15_L002_R2_001.fastq.gz p11_a_S23_L001_R2_001.fastq.gz p11_a_S23_L002_R2_001.fastq.gz p11_b_S14_L001_R2_001.fastq.gz p11_b_S14_L002_R2_001.fastq.gz p11_c_S21_L001_R2_001.fastq.gz p11_c_S21_L002_R2_001.fastq.gz p11_d_S13_L001_R2_001.fastq.gz p11_d_S13_L002_R2_001.fastq.gz \
    --tgMap /media/2/Reference/mm/vM25/txp2gene.tsv \
    --chromiumV3 \
    --threads 6\
    --output /home/yongbao/scrna-kulkarni-macrophages/01_quantification/salmon_quants/p11_quant \
    --dumpFeatures --dumpBfh \
    --numCellBootstraps 30
    
#P20