#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=2G,h_vmem=2G

cd /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files

# 1. Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys_mouse.txt
sed -i.bak -e 's/>//g' decoys_mouse.txt

# 2. Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat gencode.vM25.transcripts.mouse.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome_transcripts_mouse.fa.gz
