#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G,h_fsize=30G


# create salmon index with decoys (this process takes ~1 hour for mouse)
salmon index -t /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files/gentrome_transcripts_mouse.fa.gz \
             -d /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files/decoys_mouse.txt \
             -i /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-nodecoys \
             --gencode --threads 4 -k 31

# # create salmon index without decoys (this process takes <30 mins for mouse)
# salmon index -t /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files/gentrome_transcripts_mouse.fa.gz \
#              -i /fastscratch/myscratch/shicks1/scrna-kulkarni-macrophages/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-nodecoys \
#              --gencode --threads 4 -k 31
