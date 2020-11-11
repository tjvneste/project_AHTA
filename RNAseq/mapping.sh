#!/bin/bash 

mkdir genome
cd genome

#download hg18 genome and annotation
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip *.gz

gtf_to_fasta Homo_sapiens.GRCh38.101.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa hg18_transcriptome.fa
kallisto index -i hg18_index hg18_transcriptome.fa
cd ..
mkdir counts
cd raw

for file in $(ls trimmed*gz | sed -e "s/trimmed_//" | sed "s/_.*//" | sort -u ); { kallisto quant -i ../genome/hg18_index \
-o ../counts/$file --threads 3 "trimmed_"$file"_1.fastq.gz" "trimmed_"$file"_2.fastq.gz"; }

# Paste quantification files together to create a count file
cd ../counts
paste breastCancer/abundance.tsv control/abundance.tsv | cut -f 1,2,4,5,9,10 > counts.txt
# Now download this file to your local computer for further analysis in R
