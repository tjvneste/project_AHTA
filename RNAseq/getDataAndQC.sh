#/bin/bash

#Download data
mkdir raw
mkdir QC
cd raw

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/SRR201983/SRR201983.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/SRR201984/SRR201984.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/SRR201985/SRR201985.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/SRR201986/SRR201986.fastq.gz

fastqc -o ../QC -threads 4 *.gz
#Check html files, overrepresented sequences are present
#	Copy paste the overrepresented into a text file in bash (nano adapters.txt, copy from html and right mouse click in nano)
#	Create fasta adapter secuences from this file
cat adapters.txt | cut -f 1 > temp.txt
#find an elegant way to add >adapter \n before each sequence

#Trimming
⁫java -jar /data/student_homes/public/PracticalSession3/trimmomatic-0.36.jar PE -threads 3 breastCancer_1.fastq.gz breastCancer_2.fastq.gz \
trimmed_breastCancer_1.fastq.gz /dev/null trimmed_breastCancer_2.fastq.gz /dev/null ILLUMINACLIP:adapters.fa:2:30:10

java -jar /data/student_homes/public/PracticalSession3/trimmomatic-0.36.jar PE -threads 3 control_1.fastq.gz control_2.fastq.gz \
trimmed_control_1.fastq.gz /dev/null trimmed_control_2.fastq.gz /dev/null ILLUMINACLIP:adapters.fa:2:30:10

#QC to check if trimming was succesfull
fastqc -o ../QC -threads 2 trimmed*.gz
