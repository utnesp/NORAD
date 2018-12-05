#!/bin/bash

#PBS -d /global/work/put001/SRA/sra
#PBS -lnodes=1:ppn=1,pmem=32gb,walltime=24:00:00
#PBS -t 1-161
#PBS -e /global/work/put001/SRA/sra/log
#PBS -o /global/work/put001/SRA/sra/log

i=$PBS_ARRAYID

printf "Downloading SRA files\n"
echo "Downloading SRA files\n"
name=$(eval 'sed "${i}q;d" /home/put001/scripts2/SRA/SRA.txt')
prefetch $name

printf "Dump fastq\n"
echo "Dump fastq\n"
fastq-dump --outdir /global/work/put001/SRA/sra/fastq -I --split-files $name
cd fastq
printf "Fastqc on raw fastq\n"
fastqc --outdir fastqc $name\_1.fastq $name\_2.fastq

printf "Trim adapters\n"
echo "Trim adapters\n"
output_trim=/global/work/put001/SRA/sra/fastq/trimmed
trim_galore --paired --quality 28 --stringency 6 --output_dir $output_trim --dont_gzip $name\_1.GGG
# The "singles" file contains reads that passed filter in either the forward or reverse direction, but not the other
printf "Trim low quality bases, reads and reads containing n's\n"
sickle pe --no-fiveprime --qual-type sanger --qual-threshold 28 --length-threshold 20 --truncate-n --pe-file1 $name\_1_val_1.fq --pe-file2 $name\_2_val_2.fq --output-pe1 sickle_trimmed/$name\_1.fastq --output-pe2 sickle_trimmed/$name\_2.fastq --output-single sickle_trimmed/singles/$name\_singles.fastq

# Declare variables
bam_output_path=/global/work/put001/SRA/sra/bam
read1=/global/work/put001/SRA/sra/fastq/trimmed/sickle_trimmed/$name\_1.fastq
read2=/global/work/put001/SRA/sra/fastq/trimmed/sickle_trimmed/$name\_2.fastq
genfilework=/home/put001/data/GRCh38.p3/ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.reordered.fa
genomeloc=/global/work/put001/data/GRCh38.p3/ENSEMBL/STAR_index_GRCh38.p3_non_chr/
gtffile=/home/put001/data/GRCh38.p3/ENSEMBL/Homo_sapiens.GRCh38.82.gtf

printf "Do fastqc before starting mapping\n"
echo "Do fastqc before starting mapping\n"
fastqc --outdir /global/work/put001/SRA/sra/fastq/trimmed/sickle_trimmed/fastqc $read1 $read2

printf "Perform mapping\n"
echo "Perform mapping\n"
cd $bam_output_path
STAR_2.5.1b --clip5pNbases 13 --clip3pNbases 4 --outSAMunmapped Within --outBAMcompression 10 --outFileNamePrefix $name --genomeDir $genomeloc --outSAMtype BAM SortedByCoordinate --readFilesIn $read1 $read2 --chimSegmentMin 20 --runThreadN 1

printf "Fastqc: Check that the BAM file looks OK\n"
echo "Fastqc: Check that the BAM file looks OK\n"
fastqc --outdir /global/work/put001/SRA/sra/bam/fastqc $name\Aligned.sortedByCoord.out.bam

printf "Index the sorted bam files\n"
echo "Index the sorted bam files\n"
samtools index $name\Aligned.sortedByCoord.out.bam
