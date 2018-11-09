#!/bin/bash

#SBATCH -D /path/to/fastq/directory/
#SBATCH -N1 -n16 --mem-per-cpu=2000 -t4:00:00
#SBATCH --array=1-6
#SBATCH -e salmon-%A_%a.err
#SBATCH -o salmon-%A_%a.out

i=$SLURM_ARRAY_TASK_ID

# index was created using the command: 
# salmon index -t GRCh38.utr.fa.gz -i Homo_sapiens.GRCh38.ensembl.index
# GRCh38.utr.fa.gz was created from https://github.com/utnesp/NORAD/blob/master/getUTR.sequence.R
index_path=/global/work/put001/data/GRCh38.ensembl/
index=$index_path\Homo_sapiens.GRCh38.ensembl.index

sample=$(eval 'sed "${i}q;d" salmon.fastq.files.to.count.txt') #txt file contains basename of each fastq file (one per line)
# define path to raw fastq files 
fastq_path=/path/to/fastq/directory/
read1=$fastq_path$sample\_R1.fastq.gz
read2=$fastq_path$sample\_R2.fastq.gz

results=./salmon_results/combined/$sample\_quant
mkdir -p $results
salmon quant -i $index -l A -1 $read1 -2 $read2 -p 16 -o $results --numBootstraps 100 --seqBias --gcBias
