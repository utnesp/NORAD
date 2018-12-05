#!/bin/bash

#PBS -d /global/work/put001/SRA/sra/bam
#PBS -lnodes=1:ppn=16,pmem=8gb,walltime=99:99:99
#PBS -e /global/work/put001/SRA/sra/log/featurecounts
#PBS -o /global/work/put001/SRA/sra/log/featurecounts
#PBS -q highmem

genfilework=/home/put001/data/GRCh38.p3/ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.reordered.fa
gtffile=/home/put001/data/GRCh38.p3/ENSEMBL/Homo_sapiens.GRCh38.82.gtf

bamFiles=`awk 'BEGIN { OFS = ""} {print $1,".bam"}' /home/put001/scripts2/SRA/SRA.txt | paste -sd " " -`
featureCounts -T 16 -g gene_id -J -G $genfilework -p -s 1 -d 20 -a $gtffile -o all.gene.counts.txt $(eval 'cat $bamFiles')
featureCounts -T 16 -g exon_id -f -J -G $genfilework -p -s 1 -d 20 -a $gtffile -o all.exon.counts.txt $(eval 'cat $bamFiles')

printf "featureCounts was run on these files:\n"
printf $(eval 'cat $bam_file_to_count') | sed s/\ /\n/g
