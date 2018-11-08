# NORAD manuscript

## MYCN knock-down and overlapping gene signature with NORAD knock-down

Raw sequencing data from MYCN knock-down in SK-N-BE(2) (ArrayExpress accession [E-GEOD-84389](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-84389/, SRA [GSE84389](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84389)) was mapped to GRCh38 coordinates:

```bash
#!/bin/bash

#SBATCH -N1 -n16 --mem-per-cpu=2000 -t4:00:00
#SBATCH --array=1-9
#SBATCH -e hisat-%A_%a.err
#SBATCH -o hisat-%A_%a.out

i=$SLURM_ARRAY_TASK_ID
index=HISAT2_INDEX/grch38_snp_tran/genome_snp_tran
sample=$(eval 'sed "${i}q;d" SRA.txt') # SRA.txt contains list of SRA accessions (one per line)
sam_outdir=PRJNA329050_GSE84389/sam/
mkdir -p $sam_outdir
sam=$sam_outdir$sample\.sam
hisat2 -p $SLURM_TASKS_PER_NODE -x $index --sra-acc $sample -S $sam
```

Counts were obtained with featureCounts:

```bash
#!/bin/bash

#SBATCH -d PRJNA329050_GSE84389/sam/
#SBATCH -N1 -n4 --mem-per-cpu=8000 -t4:00:00
#SBATCH -e featureCounts.err
#SBATCH -o featureCounts.out

gtffile=HISAT2_INDEX/grch38_snp_tran/Homo_sapiens.GRCh38.94.gtf
results=all.gene.counts.txt

featureCounts -T 4 -g gene_id  -a $gtffile -o $results SRR3922065.sam SRR3922066.sam SRR3922067.sam SRR3922068.sam SRR3922069.sam SRR3922070.sam SRR3922071.sam SRR3922072.sam SRR3922073.sam
```

Code describing gene set enrichment, differential expression and overlapping gene signature is available in the script [MYCN.NORAD.knockdown.BE2.R](https://github.com/utnesp/NORAD/blob/master/MYCN.NORAD.knockdown.BE2.R). 
Differential expression from NORAD knock-down was based on data in [NORAD.knockdown.CPM.BE2.txt](https://github.com/utnesp/NORAD/blob/master/NORAD.knockdown.CPM.BE2.txt) from the neuroblastoma cell line SK-N-BE(2)c. 
Differential expression results was later compared to hg19 differential expression results obtained from the excel sheet "TEAD4-MYCN KD DGE signature" within [Table S7 RNA-seq derived signatures and pathway analysis from shTEAD4, shMYCN ans shWWTR1](http://cancerdiscovery.aacrjournals.org/highwire/filestream/43006/field_highwire_adjunct_files/6/169577_3_supp_4574645_p49g2y.xlsx) from the article [Cross-Cohort Analysis Identifies a TEAD4–MYCN Positive Feedback Loop as the Core Regulatory Element of High-Risk Neuroblastoma](http://cancerdiscovery.aacrjournals.org/content/8/5/582). Pearson's product-moment correlation comparing hg19 adjusted p-values to hg38 adjusted p-values resulted in r=0.96 (p<2.2e-16, 95% CI 0.959, 0.962).


## ENCODE project and identication of proteins, transcription factors and histones binding to NORAD
To identify proteins binding to the NORAD lncRNA, human RNA-binding data was downloaded as bed files from [ENCODE](https://www.encodeproject.org/search/?type=Experiment&assay_slims=RNA+binding&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens) with the command:

```bash
files=files.txt # file downloaded from the batch download command. Contains URL path to file.
```
An accession files was created from this file using sed with the command: 
```bash
accessions_files=accessions_files.txt

#create accession file
sed 's/.*\///' $files  | sed 's/.bed.gz//g' > $accessions_files
```


```{bash}
outfile=out.bed

IFS=$'\n'       # make newline as separator
for j in $(cat $accessions_files)
do
                wget https://www.encodeproject.org/files/$j\/@@download/$j\.bed.gz # get bed file with enriched peak (format for BED6 is: chr start stop name score strand)
 
                gzip -d $j\.bed.gz # deflate gz file
 
                # get only peaks mapping within the NORAD gene
                awk 'BEGIN{OFS="\t"} $1=="chr20" && $2>=36045622 && $3<=36150960 {print $1,$2, $3, $4="GRCh38_"FILENAME, $5, $6, $7, $8, $9, $10}' $j\.bed >> $outfile

                rm $j\.bed # remove each file that has been downloaded to avoid cramping up space on disk
done
```

The resulting outfile contains a [ENCODE narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7)  (Narrow or Point-Source Peaks) format which is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format with the following columns:
    1. chrom - Name of the chromosome (or contig, scaffold, etc.).
    2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined aschromStart=0, chromEnd=100, and span the bases numbered 0-99.
    4. name - Name given to a region (preferably unique). Use '.' if no name is assigned.
    5. score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
    6. strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
    7. signalValue - Measurement of overall (usually, average) enrichment for the region.
    8. pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
    9. qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
    10. peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
    

This file contains peaks within the NORAD gene which then was filtered using the command:
 ```bash
 awk '$1=="chr20" && $2>=36045622 && $3<=36050960 && $6=="-" && $8>=10 && $7>=2' $outfile >> NORAD.peaks.bed
 ```
 
All accessions within NORAD.peaks.bed and URL paths from ENCODE was then saved to a new file called bigwig.files.txt. Using bigWigToBedGraph, regions within the NORAD gene was downloaded with the following command:
 
```bash 
IFS=$'\n'       # make newlines the only separator
for j in $(cat ./bigwig.files.txt)
do
    bigWigToBedGraph -chrom=chr20 -start=36045622 -end=36050960 -udcDir=/cache $j\.bigWig bedGraph/$j\.bedGraph #  -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
done
```

The resulting bedGraph files was then manually inspected with the plotBedgraph function from the R package [Sushi](https://bioconductor.org/packages/release/bioc/html/Sushi.html).

The same pipeline was followed for NORAD hg19 coordinates (chrom=chr20, -start=34633544 -end=34638882). 
For ChIP-seq data, the same pipeline was also followed. However, a 500 bp upstream of NORAD transcriptional start site (TSS) was included as well as 2000 bp downstream of NORAD TSS.


## Survival analysis of lncRNAs in SEQC dataset
Survial analysis of lncRNAs follows the code supplied in [SEQC.survival.R](https://github.com/utnesp/NORAD/blob/master/SEQC.survival.R) and the function plot.surv() available in [plot.surv.R](https://github.com/utnesp/NORAD/blob/master/plot.surv.R).

## Artwork presentation
[Inkscape](https://inkscape.org/) was used in the development of the artwork associated with the manuscript.


