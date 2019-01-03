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
To identify proteins binding to the NORAD lncRNA, human RNA-binding data was downloaded as bed files from [ENCODE]. 

```R
metadata <- fread("https://www.encodeproject.org/metadata/type=Experiment&status=released&assay_slims=RNA+binding&assay_title=eCLIP/metadata.tsv")
metadata <- metadata[metadata$`File format` == "bed narrowPeak", ]

   
    
cores=detectCores(); cl <- makeCluster(cores); registerDoParallel(cl)

temp <- foreach(i = 1:nrow(metadata), .export = "fread") %dopar% {
        t <- fread(paste("wget -nc -O - ", metadata$`File download URL`[i], "| gzip -d | cat"))
        
if(ncol(t) == 10 & nrow(t) > 0) {
            colnames(t)[1:10] <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
            if(metadata$Assembly[i] == "hg19") t <- t[chromStart >= NORAD.coord$hg19.start & chromEnd <= NORAD.coord$hg19.end]
            if(metadata$Assembly[i] == "GRCh38") t <- t[chromStart >= NORAD.coord$hg38.start & chromEnd <= NORAD.coord$hg38.end]
            
                if(nrow(t) > 0) {
                    t$Assembly <- metadata$Assembly[i]
                    t$CellLine <- metadata$`Biosample term name`[i]
                    t$metadata.index <- i
                }
        }
        return(t)
    }
    stopImplicitCluster(); stopCluster(cl)
    NORAD.ENCODE <- rbindlist(temp, fill = T)
```

The [ENCODE narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7)  (Narrow or Point-Source Peaks) format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format with the following columns:
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
    

For visual inspection, all accessions within the NORAD gene was inspected using bigWigToBedGraph -chrom=chr20 -start=36045622 -end=36050960. The resulting bedGraph files was then manually inspected with the plotBedgraph function from the R package [Sushi](https://bioconductor.org/packages/release/bioc/html/Sushi.html). The same pipeline was followed for NORAD hg19 coordinates (chrom=chr20, -start=34633544 -end=34638882). 

For ChIP-seq data, the same pipeline was also followed. However, a 500 bp upstream of NORAD transcriptional start site (TSS) was included as well as 2000 bp downstream of NORAD TSS.

## Survival analysis of lncRNAs in SEQC dataset
Survial analysis of lncRNAs follows the code supplied in [SEQC.survival.R](https://github.com/utnesp/NORAD/blob/master/SEQC.survival.R) and the function plot.surv() available in [plot.surv.R](https://github.com/utnesp/NORAD/blob/master/plot.surv.R).

## Profiling 3' UTR in SK-N-BE(2)c and HCT116 after NORAD knock-down
Profiling of [effective length](https://salmon.readthedocs.io/en/latest/file_formats.html#quantification-file) within 3' UTRs was done using [Salmon](https://combine-lab.github.io/salmon/). Effective Length is is the computed effective length of the target transcript. It takes into account all factors being modeled that will effect the probability of sampling fragments from this transcript, including the fragment length distribution and sequence-specific and gc-fragment bias. Using salmon, an index was created using sequences of all transcripts as well as 5' and 3' UTRs from the fasta file generate from [getUTR.sequence.R](https://github.com/utnesp/NORAD/blob/master/getUTR.sequence.R). Raw fastq files from both SK-N-BE(2)c and HCT116 after NORAD knock-down was thereafter run using the script [salmon.sh](https://github.com/utnesp/NORAD/blob/master/salmon.sh). The effective length change was thereafter computed as the effective length in NORAD minus the effective length in control treated cell.

## Artwork presentation
[Inkscape](https://inkscape.org/) was used in the development of the artwork associated with the manuscript.

