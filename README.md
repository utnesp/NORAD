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




## ENCODE project and identificaiton of proteins binding to NORAD
To identify proteins binding to NORAD, human RNA-binding data was downloaded as bigWig files from [ENCODE](https://www.encodeproject.org/search/?type=Experiment&assay_slims=RNA+binding&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bigWig). 




## Artwork presentation
[Inkscape](https://inkscape.org/) was used in the development of the artwork associated with the manuscript.
