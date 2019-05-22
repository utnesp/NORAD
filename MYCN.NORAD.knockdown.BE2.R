library(data.table)
library(ggplot2)
library(ggrepel)
library(enrichR) # https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
library(edgeR)
library(biomaRt)
library(easybiomart)
library(easyedgeR)
library(gridExtra)

mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')

setEPS()
workingDir = "GSE84389/"
dir.create(workingDir, recursive = T)
setwd(workingDir)

# read in RNA-seq data after NORAD knock-down in SK-N-BE(2) (supplementary file 2 in manuscript)
NORAD.kd <- fread("https://raw.githubusercontent.com/utnesp/NORAD/master/NORAD.knockdown.CPM.BE2.txt")
setkey(NORAD.kd, external_gene_name)

url="https://raw.githubusercontent.com/utnesp/NORAD/master/GSE84389.hg38.txt.gz"
EGEOD84389.hg38 <- fread(paste(paste('URL=',url, sep = ""), ' ; curl "$URL" | gunzip -c', sep = ""))
EGEOD84389.hg38 <- EGEOD84389.hg38[, c(1,7:ncol(EGEOD84389.hg38)), with = F]
group = factor(c(1,1,1,2,2,2))
design <- model.matrix(~group)
t <- as.data.frame(EGEOD84389.hg38[, c(2:4,8:10), with = F])
row.names(t) <- EGEOD84389.hg38$Geneid
colnames(t) <- gsub(".sam", "", colnames(t))
names = c(paste("PR", 1:3, "_BE2_control", sep = ""), paste("PR", 7:9, "_BE2_shMYCN", sep = ""))
colnames(t) <- paste(colnames(t), names, sep = "_")

# edgeR differential expression
y <- easyedgeR::DGEList.filter(t, group, cutoff = 5, method = "RLE", filter.fun = "AEG") # removes genes with average expression in each group less than 10. Normalize using RLE (Anders, Huber)
suppressWarnings(easyedgeR::DGE.toptags(y, group = group, design = design, coef = 2, method = "glmQLFit"))
plotMDS(y, gene.selection = "pairwise") 
# limma voom differential expression
v <- voom(y, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
v.topTags <- topTable(fit, coef=2, number = Inf)
# compare edgeR and limma since limma was used in MYCN and TEAD4 paper
plot(-log10(v.topTags$adj.P.Val), -log10(y.topTags$FDR))
cor.test(-log10(v.topTags$adj.P.Val), -log10(y.topTags$FDR))
# continue with edgeR results (since limma and edgeR were effectively the same)
counts <- cpm(y, method = "RLE")
counts <- counts[match(row.names(y.topTags), row.names(counts)), ]
MYCN.kd.hg38 <- cbind(counts, y.topTags)
# round and add annotation
MYCN.kd.hg38$PValue <- signif(MYCN.kd.hg38$PValue, 3)
MYCN.kd.hg38$FDR <- signif(MYCN.kd.hg38$FDR, 3)
MYCN.kd.hg38[,1:6] <- apply(MYCN.kd.hg38[,1:6], 2, function(x) round(x, 0))
MYCN.kd.hg38$logFC <- round(MYCN.kd.hg38$logFC, 2)
MYCN.kd.hg38$F <- round(MYCN.kd.hg38$F, 0)
MYCN.kd.hg38$logCPM <- round(MYCN.kd.hg38$logCPM, 1)
MYCN.kd.hg38$ensembl_gene_id <- row.names(MYCN.kd.hg38)
MYCN.kd.hg38 <- ensg2ext_name_biotype(MYCN.kd.hg38$ensembl_gene_id, combine = T)
MYCN.kd.hg38 <- MYCN.kd.hg38[order(MYCN.kd.hg38$FDR), ]

# genes common for NORAD and MYCN knock-down 
t1 <- MYCN.kd.hg38[, c(1,4:9, 10:14)]
colnames(t1)[8:12] <- paste(colnames(t1)[8:12], "MYCN", sep = "_")
t2 <- NORAD.kd
colnames(t2)[7:11] <- paste(colnames(t2)[7:11], "NORAD", sep = "_") 
diff <- merge(t1, t2, by = "ensembl_gene_id")
diff <- as.data.table(diff)
diff <- diff[FDR_NORAD <= 0.05 | FDR_MYCN <= 0.05, ]
# plot venn
par(mfrow = c(1,2))
input = list(
NORAD.up = unlist(diff[logFC_NORAD > 0 & FDR_NORAD <= 0.05, "external_gene_name"], use.names = F), 
MYCN.up = unlist(diff[logFC_MYCN > 0 & FDR_MYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)
input = list(
NORAD.down = unlist(diff[logFC_NORAD < 0 & FDR_NORAD <= 0.05, "external_gene_name"], use.names = F), 
MYCN.down = unlist(diff[logFC_MYCN < 0 & FDR_MYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)

# use processed data from MYCN knock-down in SK-N-BE(2) (E-GEOD-84389) and correlate hg38 to hg19
EGEOD84389 <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232918_PR1_BE2_control.txt", col.names = c("external_gene_name", "GSM2232918_PR1_BE2_control"), key = "external_gene_name")
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232919_PR2_BE2_control.txt", col.names = c("external_gene_name", "GSM2232919_PR2_BE2_control"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232920_PR3_BE2_control.txt", col.names = c("external_gene_name", "GSM2232920_PR3_BE2_control"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232921_PR4_BE2_shTEAD4.txt", col.names = c("external_gene_name", "GSM2232921_PR4_BE2_shTEAD4"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232922_PR5_BE2_shTEAD4.txt", col.names = c("external_gene_name", "GSM2232922_PR5_BE2_shTEAD4"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232923_PR6_BE2_shTEAD4.txt", col.names = c("external_gene_name", "GSM2232921_PR6_BE2_shTEAD4"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232924_PR7_BE2_shMYCN.txt", col.names = c("external_gene_name", "GSM2232924_PR7_BE2_shMYCN"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232925_PR8_BE2_shMYCN.txt", col.names = c("external_gene_name", "GSM2232925_PR8_BE2_shMYCN"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-84389/E-GEOD-84389.processed.1.zip/GSM2232926_PR9_BE2_shMYCN.txt", col.names = c("external_gene_name", "GSM2232926_PR9_BE2_shMYCN"), key = "external_gene_name")
EGEOD84389 <- merge(EGEOD84389, temp)
temp <- fread("https://raw.githubusercontent.com/utnesp/NORAD/master/PMID_29510988_sheet1.txt") # supplemental table 3, sheet1 from http://cancerdiscovery.aacrjournals.org/content/8/5/582 saved as txt. PMID: 29510988. Rajbhandari P, Lopez G, Capdevila C, Salvatori B et al. Cross-Cohort Analysis Identifies a TEAD4-MYCN Positive Feedback Loop as the Core Regulatory Element of High-Risk Neuroblastoma. Cancer Discov 2018 May;8(5):582-599. 
colnames(temp)[1] <- colnames(EGEOD84389)[1]; temp$V8 <- NULL
EGEOD84389 <- merge(EGEOD84389, temp)
EGEOD84389 <- EGEOD84389[order(EGEOD84389$FDR_shMYCN), ]
EGEOD84389$range_shMYCN <- rowMeans(2^EGEOD84389[, grep("BE2_shMYCN", colnames(EGEOD84389)), with = F]) -
rowMeans(2^EGEOD84389[, grep("BE2_control", colnames(EGEOD84389)), with = F])
EGEOD84389$FC <- rowMeans(2^EGEOD84389[, grep("BE2_shMYCN", colnames(EGEOD84389)), with = F]) /
rowMeans(2^EGEOD84389[, grep("BE2_control", colnames(EGEOD84389)), with = F])
EGEOD84389$log2FC <- log2(EGEOD84389$FC)
setkey(EGEOD84389, external_gene_name)
MYCN.kd <- EGEOD84389[, c(1,2,3,4,8,9,10,14:19), with = F]
colnames(MYCN.kd) <- c("external_gene_name", "PR1_control", "PR2_control", "PR3_control", "PR7_shMYCN", "PR8_shMYCN", "PR9_shMYCN", "t", "p.value", "FDR", "range", "FC", "log2FC")
MYCN.kd$FC <- NULL
# venn diagram showing genes common and specific for TEAD4 and MYCN knock-down 
diff = EGEOD84389
diff <- diff[abs(FDR_shTEAD4) <= 0.05 | abs(FDR_shMYCN) <= 0.05, ]
par(mfrow = c(1,2))
input = list(
TEAD.up = unlist(diff[t_shTEAD4 > 0 & FDR_shTEAD4 <= 0.05, "external_gene_name"], use.names = F), 
MYCN.up = unlist(diff[t_shMYCN > 0 & FDR_shMYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)
input = list(
TEAD.up = unlist(diff[t_shTEAD4 < 0 & FDR_shTEAD4 <= 0.05, "external_gene_name"], use.names = F), 
MYCN.down = unlist(diff[t_shMYCN < 0 & FDR_shMYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)
# correlate hg38 to hg19 annotations
t1 <- ensg2ext_name_biotype(t1$ensembl_gene_id, combine = T)
t3 <- merge(EGEOD84389, t1, by = "external_gene_name")
cor.test(-log10(t3$FDR_MYCN), -log10(t3$FDR_shMYCN)); plot(-log10(t3$FDR_MYCN), -log10(t3$FDR_shMYCN))
cor.test(-log10(t3$FDR_MYCN), -log10(t3$FDR_shTEAD4)); plot(-log10(t3$FDR_MYCN), -log10(t3$FDR_shTEAD4))

sessionInfo()
