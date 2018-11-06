library(data.table)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pheatmap)
library(enrichR) # https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
library(stringr)
library(edgeR)
library(biomaRt)
library(easybiomart)
library(RColorBrewer)
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


#EGEOD84389.hg38 <- fread.sshfs("/home/put001/scripts2/181026_EGEOD84389_MYCN/results/all.gene.counts.txt", fill = T, servername = "put001@stallo.uit.no:", force = T, skip = 1, header = T)
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

t1 <- MYCN.kd.hg38[, c(1,4:9, 10:14)]
colnames(t1)[8:12] <- paste(colnames(t1)[8:12], "MYCN", sep = "_")
t2 <- NORAD.kd
colnames(t2)[7:11] <- paste(colnames(t2)[7:11], "NORAD", sep = "_") 
diff <- merge(t1, t2, by = "ensembl_gene_id")
diff <- as.data.table(diff)
diff <- diff[FDR_NORAD <= 0.05 | FDR_MYCN <= 0.05, ]

par(mfrow = c(1,2))
input = list(
NORAD.up = unlist(diff[logFC_NORAD > 0 & FDR_NORAD <= 0.05, "external_gene_name"], use.names = F), 
MYCN.up = unlist(diff[logFC_MYCN > 0 & FDR_MYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)
input = list(
NORAD.down = unlist(diff[logFC_NORAD < 0 & FDR_NORAD <= 0.05, "external_gene_name"], use.names = F), 
MYCN.down = unlist(diff[logFC_MYCN < 0 & FDR_MYCN <= 0.05, "external_gene_name"], use.names = F))
venn(input)

# process data from MYCN knock-down in SK-N-BE(2) (E-GEOD-84389) (correlate hg38 to hg19)
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

# gene set enrichment of differentially expressed genes after NORAD silencing
diff = unlist(NORAD.kd[FDR <= 0.01 & abs(logFC) > 1 & gene_biotype == "protein_coding", "external_gene_name"], use.names = F)
dbs = c("ENCODE_TF_ChIP-seq_2015", "GO_Molecular_Function_2017", "GO_Biological_Process_2017", "Reactome_2016")
enriched <- enrichr(diff, dbs)
# Transcription Factor
enriched.TF = head(enriched[["ENCODE_TF_ChIP-seq_2015"]], 10)
enriched.TF$Term <- gsub("_hg19", "", enriched.TF$Term); enriched.TF$Term <- gsub("_", " ", enriched.TF$Term); enriched.TF$Type <- "ENCODE TF ChIP-seq"
# Molecular Function
enriched.MF = head(enriched[["GO_Molecular_Function_2017"]], 10)
enriched.MF$Term <- gsub("GO.*", "", enriched.MF$Term); enriched.MF$Term <- gsub("(ubiquinone) ", "", enriched.MF$Term, fixed = T); enriched.MF$Term <- gsub(" (", "", enriched.MF$Term, fixed = T)
enriched.MF$Type <- "Molecular Function"
# Biological Process
enriched.BP = head(enriched[["GO_Biological_Process_2017"]], 10)
enriched.BP$Term <- gsub("GO.*", "", enriched.BP$Term); enriched.BP$Term <- gsub(" (", "", enriched.BP$Term, fixed = T)
enriched.BP$Type <- "Biological Process"
# Reactome
enriched.rea = head(enriched[["Reactome_2016"]], 10)
enriched.rea$Term <- gsub("_R-HSA.*", "", enriched.rea$Term)
enriched.rea$Term <- gsub("_Homo sapiens", "", enriched.rea$Term)
enriched.rea$Type <- "Reactome"
enriched.rea <- enriched.rea[enriched.rea$Term != "Cell Cycle, Mitotic", ]
# Combine
enriched <- rbind(enriched.TF[1:5,], enriched.MF[1:5,], enriched.BP[1:5,], enriched.rea[1:5,])
str_width=25
enriched$Term <- str_wrap(enriched$Term, width = str_width)    
enriched$Term <- factor(enriched$Term, levels = rev(enriched$Term))
enriched$Type <- factor(enriched$Type, levels = c("Biological Process", "ENCODE TF ChIP-seq", "Reactome", "Molecular Function"))
p1 <- ggplot(enriched, aes(Term, Combined.Score)) + theme_bw() + geom_col(fill = "red") + coord_flip() + labs(x = "", y = "Combined Score") + facet_wrap(~Type, scales = "free")

#postscript("GSEA.eps", width = 10, height = 5)
p1
#dev.off()

# generate plot over NORAD showing genes that are MYCN regulated and RNA binding
temp <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "source"), values = NORAD.kd$ensembl_gene_id, mart = mart)
temp <- temp[temp$source == "ensembl_havana", ]
NORAD.plot <- NORAD.kd[ensembl_gene_id %in% temp$ensembl_gene_id, ]
NORAD.plot$logFDR <- ifelse(NORAD.plot$logFC < 0, - -log10(NORAD.plot$FDR), -log10(NORAD.plot$FDR))
NORAD.plot <- NORAD.plot[order(logFDR), ] 
NORAD.plot$type <- ifelse(NORAD.plot$external_gene_name %in% unlist(strsplit(enriched.MF$Genes[1], ";")), "RNA binding", "")
topgenes <- c(unlist(NORAD.plot[type == "RNA binding", "external_gene_name"][1:10]), rev(unlist(NORAD.plot[type == "RNA binding", "external_gene_name"]))[1:10])
NORAD.plot$name <- ifelse(NORAD.plot$external_gene_name %in% topgenes, NORAD.plot$external_gene_name, "")
NORAD.plot <- NORAD.plot[order(NORAD.plot$type), ]
p2 <- ggplot(NORAD.plot, aes(logFC, -log10(FDR))) + geom_point(aes(color = type), size = 0.5, show.legend = T) + theme_bw() + theme(panel.grid = element_blank(), legend.position = c(0.15,0.87), legend.title =  element_blank()) + scale_color_manual(values = c("#999999", "#56B4E9")) + geom_text_repel(alpha = 0.7, label = NORAD.plot$name, segment.alpha = 0.3, size = 2, force = 0.4)
t <- merge(MYCN.kd.hg38[MYCN.kd.hg38$FDR < 0.05, ], NORAD.kd[FDR < 0.05, ], by = "external_gene_name")
t <- t[(t$logFC.x * t$logFC.y) > 0, ]
NORAD.plot$type <- ifelse(NORAD.plot$external_gene_name %in% t$external_gene_name, "MYCN BE2 RNA-seq", "")
NORAD.plot <- NORAD.plot[order(logFDR), ] 
topgenes <- c(unlist(NORAD.plot[type == "MYCN BE2 RNA-seq", "external_gene_name"][1:10]), rev(unlist(NORAD.plot[type == "MYCN BE2 RNA-seq", "external_gene_name"]))[1:10])
NORAD.plot$name <- ifelse(NORAD.plot$external_gene_name %in% topgenes, NORAD.plot$external_gene_name, "")
NORAD.plot <- NORAD.plot[order(NORAD.plot$type), ]
p3 <- ggplot(NORAD.plot, aes(logFC, -log10(FDR))) + geom_point(aes(color = type), size = 0.5) + theme_bw() + theme(panel.grid = element_blank(), legend.position = c(0.15,0.87), legend.title =  element_blank()) + scale_color_manual(values = c("#999999", "#E69F00")) + geom_text_repel(alpha = 0.7, label = NORAD.plot$name, segment.alpha = 0.3, size = 2, force = 0.2)
p4 <- grid.arrange(p3, p2, ncol = 1)

#svg("GSEA.volcano.svg", width = 12, height = 7)
grid.arrange(p1, p4, ncol = 2)
#dev.off()

sessionInfo()
