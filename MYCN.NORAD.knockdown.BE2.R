library(data.table)
library(ggplot2)
library(gplots)
library(pheatmap)

# read in RNA-seq data after NORAD knock-down in SK-N-BE(2) (supplementary file 2 in manuscript)
NORAD.kd <- fread("https://raw.githubusercontent.com/utnesp/NORAD/master/NORAD.knockdown.CPM.BE2.txt")
setkey(NORAD.kd, external_gene_name)

# process data from MYCN knock-down in SK-N-BE(2) (E-GEOD-84389)
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

# venn diagram showing genes common and specific for NORAD and MYCN knock-down 
fc_threshold = 0.5
FDR_threshold = 0.05
range_threshold = 100
diff = merge(NORAD.kd, MYCN.kd, by = "external_gene_name") 
diff <- diff[abs(log2FC) > fc_threshold & abs(logFC) > fc_threshold, ]
par(mfrow = c(1,2))
input = list(
NORAD.up = unlist(diff[logFC > fc_threshold & FDR.x <= FDR_threshold, "external_gene_name"], use.names = F), 
MYCN.up = unlist(diff[log2FC > fc_threshold & FDR.y <= FDR_threshold, "external_gene_name"], use.names = F))
venn(input)
input = list(
NORAD.down = unlist(diff[logFC < -fc_threshold & FDR.x <= FDR_threshold, "external_gene_name"], use.names = F), 
MYCN.down = unlist(diff[log2FC < -fc_threshold & FDR.y <= FDR_threshold, "external_gene_name"], use.names = F))
venn(input)

## 2d-plot showing both common and specific genes after NORAD and MYCN knock-down (not included in manuscript)
diff = merge(NORAD.kd[abs(logFC) > fc_threshold], 
MYCN.kd[abs(log2FC) > fc_threshold], by = "external_gene_name") 
diff$NORAD.FDR <- ifelse(diff$range.x < 0, - -log10(diff$FDR.x), -log10(diff$FDR.x))
diff$MYCN.FDR <- ifelse(diff$range.y < 0, - -log10(diff$FDR.y), -log10(diff$FDR.y))
diff$color <- ifelse(diff$FDR.x <= 0.05 & diff$FDR.y <= 0.05, "black", "grey")
threshold = -log10(0.05)
diff$name <- ifelse((abs(diff$MYCN.FDR) >= threshold & abs(diff$NORAD.FDR) >= threshold), diff$external_gene_name, "")
ggplot(diff, aes(NORAD.FDR, MYCN.FDR)) + theme_bw() + geom_point(color = "lightgrey") + geom_text(label=diff$name, size = 2)


# heatmap showing common genes after NORAD and MYCN knock-down 
diff <- diff[(diff$log2FC * diff$logFC) > 0, ]
diff <- diff[abs(log2FC) > fc_threshold & abs(logFC) > fc_threshold]
diff <- diff[abs(diff$NORAD.FDR) >= -log10(FDR_threshold) & abs(diff$MYCN.FDR) >= -log10(FDR_threshold), ]
t <- data.frame(NORAD = diff$NORAD.FDR, MYCN = diff$MYCN.FDR, row.names = diff$external_gene_name)
t <- t[order(rowMeans(t)), ]
setEPS()
postscript("MYCN.NORAD.eps", height = 14, width = 2)
pheatmap(t, scale = "none", cluster_rows = F, cluster_cols = F, border_color = NA, cellwidth = 20, cellheight = 6, fontsize_row =  8)
dev.off()

# combine and edit heatmap and venn diagrams in Inkscape
