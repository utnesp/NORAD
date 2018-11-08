library(survival)
library(survminer)
library(data.table)
library(biomaRt)
library(easybiomart)
library(gridExtra)
library(doParallel)
library(foreach)
library(ggrepel)
library(ggsci)

mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', ensemblRedirect = T)

source_https <- function(url, ...) {
  # load package
  require(RCurl)
 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://raw.githubusercontent.com/utnesp/NORAD/master/plot.surv.R")
source_https("https://raw.githubusercontent.com/utnesp/parallell_functions/master/R/coxph.parallell")
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62564/suppl/GSE62564%5FSEQC%5FNB%5FRNA%2DSeq%5Flog2RPM%2Etxt%2Egz" # read GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt.gz from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564
GSE62564_SEQC_NB_RNA.Seq_log2RPM <- fread(paste(paste('URL=',url, sep = ""), ' ; curl "$URL" | gunzip -c', sep = ""))
SEQC.log <- as.data.frame(GSE62564_SEQC_NB_RNA.Seq_log2RPM)
SEQC.log[, 2:ncol(SEQC.log)] <- 2^SEQC.log[, 2:ncol(SEQC.log)]
SEQC.log[, 2:ncol(SEQC.log)] <- log2((SEQC.log[, 2:ncol(SEQC.log)]+1))
SEQC.log <- SEQC.log[rowMeans(SEQC.log[, 2:ncol(SEQC.log)]) > 1, ]
SEQC.ClincicalAttributes <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/SEQC.SampleCharacteristics.txt")
SEQC.ClincicalAttributes <- as.data.frame(SEQC.ClincicalAttributes)
SEQC.log <- as.data.table(SEQC.log)
coxph.parallell(SEQC.log, SEQC.ClincicalAttributes[, c("os_bin", "os_day")], file = "test/coxph.res.txt") 
SEQC.genes.coxph <- fread("test/coxph.res.txt", col.names = c("var", "n", "estimate", "se.coef", "hazard.ratio", "zvalue", "p.wald", "p.lrt", "p.logrank", "rsq", "concordance", "concordance.se", "CI.lower", "CI.upper", "chi.sq.p"))
SEQC.genes.coxph <- refseq2ensg_ext_name(SEQC.genes.coxph$var, combine = T)
SEQC.genes.coxph$external_gene_name <- NULL
SEQC.genes.coxph <- ensg2ext_name_biotype(SEQC.genes.coxph$ensembl_gene_id, combine = T, all = T)
SEQC.genes.coxph$var <- gsub("\\..*","", SEQC.genes.coxph$var)
SEQC.genes.coxph <- ensg2ext_name_biotype(SEQC.genes.coxph$var, combine = T)
SEQC.genes.coxph <- as.data.table(SEQC.genes.coxph)
SEQC.genes.coxph$p.lrt <- ifelse(SEQC.genes.coxph$p.lrt == 0, 1.11e-16, SEQC.genes.coxph$p.lrt)
SEQC.genes.coxph$p.wald <- ifelse(SEQC.genes.coxph$p.wald == 0, 1.11e-16, SEQC.genes.coxph$p.wald)
SEQC.genes.coxph$p.logrank <- ifelse(SEQC.genes.coxph$p.logrank == 0, 1.11e-16, SEQC.genes.coxph$p.logrank)
SEQC.genes.coxph <- SEQC.genes.coxph[order(SEQC.genes.coxph$p.lrt), ]

t <- SEQC.genes.coxph[!duplicated(paste(SEQC.genes.coxph$external_gene_name, SEQC.genes.coxph$gene_biotype, SEQC.genes.coxph$hazard.ratio, SEQC.genes.coxph$zvalue)), ]
t <- t[t$gene_biotype != "protein_coding", ]
t <- t[t$chi.sq.p > 0.05, ]
means <- rowMeans(SEQC.log[,2:ncol(SEQC.log)])
names(means) <- SEQC.log$RefSeqID
t$means <-  means[match(t$refseq, names(means))]
t$name <- ifelse(abs(log2(t$hazard.ratio)) > 1.8 & t$p.logrank < 0.01 & t$means > 3, t$external_gene_name, "")

ggplot(t, aes(log2(hazard.ratio), -log10(t$p.logrank))) +  geom_point(aes(color = means)) + theme_bw() + labs(x= "Hazard ratio (log2)", y = "P value (-log10)") + geom_text_repel(label = t$name) + scale_color_continuous(low = "grey", high = "red", breaks = c(0,2,4,6,8,10), name = "Average\nexpression") + theme(legend.position = c(0.9,0.2), legend.box.background = element_rect(color = "black"))

t <- t[t$name != "", ]
t <- t[order(abs(log2(t$hazard.ratio))), ]
t <- t[!duplicated(t$external_gene_name), ]
t <- t[order(ifelse(t$hazard.ratio < 1, -t$means, t$means)), ]
t <- t[t$gene_biotype !=  "processed_transcript", ]
t <- t[c(1,2,3,(nrow(t)-2):nrow(t)), ]

t2 <- SEQC.log
t2 <- t2[t2$RefSeqID %in% t$refseq, ]
t2$RefSeqID <-  unlist(t[match(t2$RefSeqID, t$refseq), c("ensembl_gene_id")])
colnames(t2)[1] <- "geneid" 
t2 <- as.data.frame(t2)
SEQC.ClincicalAttributes <- as.data.table(SEQC.ClincicalAttributes)
genenames = t$external_gene_name

for (i in 1:length(genenames)) {
    p <- plot.surv(genenames[i], SEQC.ClincicalAttributes[, c("os_bin", "os_day")], t2, cutoff.modus = "scan", log.transform = F)
    pval = ifelse(surv.coxph$p.logrank == 0, "<2e-16", surv.coxph$p.logrank)
    p <- p + ylab("Overall survival probabilty") + ggtitle(paste(genenames[i], "    p=", pval, paste("    LOW=",survdiff$n[1], sep = ""), paste("    HIGH=",survdiff$n[2], sep = ""), sep = ""))
    p <- ggpar(p, legend.title = "", legend = c(0.9,0.2), )
    assign(paste("p", i, sep = ""), p)
}

grid.arrange(p1$plot,p4$plot,
             p2$plot,p5$plot,
             p3$plot,p6$plot,ncol = 2)

sessionInfo()

