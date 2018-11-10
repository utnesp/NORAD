library(biomaRt)

enst2ensg2_ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", 
    by.y = "ensembl_transcript_id", all = F) 
{
    df <- getBM(filters = "ensembl_transcript_id", attributes = c("ensembl_transcript_id", 
        "ensembl_gene_id", "transcript_biotype", "external_transcript_name", 
        "external_gene_name", "gene_biotype"), values = ensg, 
        mart = biomart)
    if (combine == T) {
        if (df2 == "") {
            if (grepl("$", strsplit(as.character(match.call()), 
                "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*", "", strsplit(as.character(match.call()), 
                  "=")[[2]]))
                by.y = gsub(".*\\$", "", strsplit(as.character(match.call()), 
                  "=")[[2]])
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                  all = all)
            }
            else {
                df2 <- get(strsplit(as.character(match.call()), 
                  "=")[[2]])
                by.y = "row.names"
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                  all = all)
            }
            return(df.anno)
        }
        else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                all = all)
            return(df.anno)
        }
    }
    else {
        return(df)
    }
}



enst2utr <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", 
    by.y = "ensembl_transcript_id", all = F) 
{
    df <- getBM(filters = "ensembl_transcript_id", attributes = c("ensembl_transcript_id", 
        "3_utr_start", "3_utr_end", "5_utr_start", "5_utr_end"), 
        values = ensg, mart = biomart)
    if (combine == T) {
        if (df2 == "") {
            if (grepl("$", strsplit(as.character(match.call()), 
                "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*", "", strsplit(as.character(match.call()), 
                  "=")[[2]]))
                by.y = gsub(".*\\$", "", strsplit(as.character(match.call()), 
                  "=")[[2]])
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                  all = all)
            }
            else {
                df2 <- get(strsplit(as.character(match.call()), 
                  "=")[[2]])
                by.y = "row.names"
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                  all = all)
            }
            return(df.anno)
        }
        else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, 
                all = all)
            return(df.anno)
        }
    }
    else {
        return(df)
    }
}

# get all transcripts with ids, name and biotype for those transcripts that contains UTRs
transcipts <- enst2ensg2_ext_name_biotype("")
transcipts_utrs <- enst2utr(transcipts$ensembl_transcript_id, combine = T)
transcipts_utrs <- transcipts_utrs[!is.na(transcipts_utrs$`3_utr_start`) | !is.na(transcipts_utrs$`3_utr_end`) | !is.na(transcipts_utrs$`5_utr_start`) | !is.na(transcipts_utrs$`5_utr_end`), ]
transcipts_utrs <- transcipts_utrs[!duplicated(transcipts_utrs$ensembl_transcript_id),]
transcipts_utrs <- transcipts_utrs[, c(1,6:9)]

# get 5utr and 3utr sequence for each transcript, and append to file
fasta.output="GRCh38.utr.fa"
fasta.temp="GRCh38.utr.temp.fa"
# split up queries so that we dont heap up memory, and not to small queries so it does not take too long
intervals=5000
remainder=nrow(transcipts_utrs) %% intervals
number_of_times <- (nrow(transcipts_utrs) - remainder)/intervals
# query
for(i in 1:number_of_times) {
    end_row=i*intervals
    start_row=end_row-(intervals-1)
    df=transcipts_utrs[start_row:end_row, ]
    df$name =  paste(df$ensembl_transcript_id, df$ensembl_gene_id, df$transcript_biotype, df$external_transcript_name, df$external_gene_name, sep = "|")
    # 5' UTR
    seq_5utr = getSequence(id=df$ensembl_transcript_id, type="ensembl_transcript_id", seqType="5utr", mart = mart)
    seq_5utr$ensembl_transcript_id <- paste(seq_5utr$ensembl_transcript_id, df[match(seq_5utr$ensembl_transcript_id, df$ensembl_transcript_id), "name"], "5utr", sep = "|")
    seq_5utr <- seq_5utr[nchar(seq_5utr$`5utr`) > 31, ]
    exportFASTA(seq_5utr, fasta.temp)
    # 3' UTR
    seq_3utr = getSequence(id=df$ensembl_transcript_id, type="ensembl_transcript_id", seqType="3utr", mart = mart)
    seq_3utr$ensembl_transcript_id <- paste(seq_3utr$ensembl_transcript_id, df[match(seq_3utr$ensembl_transcript_id, df$ensembl_transcript_id), "name"], "3utr", sep = "|")
    seq_3utr <- seq_3utr[nchar(seq_3utr$`3utr`) > 31, ]
    exportFASTA(seq_3utr, fasta.temp)
}
# Process fasta file
# remove empty lines in fasta
command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.temp, " > ", fasta.output, sep = ""))
system(command)
# compress file using gzip 
system(paste("gzip ", fasta.output, sep=""))
