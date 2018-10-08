getSeriesMatrixCharacteristics <- function(ftpURLof_series_matrix.txt.gz) {
    file = gsub(".*/", "", ftpURLof_series_matrix.txt.gz)
    file = paste(getwd(), file, sep = .Platform$file.sep)
    system(paste("wget ", ftpURLof_series_matrix.txt.gz, " -O ", file, sep = ""))
    system(paste("gzip -f -d ", file, sep = ""))
    file = gsub(".gz", "", file)
    temp.file = readLines(file)
    headers <- read.delim2(file, header = F, skip = grep('!Sample_title', temp.file, fixed = T)-1, nrows = 1)
    headers <- unlist(headers, recursive = F, use.names = F)
    headers <- gsub("\\ .*", "", headers)
    skip=grep('!Sample_characteristics_ch1', temp.file)
    series_matrix <- read.delim2(file, header = F, skip = skip[1]-1, nrows = length(skip), col.names = headers)
    colnames(series_matrix) <- c("characteristic", colnames(SEQC)[2:ncol(SEQC)])
    series_matrix$characteristic <- gsub("\\:.*","", series_matrix$SEQC_NB001)
    n <- colnames(series_matrix)
    series_matrix <- as.data.table(t(series_matrix))
    series_matrix <- cbind(n, series_matrix)
    
    colnames(series_matrix) <- unlist(series_matrix[1])
    series_matrix <- series_matrix[2:nrow(series_matrix), ]
    series_matrix <- as.data.frame(series_matrix)
    row.names(series_matrix) <- series_matrix$characteristic
    series_matrix$characteristic <- NULL
    
    for (i in 1:nrow(series_matrix)) {
        temp <- unlist(series_matrix[i, ], use.names = F)[match(colnames(series_matrix), gsub("\\:.*","", unlist(series_matrix[i, ], use.names = F)))]
        if(i == 1) temp.all <- temp
        if(i != 1) temp.all <- rbind(temp.all, temp)
    }
    
    row.names(temp.all) <- row.names(series_matrix)
    colnames(temp.all) <- colnames(series_matrix)
    series_matrix <- temp.all
    series_matrix <- apply(series_matrix, 2, function(x) gsub(".*\\:","", x))
    series_matrix <- cbind(headers[2:length(headers)], series_matrix)
    colnames(series_matrix)[1] <- "sample"
    series_matrix <- as.data.table(series_matrix)
    colnames(series_matrix) <- gsub(" ", "_", colnames(series_matrix))
    series_matrix <- apply(series_matrix, 2, function(x) gsub(" ", "", x))
    series_matrix <- as.data.frame(series_matrix, stringAsFactors = F)

    return(series_matrix)
}
