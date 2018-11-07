plot.surv <- function(genename, clindata, counts, sub.idx = NULL, textsize = 18, count.index = NA, log.transform = T, cutoff.modus = c("mean", "scan", "median", "Kmeans", "user.defined"),scan.nr = 5, user.defined = NULL, max.time = NA, min.time = NA, ...) {
        if(log.transform == T) cat("Log transforming data with pseudocount 1\n")
        if(log.transform == F) cat("Using raw data. Not log transforming\n")
         
        if(length(cutoff.modus) == 5) cutoff.modus = "scan"
        if(!is.null(user.defined)) cutoff.modus = "user.defined"
        cat(paste("Using", cutoff.modus, "as cutoff.modus\n"))
        cat("Cutoff modus may be mean, median, Kmeans, scan or user.defined\n")
        
        cat(paste("Testing", genename, "\n"))
        colnames(counts)[1] <- "ensembl_gene_id"
        
        if (is.na(count.index)) {
            survdata <- cbind(clindata, t(counts[counts$ensembl_gene_id == hugo2ensg(genename)[1,2], 2:ncol(counts)])) 
        } else {
            survdata <- cbind(clindata, t(counts[count.index, 2:ncol(counts)]))
        }
        
        survdata <- survdata[, 1:3]
        colnames(survdata) <- c("sbin", "time", "var")

         if (cutoff.modus == "user.defined") {
                survdata$bin <- as.factor(user.defined)
         } 
        
        
        if (!is.null(sub.idx)) {
            survdata <- survdata[sub.idx, ]    
            survdata$bin <- factor(as.character(survdata$bin)) # redo this to remove any levels that should not be present
        }
        
        if (!is.na(max.time)) survdata <- survdata[survdata$time <= max.time, ]
        if (!is.na(min.time)) survdata <- survdata[survdata$time >= min.time, ]
        
        # order survdata in increasing expression
        survdata <- survdata[order(survdata$var), ]
        
        if(log.transform == T) survdata$var <- log2(survdata$var+1)
        
            # categorize data based on of the function scan, kmeans, median or mean
            if(cutoff.modus == "scan") {
                survdata$bin = 0
                survdata$bin[1:scan.nr] <- 1 # set 1 for the first lowly expressed samples
        
                pval <- list(pval = rep(1,length(survdata$bin)))
                raw.pval <- list(pval = rep(1,length(survdata$bin)))
                
                for(i in 1:nrow(survdata)) {
                    res.cox <- suppressWarnings(coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
                    
                    s <- summary(res.cox)
                    raw.pval$pval[i] <- ifelse(s$logtest[3] == 0, 2e-16, s$logtest[3])
                    pval$pval[i] <- p.adjust(raw.pval$pval[i], method = "bonferroni", n = i)
                    survdata$bin[i] <- 1
                }
                
                # we would not like to select among first and last samples 
                choose.from = (scan.nr+1):(length(raw.pval$pval)-scan.nr-1)
                pvalues = which(pval$pval[choose.from] == min(pval$pval[choose.from]))
                # if all pvalues were adjusted to 1, we choose one of from raw pvalues
                use.raw.pvalues=F
                if(length(pvalues) > 1) use.raw.pvalues = T
                if(use.raw.pvalues == T) pvalues = which(raw.pval$pval[choose.from] == min(raw.pval$pval[choose.from]))
                pvalues <- pvalues+(scan.nr)
                
                survdata$bin <- 2
                survdata$bin[1:(pvalues[length(pvalues)]-1)] <- 1
                col =  c(rep("grey", table(survdata$bin)[1]), rep("black", nrow(survdata) - table(survdata$bin)[1]))
                
                if(use.raw.pvalues == F) df.bar <- barplot(-log10(pval$pval), ylab = "-log10 PValue", xlab = "sample", col = col, border = col)
                if(use.raw.pvalues == T) df.bar <- barplot(-log10(raw.pval$pval), ylab = "-log10 PValue", xlab = "sample", col = col, border = col)
                # plot points of raw values (not log transformed)
                if(log.transform == T) point.exps = 2^(survdata$var-1)
                if(log.transform == F) point.exps = survdata$var
                if(use.raw.pvalues == F) point.exps <- point.exps / (max(point.exps) / max(-log10(pval$pval)))
                if(use.raw.pvalues == T) point.exps <- point.exps / (max(point.exps) / max(-log10(raw.pval$pval)))
                points(df.bar, point.exps, col = "orange", pch = 15)
                title(paste("Scan: logrank iterative test\n", "P=", signif(min(raw.pval$pval),2),"P.bonf=", signif(min(pval$pval),2)))
            }
            
            if (cutoff.modus == "Kmeans") {
                k.fit <- kmeans(survdata[,3], 2, iter.max = 100, nstart = 50) 
                # get cluster means
                agg <- aggregate(survdata[,3],by=list(k.fit$cluster),FUN=mean)
                # append cluster assignment
                cl <- data.frame(survdata[,3], k.fit$cluster)
                cl <- cl[order(cl$k.fit.cluster),]
                cl <- cl[order(cl[,1]), ]
                survdata <- cbind(survdata, cl[, 2])
                colnames(survdata)[4] <- "bin"
            }
            
            if (cutoff.modus == "median") {
                if (log.transform == T) median.sample <- which(round(survdata[,3], 1) == median( unlist(round(survdata[,3], 1))))
                if (log.transform == F) median.sample <- which(round(survdata[,3], 0) == median( unlist(round(survdata[,3], 0))))
                
                survdata$bin <- 2
                survdata$bin[1:length(median.sample)] <- 1
            }
            
            if (cutoff.modus == "mean") {
                mean.sample <- which(survdata[,3] <= mean(unlist(survdata[,3])))
                survdata$bin <- 2
                survdata$bin[1:length(mean.sample)] <- 1
            }
        
                    if (cutoff.modus != "user.defined") {
                        # Categorize into LOW or HIGH
                        if(mean(unlist(survdata[survdata$bin == 1, "var"])) > mean(unlist(survdata[survdata$bin == 2, "var"]))) {
                            survdata$bin <- gsub("1", "HIGH", survdata$bin) 
                            survdata$bin <- gsub("2", "LOW", survdata$bin)
                        } else {
                            survdata$bin <- gsub("1", "LOW", survdata$bin) 
                            survdata$bin <- gsub("2", "HIGH", survdata$bin)
                        }
                    
                        survdata$bin <- factor(survdata$bin, levels = c("LOW", "HIGH"))
                    }
        
        print(table(survdata$bin))
        assign("survdata", survdata, envir = .GlobalEnv)
        
        fit <- survfit(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        assign("survfit", fit, envir = .GlobalEnv)
        print(coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
        
        survdiff <- survdiff(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        assign("survdiff", survdiff, envir = .GlobalEnv)
        print(survdiff(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
        
        case.table = table(survdata$bin)
        
        res.cox <- coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        chi.sq.p <- cox.zph(res.cox) # test the proportional hazards assumption of cox regression
                
        s <- summary(res.cox)
                
                res <- list(var=genename, 
                    n=sum(case.table), 
                    yes.percent=round(case.table[2]/sum(case.table),2)*100, 
                    no.percent=round(case.table[1]/sum(case.table),2)*100, 
                    estimate = round(s$coefficients[1,1],2), 
                    se.coef=round(s$coefficients[1,3],2), 
                    hazard.ratio=round(s$coefficients[1,2],2), 
                    zvalue=round(s$coefficients[1,4],2), 
                    p.wald=signif(s$coefficients[1,5], 3),  
                    p.lrt = signif(s$logtest[3], 3),
                    p.logrank = signif(s$sctest[3], 3),
                    rsq=round(s$rsq[1], 3),
                    concordance = round(s$concordance[1], 2),
                    concordance.se = round(s$concordance[2], 2),
                    CI.lower = round(s$conf.int[3],2), 
                    CI.upper = round(s$conf.int[4],2),
                    chi.sq.p = signif(chi.sq.p$table[3], 3))
        
        assign("surv.coxph", res, envir = .GlobalEnv)  
        assign("surv.coxph.summary", s, envir = .GlobalEnv)  
        return(ggsurvplot(fit, ...))
}
