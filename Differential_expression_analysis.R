do_deseq2_new <- function(deseq2_exp_mat, group_df, g1, g2, to_int=FALSE) {
    if(to_int==TRUE){
        rn = rownames(deseq2_exp_mat)
        deseq2_exp_mat = sapply(deseq2_exp_mat,round)
        rownames(deseq2_exp_mat) = rn
    }

    cat(sprintf('foldchange %s : %s',g2,g1))

    group_df = rbind(group_df%>%filter(group%in%c(g1)),group_df%>%filter(group%in%c(g2)))
    group_df = group_df%>%filter(sample%in%colnames(deseq2_exp_mat))
    
    deseq2_exp_mat = deseq2_exp_mat[,group_df$sample]

    group_df[group_df$group==g1,]$group = '1'   
    group_df[group_df$group==g2,]$group = '2'

    condition_table = group_df$group
    
    library(DESeq2)
    # deseq2_exp_mat

    dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                                DataFrame(condition_table), 
                                design= ~ condition_table)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    # acquire the results using function results(), and assign to res
    res <- results(dds2)
    # view the summary of results
    summary(res)
    resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=TRUE)

    vsd <- vst(dds2, blind=FALSE)
    resdata_vst <- merge(as.data.frame(res),SummarizedExperiment::assay(vsd),by="row.names",sort=TRUE)

    # output
    resdata = resdata[order(resdata$padj),]
    resdata_vst = resdata_vst[order(resdata_vst$padj),]

    colnames(resdata)[1] = 'row.names'
    colnames(resdata_vst)[1] = 'row.names'

    return(list(resdata=resdata, resdata_vst=resdata_vst, dds=dds2))
}

