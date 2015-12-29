#output boxplots of expression values
boxplot(exprs(GSE51024.qc.rma)["204975_at" , ]~ GSE51024.qc.rma$status, col="lightgreen", xlab='EMP2', ylab='Normalized Expression')
boxplot(exprs(GSE51024.qc.rma)["202286_s_at" , ]~ GSE51024.qc.rma$status, col="lightgreen", xlab='TACSTD2', ylab='Normalized Expression')
boxplot(exprs(GSE51024.qc.rma)["206069_s_at" , ]~ GSE51024.qc.rma$status, col="lightgreen", xlab='ACADL', ylab='Normalized Expression')

boxplot(exprs(GSE51024.qc.rma)["209875_s_at" , ]~ GSE51024.qc.rma$status, col="orange", xlab='SPP1', ylab='Normalized Expression')
boxplot(exprs(GSE51024.qc.rma)["201291_s_at" , ]~ GSE51024.qc.rma$status, col="orange", xlab='TOP2A', ylab='Normalized Expression')


#for obtaining gene names for probe
Tumor_results$SYMBOL[Tumor_results$ID=='201291_s_at']
head(Tumor_results)


#for outputting list for network analysis of significant genes and processing 
l_results <- (Tumor_results$SYMBOL)
l_results <- l_results[!is.na(l_results)]
l_results <- unique(l_results)
length(l_results)
lapply(l_results, write, "GEO_genes.txt", append=TRUE, ncolumns=100)


