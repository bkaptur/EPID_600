untar("./GSE51024/GSE51024_RAW.tar", exdir="./GSE51024/data")


library(affy)
library(limma)
library(hgu133plus2.db)

#load and label raw data
raw.data <- read.affybatch(list.celfiles("./GSE51024/data/", full.names=TRUE), compress=TRUE)
pData(raw.data)

pData(raw.data)$status <- c(rep(c("Normal","Tumor"),23),"Tumor","Tumor",rep(c("Normal","Tumor"), 5),"Tumor","Tumor",rep(c("Normal","Tumor"), 4),"Tumor", "Tumor", rep(c("Normal","Tumor"),3),rep("Tumor",8),rep(c("Normal","Tumor"),6))

pData(raw.data)



#boxplot(exprs(raw.data), col="red",main="Raw Probe Intensities")
#boxplot(raw.data, col="red",main="Raw Probe Intensities")

#quality control
GSE51024.rma <- rma(raw.data)

#boxplot(exprs(GSE51024.rma), col="blue", main="RMA Expression Values")

GSE51024.qc <- raw.data[, !sampleNames(raw.data) %in% c("GSM1235084_Tumor42.CEL.gz", "GSM1235085_Tumor43.CEL.gz",
                                                        "GSM1235092_Tumor44.CEL.gz", "GSM1235093_Tumor45.CEL.gz",
                                                        "GSM1235094_Tumor46.CEL.gz", "GSM1235095_Tumor47.CEL.gz",
                                                        "GSM1235096_Tumor48.CEL.gz", "GSM1235097_Tumor49.CEL.gz",
                                                        "GSM1235098_Tumor50.CEL.gz", "GSM1235099_Tumor51.CEL.gz",
                                                        "GSM1235062_Tumor52.CEL.gz", "GSM1235063_Tumor53.CEL.gz",
                                                        "GSM1235074_Tumor54.CEL.gz", "GSM1235075_Tumor55.CEL.gz") ]


#use sibship pairs for paired study design
pData(GSE51024.qc)$SibShip <- ceiling(c(1:82)/2)
GSE51024.qc.rma <- rma(GSE51024.qc)
SibShip <- factor(pData(GSE51024.qc.rma)$SibShip)
Status <- factor(pData(GSE51024.qc.rma)$status, levels = c("Normal","Tumor"))

#design incorporates sipship and status variables
design <- model.matrix(~-1+SibShip+Status)

fit <- lmFit(GSE51024.qc.rma, design)

#adjust fit coefficients using empirical Bayes moderation of standard errors
fit2 <- eBayes(fit)

#output hypothesis test results
Tumor_results <- topTable(fit2, coef ="StatusTumor" , adjust="BH", num=100, p.value=0.05)
head(Tumor_results)

#load reference gene IDs
library(hgu133plus2.db)
Tumor_results$ID = row.names(Tumor_results)
Tumor_results$SYMBOL <- lapply(Tumor_results$ID, function(x) mget(x, env=hgu133plus2SYMBOL, ifnotfound=NA)[[1]])
head(Tumor_results)
cat('There are', nrow(Tumor_results), 'significant probes.')

#output heatmap with status label at top of plot
top.eset <- GSE51024.qc.rma[row.names(exprs(GSE51024.qc.rma)) %in% row.names(Tumor_results)]
treatment.colors <- unlist(lapply(GSE51024.qc.rma$status, function(x){if (x=="Tumor") "red" else "blue"}))
heatmap(exprs(top.eset), col=topo.colors(100), ColSideColors=treatment.colors)
