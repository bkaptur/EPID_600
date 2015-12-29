library(survival)
library(dplyr)
library(limma)
#code adapted from BioStars tutorial here: https://www.biostars.org/p/153013/

#load data from saved files
rna <- read.table("Meso_RNA/MESO.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T,row.names=1,sep="\t")
rna <- rna[-1,]
clinical <- data.frame(t(read.table("Meso_Clinical/MESO.merged_only_clinical_clin_format.txt",header=T, row.names=1, sep="\t")))

#clean RNA file
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]
vm <- function(x){
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x)
  return(ex$E)
}
rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))

#look at rna distribution pre-normalization
hist(rna_vm)

#normalize RNA
scal <- function(x){
  mean_n <- rowMeans(x)
  sd_n <- apply(x,1,sd)
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm)

#look at RNA distribution post-normalization
hist(z_rna)

#set row names to gene names
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,"\\|"))[[1]])

#match patient info, output number of useable patients
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(clinical$IDs %in% colnames(z_rna))

#retain useful columns
ind_keep <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))

#condense followups
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if(sum(is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- max(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,"NA")
  }
}
ind_keep <- grep("days_to_death",colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if(sum(is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,"NA")
  }
}
ind_keep <- grep("days_to_last_followup",colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if(sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- min(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,"NA")
  }
}
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")

#create vectors for clinical parameters
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                  as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

#output number alive and dead patients
table(clinical$patient.vital_status)

all_clin$death_event <- ifelse(clinical$patient.vital_status == "alive", 0,1)

#add clinical IDs as row names
rownames(all_clin) <- clinical$IDs

#use t-score of 1 to define low expression levels
event_rna <- t(apply(z_rna, 1, function(x) ifelse(x < -1,1,0)))
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))


#create list of genes of interest
gene_list <- c("LMNB2", "KPNA2", "UHRF1", "GLT25D1","MYBL2")
gene_list_length <- length(gene_list)

#use loop to iterate through genes of interest
for (i in 1:gene_list_length){
  #define gene of interest for graph
  ind_gene <- which(rownames(z_rna) == gene_list[i])
  
  #check how many samples are altered
  sample_dist <- table(event_rna[ind_gene,])
  
  #perform survival analysis
  s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
  s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))
  #find p-value
  p_val <- ifelse(is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]
  p_val
  
  #graph survival curves
  plot(survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]),
       col=c(1:3), frame=F, lwd=2,main=paste("MESO",rownames(z_rna)[ind_gene],sep="\n"))
  x1 <- ifelse(is.na(as.numeric(summary(s)$table[,'median'][1])),"NA",as.numeric(summary(s)$table[,'median'][1]))
  x2 <- as.numeric(summary(s)$table[,'median'][2])
  if(x1 != "NA" & x2 != "NA"){
    lines(c(0,max(x1,x2)),c(0.5,0.5),col="blue")
    lines(c(x1,x1),c(0,0.5),col="black")
    lines(c(x2,x2),c(0,0.5),col="red")
  }
  
  #add legend to plot
  legend(max(as.numeric(as.character(all_clin$death_days)[ind_clin]),na.rm = T)*0.3,1,
         legend=c(paste("Normal, n =", sample_dist[1]),paste("Low, n =", sample_dist[2])),bty="n",cex=1.4,lwd=3,col=c("black","red"))
  title(xlab = "Days", ylab = "Overall Survival")
}