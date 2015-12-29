#load GEO file for study
library(GEOquery)
getGEOSuppFiles("GSE51024")
untar("./GSE51024/GSE51024_RAW.tar", exdir="./GSE51024/data")
