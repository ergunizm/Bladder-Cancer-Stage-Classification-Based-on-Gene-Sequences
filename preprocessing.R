if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)

my.gse <- "GSE244266"

my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./dataset", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
my.geo.gse
my.geo.gse <- my.geo.gse[[1]]

untar(paste0("dataset/",my.gse,"_RAW.tar"), exdir=paste0("dataset/CEL"))
my.cels <- list.files(paste0("dataset/CEL"))

my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
head(my.pdata)
dim(my.pdata)
colnames(my.pdata)

head(my.cels)
head(rownames(my.pdata))

rownames(my.pdata) <- my.cels
table(rownames(my.pdata) == my.cels)

write.table(my.pdata, file=paste0("dataset/CEL/",my.gse,"_PhenoData.txt"), sep="\t", quote=F)

BiocManager::install("affy")
library(affy)

cel.path <- paste0("dataset/CEL")
my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse,"_PhenoData.txt"), sep="/"))
show(my.affy)
exprs(my.affy)
head(exprs(my.affy))

colnames(pData(my.affy))
pData(my.affy)
pData(my.affy)$title

my.rma <- rma(my.affy, normalize=T, background=T)
head(exprs(my.rma))

my.rma@annotation
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
BiocManager::install("annotate")
library(annotate)
BiocManager::install("R2HTML")
library(R2HTML)

ID <- featureNames(my.rma)
Symbol<-getSYMBOL(ID,"hgu133plus2.db")
sym <- as.data.frame(Symbol)

data <- as.data.frame(exprs(my.rma)) 
data <- cbind(sym,data)

i <- which(is.na(data$Symbol) == TRUE)
data<-data[-c(i),]

library(data.table)

X <- data.table::as.data.table(data)
final_data <- X[,lapply(.SD,mean),"Symbol"]
final_data <- as.data.frame(final_data)
rownames(final_data) <- final_data[,1] 
final_data <- final_data[,-c(1)]

saveRDS(final_data,"dataset/GSE244266_raw.RDS")
final_data = readRDS("dataset/GSE244266_raw.RDS")

final_data  = as.data.frame(t(final_data))
metadata = pData(my.affy)
table (rownames(final_data) == rownames(metadata))

final_data$stage = substr(metadata$clinical_stage_strat_factor.ch1, 18, 18)
final_data$stage = as.numeric(final_data$stage)

write.csv(final_data, file="dataset/GSE244266.csv")