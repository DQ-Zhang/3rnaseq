args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(dplyr)
library(data.table)

if (args[1]=="h"){
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
}else if(args[1]=="m"){
  library(org.Mm.eg.db)
  orgdb <- org.Mm.eg.db
}else if(args[1]=="r"){
  library(org.Rn.eg.db)
  orgdb <- org.Rn.eg.db
}else{
  stop("species not supported! {h/m/r}")
}

##DATA PRE-PROCESSING
#compile counts matrix from STAR alignment results
files <- list.files(path = "output/aligned", pattern = "*_counts.tab$", full.names = T, recursive = T)
sample.dirs <- list.dirs(path = "output/aligned", full.names = TRUE, recursive = F)
sample.names <- gsub(paste0("output/aligned/"), "", sample.dirs)

counts <- data.frame(fread(files[1]))[c(1,2)]
for(i in 2:length(files)) {
  counts <- cbind(counts, data.frame(fread(files[i]))[2])
}

#Seperate last 5 lines, which contains count QC metrics, from the rest of the data.frame
row_num<-dim(counts)[1]
counts_qc <- counts[(row_num-4):row_num,]
counts <- counts[1:(row_num-4),]
#Assign column names to counts matrix. The first column are EMSEMBL IDs, remaining columns are counts data of samples.
colnames(counts) <- c("ENSEMBL", sample.names)
#Assign column names counts qc metrics info. The first column are metrics, remaining columns are counts data of samples.
colnames(counts_qc) <- c("metrics", sample.names)
#Convert EMSEMBL ID to gene SYMBOL

counts$Gene <- mapIds(orgdb, keys=as.vector(gsub("\\..*","",counts[["ENSEMBL"]])), 
                      column="SYMBOL", keytype="ENSEMBL",multiVals = "first")
counts <- na.omit(counts)
#Collapse rows with the same gene SYMBOL
counts.matrix <- aggregate(counts[,2:(ncol(counts)-1)],by=list(counts$Gene),FUN=sum)
#Remove "undetermined" sample
if ("undetermined" %in% colnames(counts.matrix)){
  counts.matrix <- subset(counts.matrix, select=-c(undetermined))
}
#Make gene SYMBOL row names
rownames(counts.matrix)=counts.matrix[,1]
counts.matrix <- counts.matrix[,c(2:ncol(counts.matrix))]
#Save counts matrix
dir.create("output/analysis", showWarnings = FALSE)
dir.create("output/analysis/counts", showWarnings = FALSE)
write.table(counts.matrix, file="output/analysis/counts/counts.matrix.txt", quote=F,sep="\t",row.names=T, col.names=T)
write.table(counts_qc, file="output/analysis/counts/counts_qc.txt", quote=F,sep="\t",row.names=T, col.names=T)
