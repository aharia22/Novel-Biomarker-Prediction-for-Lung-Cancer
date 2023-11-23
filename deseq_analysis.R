if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#load the DESeq2 package
library(DESeq2)

#load RNA-seq count data
counts <- read.delim("GSE103584_R01_NSCLC_RNAseq.txt", header=TRUE, row.names=1)

#replace NAs with zeros
counts[is.na(counts)] <- 0
#apply the integer conversion to all columns except the first
counts[, -1] <- lapply(counts[, -1], as.integer)
counts$`R01.023` <- as.integer(counts$`R01.023`)
counts$`R01.024` <- as.integer(counts$`R01.024`)


#read in metadata (sample information)
colData <- read.delim("sample_metadata.txt", header=TRUE, row.names=1)
rownames(colData) <- gsub("-", ".", rownames(colData))

#then we use these row names as the column names for counts
colnames(counts) <- rownames(colData)
#Convert the treatment column in colData to a factor
colData$treatment <- as.factor(colData$treatment)

#check for NA values in the counts matrix
if (any(is.na(counts))) {
  cat("NA values found in the counts matrix.\n")
  # Replace NA values with zeros
  counts[is.na(counts)] <- 0
}

#confirm that there are no more NA values
if (any(is.na(counts))) {
  cat("NA values are still present after replacement.\n")
} else {
  cat("No NA values are present. Proceeding with DESeq2.\n")
}

#proceed with DESeq2 if no NAs are present
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ treatment)

#run the DESeq pipeline
dds <- DESeq(dds)

#running the differential analysis
res <- results(dds, contrast=c("treatment", "male", "female"))

#ordering results by significance
resOrdered <- res[order(res$pvalue),]

#exporting results to csv
write.csv(as.data.frame(resOrdered), file="DESeq2_results.csv", row.names=TRUE)
