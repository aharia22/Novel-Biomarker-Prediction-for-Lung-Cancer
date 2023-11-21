if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#load the DESeq2 package
library(DESeq2)

#load your RNA-seq count data
counts <- read.delim("GSE103584_R01_NSCLC_RNAseq.txt", header=TRUE, row.names=1)

#replace NAs with zeros
counts[is.na(counts)] <- 0
counts <- as.data.frame(lapply(counts, as.integer))


#read in metadata (sample information)
colData <- read.delim("sample_metadata.txt", header=TRUE, row.names=1)

#making sure that the column names of count data match the row names of colData
colnames(counts) <- colData$sample
colData$treatment <- as.factor(colData$treatment)


#create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ treatment)

#run the DESeq pipeline
dds <- DESeq(dds)

#running the differential analysis (adjust the contrast as per your experimental design)
res <- results(dds, contrast=c("treatment", "male", "female"))

#ordering results by significance
resOrdered <- res[order(res$pvalue),]

#exporting results to csv
write.csv(as.data.frame(resOrdered), file="DESeq2_results.csv")
