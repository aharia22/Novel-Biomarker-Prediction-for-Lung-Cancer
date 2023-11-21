if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Load the DESeq2 package
library(DESeq2)

# Load your RNA-seq count data
# Replace 'path/to/GSE103584_R01_NSCLC_RNAseq.txt' with the actual path to your file
counts <- read.delim("GSE103584_R01_NSCLC_RNAseq.txt", header=TRUE, row.names=1)

# Replace NAs with zeros
counts[is.na(counts)] <- 0
counts <- as.data.frame(lapply(counts, as.integer))


# Read in your metadata (sample information)
colData <- read.delim("sample_metadata.txt", header=TRUE, row.names=1)

# Ensure the column names of count data match the row names of colData
# This might involve renaming the columns of counts or the rows of colData
# For example:
colnames(counts) <- colData$sample
colData$treatment <- as.factor(colData$treatment)


# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ treatment)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Running the differential analysis (adjust the contrast as per your experimental design)
res <- results(dds, contrast=c("treatment", "male", "female"))

# Ordering results by significance
resOrdered <- res[order(res$pvalue),]

# Exporting results
write.csv(as.data.frame(resOrdered), file="DESeq2_results.csv")
