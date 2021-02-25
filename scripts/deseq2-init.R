#modification for working in hipergator cluster environment
#Daniel ence
#October 22, 2020
args=commandArgs(trailingOnly=TRUE)
counts_file=args[1]
output_file=args[2]
sample_file=args[3]
log_file=args[4]
threads=args[5]



log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(counts_file, header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(sample_file, header=TRUE, row.names="sample", check.names=FALSE)

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output_file)
