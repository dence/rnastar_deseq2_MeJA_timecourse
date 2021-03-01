#modification for working in hipergator cluster_jobid
#Daniel Ence
#February 4, Jan_2021_units

args=commandArgs(trailingOnly=TRUE)
rds=args[1]
output_table=args[2]
output_ma_plot=args[3]
params.config_yaml=args[4]
log=args[5]
threads=args[6]

log <- file(log, open="wt")
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

dds <- readRDS(rds)

library("yaml")
config_yaml <- read_yaml(params.config_yaml)
contrast <- c("tissue", as.vector(config_yaml$diffexp$contrasts$"pWVK305-vs-17544Ri"))
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res,type="normal")
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(output_ma_plot)
plotMA(res, ylim=c(-8,8))
dev.off()

write.table(as.data.frame(res), file=output_table)
