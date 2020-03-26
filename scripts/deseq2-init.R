log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

#### From https://github.com/snakemake-workflows/rna-seq-star-deseq2 

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix -- This not allways work
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

#now they are forced to have the same order
coldata$samples <- rownames(coldata)
coldata <- coldata[match(colnames(cts), coldata$samples), ]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
© 2020 GitHub, Inc.
