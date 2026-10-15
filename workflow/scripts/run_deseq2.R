#!/usr/bin/env Rscript
.libPaths(c("/content/rnaseq_env/lib/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
})

option_list <- list(
  make_option(c("--samples"), type="character", help="Path to samples.tsv"),
  make_option(c("--tx2gene"), type="character", help="Path to tx2gene.tsv"),
  make_option(c("--ref"), type="character", default="control", help="Reference condition"),
  make_option(c("--target"), type="character", default="treated", help="Target condition"),
  make_option(c("--out-results"), dest="out_results", type="character", help="Output DESeq2 TSV"),
  make_option(c("--out-counts"), dest="out_counts", type="character", help="Output normalized counts TSV")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load metadata and transcript-to-gene map
samples <- read.delim(opt$samples, stringsAsFactors = FALSE)
rownames(samples) <- samples$sample_id
tx2gene <- read.delim(opt$tx2gene, stringsAsFactors = FALSE)

# 2. Point to Salmon quant.sf files and verify existence
files <- file.path("results", "salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
stopifnot(all(file.exists(files)))

# 3. Summarize transcript quantifications to gene level
# Check if transcript IDs in quant.sf match tx2gene directly or after stripping version
first_quant <- read.delim(files[1], nrows = 5, stringsAsFactors = FALSE)
has_direct_match <- any(first_quant$Name %in% tx2gene$TXNAME)
if (!has_direct_match) {
  tx2gene$TXNAME <- sub("\\..*$", "", tx2gene$TXNAME)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, dropInfReps = TRUE)
} else {
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, dropInfReps = TRUE)
}

# 4. Construct DESeqDataSet & explicitly set reference factor level
samples$condition <- factor(samples$condition)
samples$condition <- relevel(samples$condition, ref = opt$ref)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# Pre-filter low-count genes (at least 10 reads across all samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 5. Run negative binomial Wald test
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", opt$target, opt$ref))
res_Ordered <- as.data.frame(res[order(res$pvalue), ])
res_Ordered$gene_id <- rownames(res_Ordered)
res_Ordered <- res_Ordered[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# 6. Export DESeq2 table and size-factor normalized counts
dir.create(dirname(opt$out_results), recursive = TRUE, showWarnings = FALSE)
write.table(res_Ordered, file = opt$out_results, sep = "\t", quote = FALSE, row.names = FALSE)

norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$gene_id <- rownames(norm_counts)
norm_counts <- norm_counts[, c("gene_id", samples$sample_id)]
write.table(norm_counts, file = opt$out_counts, sep = "\t", quote = FALSE, row.names = FALSE)

message("DESeq2 analysis complete. Tested genes: ", nrow(res_Ordered))
