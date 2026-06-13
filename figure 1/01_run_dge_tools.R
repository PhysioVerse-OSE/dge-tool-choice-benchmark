suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- ifelse(length(args) >= 1, args[1], "figure1_inputs")
out_dir <- ifelse(length(args) >= 2, args[2], "figure1_dge_results")

dir.create(file.path(out_dir, "edgeR"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DESeq2"), recursive = TRUE, showWarnings = FALSE)

load_counts <- function(filepath) {
  dat <- read.csv(filepath, header = TRUE, check.names = FALSE)
  gene_symbols <- dat[[1]]
  counts <- dat[, -1, drop = FALSE]
  rownames(counts) <- gene_symbols
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  return(counts)
}

detect_groups <- function(column_names) {
  sub(" \\(.*\\)", "", column_names)
}

run_edgeR <- function(counts) {
  groups <- factor(detect_groups(colnames(counts)))
  groups <- relevel(groups, ref = "Control")
  dge <- DGEList(counts = counts, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = dge$samples)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- topTags(qlf, n = Inf)$table
  tab$GeneSymbol <- rownames(tab)
  tab$Bonferroni_pvalue <- p.adjust(tab$PValue, method = "bonferroni")
  tab <- tab[, c("GeneSymbol", setdiff(colnames(tab), "GeneSymbol"))]
  return(tab)
}

run_deseq2 <- function(counts) {
  groups <- factor(detect_groups(colnames(counts)))
  groups <- relevel(groups, ref = "Control")
  coldata <- data.frame(row.names = colnames(counts), group = groups)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "RSVB", "Control"))
  tab <- as.data.frame(res[order(res$pvalue), ])
  tab$GeneSymbol <- rownames(tab)
  colnames(tab)[colnames(tab) == "log2FoldChange"] <- "logFC"
  colnames(tab)[colnames(tab) == "pvalue"] <- "PValue"
  colnames(tab)[colnames(tab) == "padj"] <- "FDR"
  tab$Bonferroni_pvalue <- p.adjust(tab$PValue, method = "bonferroni")
  tab <- tab[, c("GeneSymbol", setdiff(colnames(tab), "GeneSymbol"))]
  return(tab)
}

manifest_path <- file.path(input_dir, "manifest.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
all_files <- unique(manifest$file)

for (f in all_files) {
  message("Processing ", f)
  counts <- load_counts(f)
  stem <- tools::file_path_sans_ext(basename(f))
  edge_out <- file.path(out_dir, "edgeR", paste0("edgeR_", stem, ".csv"))
  deseq_out <- file.path(out_dir, "DESeq2", paste0("DESeq2_", stem, ".csv"))
  write.csv(run_edgeR(counts), edge_out, row.names = FALSE)
  write.csv(run_deseq2(counts), deseq_out, row.names = FALSE)
}

write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
message("Finished DGE analyses. Results written to ", out_dir)
