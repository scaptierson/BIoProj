library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(grid)

dir.create("results/DE", recursive = TRUE, showWarnings = FALSE)
dir.create("results/DE/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/DE/enrichment", recursive = TRUE, showWarnings = FALSE)

prepared_data <- readRDS("data/processed/hcv_data.rds")
count_matrix <- prepared_data$counts
sample_table <- prepared_data$sample_table
gene_annotation <- prepared_data$gene_annotation

stopifnot(identical(colnames(count_matrix), sample_table$Run))
stopifnot(all(count_matrix >= 0), all(count_matrix == round(count_matrix)))

sample_table$sex <- factor(sample_table$`Sample Characteristic[sex]`)
sample_table$fibrosis_group <- factor(
  sample_table$fibrosis_group,
  levels = c("F0-F2", "F3-F4")
)

design_matrix <- model.matrix(~ sex + fibrosis_group, data = sample_table)
stopifnot(qr(design_matrix)$rank == ncol(design_matrix))

y <- DGEList(counts = count_matrix)
keep_genes <- filterByExpr(y, design = design_matrix)
y <- y[keep_genes, , keep.lib.sizes = FALSE]
y <- normLibSizes(y, method = "TMM")

filtered_gene_ids <- rownames(y)
write.csv(
  data.frame(Gene_ID = filtered_gene_ids),
  "results/DE/filtered_genes_universe.csv",
  row.names = FALSE
)

v <- voom(y, design_matrix, plot = FALSE)
png("results/DE/plots/voom_mean_variance.png", width = 1800, height = 1400, res = 200)
voom(y, design_matrix, plot = TRUE)
dev.off()

fit_voom <- eBayes(lmFit(v, design_matrix))
results_voom <- topTable(
  fit_voom,
  coef = "fibrosis_groupF3-F4",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
results_voom$Gene_ID <- rownames(results_voom)
results_voom$FDR <- results_voom$adj.P.Val

write.csv(
  results_voom,
  "results/DE/limma_voom_F3F4_vs_F0F2.csv",
  row.names = FALSE
)

y_edger <- estimateDisp(y, design_matrix)
fit_edger <- glmQLFit(y_edger, design_matrix)
test_edger <- glmQLFTest(fit_edger, coef = "fibrosis_groupF3-F4")
results_edger <- topTags(
  test_edger,
  n = Inf,
  adjust.method = "BH",
  sort.by = "P"
)$table
results_edger$Gene_ID <- rownames(results_edger)

write.csv(
  results_edger,
  "results/DE/edgeR_F3F4_vs_F0F2.csv",
  row.names = FALSE
)

col_data <- sample_table
rownames(col_data) <- col_data$Run
col_data$sex <- droplevels(col_data$sex)
col_data$fibrosis_group <- droplevels(col_data$fibrosis_group)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix[keep_genes, , drop = FALSE]),
  colData = col_data,
  design = ~ sex + fibrosis_group
)
dds <- DESeq(dds, quiet = TRUE)

results_deseq <- results(
  dds,
  contrast = c("fibrosis_group", "F3-F4", "F0-F2"),
  alpha = 0.05
)
results_deseq <- as.data.frame(results_deseq)
results_deseq$Gene_ID <- rownames(results_deseq)
results_deseq$FDR <- results_deseq$padj

write.csv(
  results_deseq,
  "results/DE/DESeq2_F3F4_vs_F0F2.csv",
  row.names = FALSE
)

threshold_fdr <- 0.05
threshold_log2fc <- 1

sig_voom <- results_voom[
  !is.na(results_voom$FDR) &
    results_voom$FDR < threshold_fdr &
    abs(results_voom$logFC) > threshold_log2fc,
]
sig_edger <- results_edger[
  !is.na(results_edger$FDR) &
    results_edger$FDR < threshold_fdr &
    abs(results_edger$logFC) > threshold_log2fc,
]
sig_deseq <- results_deseq[
  !is.na(results_deseq$FDR) &
    results_deseq$FDR < threshold_fdr &
    abs(results_deseq$log2FoldChange) > threshold_log2fc,
]

cat("\nДифференциальная экспрессия (FDR < 0.05 и |log2FC| > 1):\n")
cat("Исследовано генов:", length(filtered_gene_ids), "\n")
cat("limma-voom:", nrow(sig_voom), "\n")
cat("edgeR:", nrow(sig_edger), "\n")
cat("DESeq2:", nrow(sig_deseq), "\n")

genes_voom <- sig_voom$Gene_ID
genes_edger <- sig_edger$Gene_ID
genes_deseq <- sig_deseq$Gene_ID

common_genes <- Reduce(intersect, list(genes_voom, genes_edger, genes_deseq))
all_significant_genes <- Reduce(union, list(genes_voom, genes_edger, genes_deseq))

method_comparison <- data.frame(
  Gene_ID = all_significant_genes,
  limma_voom = all_significant_genes %in% genes_voom,
  edgeR = all_significant_genes %in% genes_edger,
  DESeq2 = all_significant_genes %in% genes_deseq
)
method_comparison$Methods_Count <- rowSums(
  method_comparison[, c("limma_voom", "edgeR", "DESeq2")]
)

write.csv(method_comparison, "results/DE/method_comparison.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = common_genes),
          "results/DE/common_genes_all_methods.csv", row.names = FALSE)

cat("\nРаспределение генов по количеству методов:\n")
print(table(method_comparison$Methods_Count))
cat("Общие гены всех трёх методов:", length(common_genes), "\n")

voom_fc <- setNames(results_voom$logFC, results_voom$Gene_ID)
edger_fc <- setNames(results_edger$logFC, results_edger$Gene_ID)
deseq_fc <- setNames(results_deseq$log2FoldChange, results_deseq$Gene_ID)

common_results <- data.frame(
  Gene_ID = common_genes,
  limma_voom_log2FC = unname(voom_fc[common_genes]),
  edgeR_log2FC = unname(edger_fc[common_genes]),
  DESeq2_log2FC = unname(deseq_fc[common_genes])
)
common_results$Same_Direction <- with(
  common_results,
  sign(limma_voom_log2FC) == sign(edgeR_log2FC) &
    sign(edgeR_log2FC) == sign(DESeq2_log2FC)
)

consistent_genes <- common_results[
  !is.na(common_results$Same_Direction) & common_results$Same_Direction,
]
consistent_genes$Direction <- ifelse(
  consistent_genes$limma_voom_log2FC > 0,
  "Upregulated",
  "Downregulated"
)

write.csv(common_results, "results/DE/common_genes_log2FC.csv", row.names = FALSE)
write.csv(consistent_genes, "results/DE/consistent_genes.csv", row.names = FALSE)
write.csv(
  consistent_genes[consistent_genes$Direction == "Upregulated", ],
  "results/DE/upregulated_genes.csv",
  row.names = FALSE
)
write.csv(
  consistent_genes[consistent_genes$Direction == "Downregulated", ],
  "results/DE/downregulated_genes.csv",
  row.names = FALSE
)

gene_annotation$`Gene ID` <- as.character(gene_annotation$`Gene ID`)
consistent_genes$Gene_ID <- as.character(consistent_genes$Gene_ID)

consistent_genes_annotated <- merge(
  consistent_genes,
  gene_annotation,
  by.x = "Gene_ID",
  by.y = "Gene ID",
  all.x = TRUE,
  sort = FALSE
)
write.csv(
  consistent_genes_annotated,
  "results/DE/consistent_genes_annotated.csv",
  row.names = FALSE
)

pvalue_data <- rbind(
  data.frame(Method = "limma-voom", PValue = results_voom$P.Value),
  data.frame(Method = "edgeR", PValue = results_edger$PValue),
  data.frame(Method = "DESeq2", PValue = results_deseq$pvalue)
)
pvalue_data <- pvalue_data[!is.na(pvalue_data$PValue), ]

p_hist <- ggplot(pvalue_data, aes(x = PValue)) +
  geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
  facet_wrap(~ Method, ncol = 1) +
  labs(
    title = "Распределение p-значений",
    x = "p-value",
    y = "Количество генов"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "results/DE/plots/pvalue_histograms.png",
  p_hist,
  width = 9,
  height = 10,
  dpi = 300
)

make_volcano <- function(data, logfc_col, fdr_col, method_name, output_file) {
  plot_data <- data.frame(
    Gene_ID = data$Gene_ID,
    log2FC = data[[logfc_col]],
    FDR = data[[fdr_col]]
  )
  plot_data$FDR_plot <- pmax(plot_data$FDR, .Machine$double.xmin, na.rm = TRUE)
  plot_data$Status <- "Not significant"
  plot_data$Status[
    !is.na(plot_data$FDR) &
      plot_data$FDR < threshold_fdr &
      abs(plot_data$log2FC) > threshold_log2fc
  ] <- "Significant"
  
  p <- ggplot(plot_data, aes(x = log2FC, y = -log10(FDR_plot), colour = Status)) +
    geom_point(alpha = 0.65, size = 1) +
    geom_vline(xintercept = c(-threshold_log2fc, threshold_log2fc),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(threshold_fdr), linetype = "dashed") +
    scale_colour_manual(values = c("Not significant" = "grey70",
                                   "Significant" = "firebrick")) +
    labs(title = paste("Volcano plot:", method_name),
         x = "log2 fold change (F3-F4 vs F0-F2)",
         y = "-log10(FDR)",
         colour = "Статус") +
    theme_minimal(base_size = 12)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}

make_volcano(results_voom, "logFC", "FDR", "limma-voom",
             "results/DE/plots/volcano_limma_voom.png")
make_volcano(results_edger, "logFC", "FDR", "edgeR",
             "results/DE/plots/volcano_edgeR.png")
make_volcano(results_deseq, "log2FoldChange", "FDR", "DESeq2",
             "results/DE/plots/volcano_DESeq2.png")

make_ma <- function(data, mean_col, logfc_col, fdr_col, method_name, output_file) {
  plot_data <- data.frame(
    MeanExpression = data[[mean_col]],
    log2FC = data[[logfc_col]],
    FDR = data[[fdr_col]]
  )
  plot_data$Status <- ifelse(
    !is.na(plot_data$FDR) &
      plot_data$FDR < threshold_fdr &
      abs(plot_data$log2FC) > threshold_log2fc,
    "Significant",
    "Not significant"
  )
  
  p <- ggplot(plot_data, aes(x = MeanExpression, y = log2FC, colour = Status)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("Not significant" = "grey70",
                                   "Significant" = "firebrick")) +
    labs(title = paste("MA plot:", method_name),
         x = "Средний уровень экспрессии",
         y = "log2 fold change",
         colour = "Статус") +
    theme_minimal(base_size = 12)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}

make_ma(results_voom, "AveExpr", "logFC", "FDR", "limma-voom",
        "results/DE/plots/MA_limma_voom.png")
make_ma(results_edger, "logCPM", "logFC", "FDR", "edgeR",
        "results/DE/plots/MA_edgeR.png")
make_ma(results_deseq, "baseMean", "log2FoldChange", "FDR", "DESeq2",
        "results/DE/plots/MA_DESeq2.png")

venn_plot <- venn.diagram(
  x = list(
    "limma-voom" = genes_voom,
    "edgeR" = genes_edger,
    "DESeq2" = genes_deseq
  ),
  filename = NULL,
  fill = c("#4E79A7", "#F28E2B", "#59A14F"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.1,
  main = "Пересечение значимых генов"
)
png("results/DE/plots/venn_methods.png", width = 1800, height = 1600, res = 200)
grid.newpage()
grid.draw(venn_plot)
dev.off()

top50_ids <- head(
  results_voom$Gene_ID[order(results_voom$P.Value, na.last = NA)],
  50
)
heatmap_matrix <- v$E[top50_ids, , drop = FALSE]
sample_annotation <- data.frame(
  Fibrosis = sample_table$fibrosis_group
)
rownames(sample_annotation) <- sample_table$Run

gene_labels <- gene_annotation$`Gene Name`[
  match(rownames(heatmap_matrix), gene_annotation$`Gene ID`)
]
gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(heatmap_matrix)[
  is.na(gene_labels) | gene_labels == ""
]
rownames(heatmap_matrix) <- make.unique(gene_labels)

png("results/DE/plots/heatmap_top50.png", width = 2200, height = 2400, res = 220)
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = sample_annotation,
  clustering_method = "complete",
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "Топ-50 генов по p-value (limma-voom)"
)
dev.off()

ensembl_ids <- sub("\\.[0-9]+$", "", filtered_gene_ids)
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids),
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
gene_map <- gene_map[!is.na(gene_map$ENTREZID), ]
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

write.csv(gene_map, "results/DE/ensembl_to_entrez_symbol.csv", row.names = FALSE)

common_ensembl <- sub("\\.[0-9]+$", "", consistent_genes$Gene_ID)
common_map <- gene_map[match(common_ensembl, gene_map$ENSEMBL), ]
target_entrez <- unique(na.omit(common_map$ENTREZID))

universe_map <- gene_map[
  gene_map$ENSEMBL %in% unique(ensembl_ids),
]
universe_entrez <- unique(na.omit(universe_map$ENTREZID))

if (length(target_entrez) > 0 && length(universe_entrez) > 0) {
  go_bp <- enrichGO(
    gene = target_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  go_mf <- enrichGO(
    gene = target_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  go_cc <- enrichGO(
    gene = target_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  write.csv(as.data.frame(go_bp), "results/DE/enrichment/GO_BP.csv", row.names = FALSE)
  write.csv(as.data.frame(go_mf), "results/DE/enrichment/GO_MF.csv", row.names = FALSE)
  write.csv(as.data.frame(go_cc), "results/DE/enrichment/GO_CC.csv", row.names = FALSE)
  
  if (nrow(as.data.frame(go_bp)) > 0) {
    png("results/DE/enrichment/GO_BP_dotplot.png", width = 2000, height = 1600, res = 220)
    print(dotplot(go_bp, showCategory = 20) + ggtitle("GO Biological Process"))
    dev.off()
  }
} else {
  warning("Недостаточно сопоставленных Entrez ID для GO enrichment.")
}

cat("\nАнализ завершён.\n")
cat("Результаты сохранены в results/DE и вложенных папках.\n")
