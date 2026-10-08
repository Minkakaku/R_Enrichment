# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# GSVA 独立分析入口
# 输入：注释后 Seurat RDS 和 GMT；先按样本×细胞类型汇总表达，再计算 GSVA。

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：从样本名推断 CON/EXP；已有分组列时不调用。
infer_group <- function(samples) {
  upper <- toupper(as.character(samples))
  group <- rep(NA_character_, length(upper))
  group[grepl("(^|[_\\.-])CON(TROL)?([_\\.-]|$)", upper)] <- "CON"
  group[grepl("(^|[_\\.-])EXP(ERIMENT)?([_\\.-]|$)", upper)] <- "EXP"
  group
}

# 功能：把 ggplot 同时保存为 PDF 和 600 dpi PNG。
save_plot <- function(plot, stem, width = 10, height = 7) {
  plot <- plot + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
}

# 功能：读取 GMT 并返回命名基因集列表。
read_gene_sets <- function(gmt_file) {
  collection <- GSEABase::getGmt(gmt_file)
  stats::setNames(GSEABase::geneIds(collection), GSEABase::setName(collection))
}

# 功能：从 Seurat v5 提取非中心化表达层。
seurat_expression <- function(obj, assay, layer) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay)
  SeuratObject::LayerData(obj, assay = assay, layer = layer)
}

# 功能：用稀疏设计矩阵计算每个样本×细胞类型的平均表达。
aggregate_expression <- function(expr, metadata, sample_col, group_col, celltype_col) {
  group_id <- interaction(
    metadata[[sample_col]], metadata[[celltype_col]],
    drop = TRUE, lex.order = TRUE, sep = "|||"
  )
  group_index <- as.integer(group_id)
  design <- Matrix::sparseMatrix(
    i = seq_along(group_index), j = group_index, x = 1,
    dims = c(length(group_index), nlevels(group_id))
  )
  counts <- as.numeric(table(group_id))
  averaged <- expr %*% design
  averaged <- t(t(averaged) / counts)
  colnames(averaged) <- levels(group_id)

  first_cell <- match(levels(group_id), as.character(group_id))
  sample_table <- data.frame(
    aggregate_id = levels(group_id),
    sample = as.character(metadata[[sample_col]][first_cell]),
    group = as.character(metadata[[group_col]][first_cell]),
    cell_type = as.character(metadata[[celltype_col]][first_cell]),
    n_cells = counts,
    stringsAsFactors = FALSE
  )
  list(expression = averaged, metadata = sample_table)
}

# 功能：使用当前 GSVA 参数对象 API 计算样本级通路分数。
run_gsva <- function(expression, gene_sets, min_size, max_size, workers) {
  parameter <- GSVA::gsvaParam(
    as.matrix(expression), gene_sets,
    kcdf = "Gaussian", minSize = min_size, maxSize = max_size
  )
  GSVA::gsva(
    parameter,
    verbose = TRUE,
    BPPARAM = BiocParallel::SnowParam(workers = workers, type = "SOCK")
  )
}

# 功能：在每个细胞类型和通路内以生物学样本为单位比较 EXP 与 CON。
compare_scores <- function(score_long) {
  split_data <- split(score_long, interaction(score_long$cell_type, score_long$pathway, drop = TRUE))
  records <- lapply(split_data, function(table) {
    con <- table$score[table$group == "CON"]
    exp <- table$score[table$group == "EXP"]
    pvalue <- if (length(con) >= 2 && length(exp) >= 2) stats::wilcox.test(exp, con)$p.value else NA_real_
    data.frame(
      cell_type = table$cell_type[1], pathway = table$pathway[1],
      n_CON = length(con), n_EXP = length(exp),
      mean_CON = mean(con), mean_EXP = mean(exp),
      delta_EXP_minus_CON = mean(exp) - mean(con), pvalue = pvalue
    )
  })
  result <- do.call(rbind, records)
  result$p_adjust_BH <- stats::p.adjust(result$pvalue, method = "BH")
  result[order(result$p_adjust_BH, result$pvalue), ]
}

# 功能：绘制最显著通路的样本级箱线图，保留每个样本点。
plot_score_boxplots <- function(score_long, statistics, output_dir, top_n) {
  selected <- head(unique(statistics$pathway[order(statistics$p_adjust_BH, statistics$pvalue)]), top_n)
  table <- score_long[score_long$pathway %in% selected, ]
  table$group <- factor(table$group, levels = c("CON", "EXP"))
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = cell_type, y = score, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.32, position = ggplot2::position_dodge(0.78)) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78),
      size = 1.4
    ) +
    ggplot2::facet_wrap(~pathway, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(CON = "#E64B35", EXP = "#4DBBD5"), drop = FALSE) +
    ggplot2::scale_color_manual(values = c(CON = "#E64B35", EXP = "#4DBBD5"), drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Sample-level GSVA score") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1), legend.position = "top")
  save_plot(plot, file.path(output_dir, "sample_level_GSVA_boxplot"), 14, 10)
}

# 功能：绘制 EXP-CON 通路差值热图；红色表示 EXP 高，蓝色表示 CON 高。
plot_delta_heatmap <- function(statistics, output_dir) {
  matrix <- xtabs(delta_EXP_minus_CON ~ cell_type + pathway, statistics)
  colors <- grDevices::colorRampPalette(c("#3C5488", "white", "#E64B35"))(101)
  grDevices::pdf(file.path(output_dir, "GSVA_delta_EXP_minus_CON_heatmap.pdf"), width = 11, height = 8)
  pheatmap::pheatmap(matrix, color = colors, border_color = NA)
  grDevices::dev.off()
  grDevices::png(file.path(output_dir, "GSVA_delta_EXP_minus_CON_heatmap.png"), width = 11, height = 8, units = "in", res = 600)
  pheatmap::pheatmap(matrix, color = colors, border_color = NA)
  grDevices::dev.off()
}

# 功能：串联样本级表达汇总、GSVA、统计、绘图和中文说明。
run_workflow <- function(
  input_rds, gmt_file, output_dir,
  assay, layer, sample_col, group_col, celltype_col,
  min_size, max_size, workers, top_n
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(input_rds)
  metadata <- obj[[]]
  if (!group_col %in% colnames(metadata)) metadata[[group_col]] <- infer_group(metadata[[sample_col]])
  expr <- seurat_expression(obj, assay, layer)
  metadata <- metadata[colnames(expr), , drop = FALSE]
  aggregated <- aggregate_expression(expr, metadata, sample_col, group_col, celltype_col)
  scores <- run_gsva(aggregated$expression, read_gene_sets(gmt_file), min_size, max_size, workers)

  score_long <- merge(
    aggregated$metadata,
    data.frame(
      aggregate_id = rep(colnames(scores), each = nrow(scores)),
      pathway = rep(rownames(scores), times = ncol(scores)),
      score = as.vector(scores), stringsAsFactors = FALSE
    ),
    by = "aggregate_id"
  )
  write.csv(score_long, file.path(output_dir, "GSVA_scores_sample_celltype.csv"), row.names = FALSE)
  statistics <- compare_scores(score_long)
  write.csv(statistics, file.path(output_dir, "GSVA_CON_EXP_statistics.csv"), row.names = FALSE)
  plot_score_boxplots(score_long, statistics, output_dir, top_n)
  plot_delta_heatmap(statistics, output_dir)

  writeLines(
    c(
      "# GSVA 结果说明", "",
      paste0("- 输入 RDS：", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
      paste0("- GMT：", normalizePath(gmt_file, winslash = "/", mustWork = FALSE)),
      paste0("- 样本/分组/细胞类型列：", sample_col, "/", group_col, "/", celltype_col, "。"),
      "- 先在样本×细胞类型层面汇总表达，再计算 GSVA；不把单细胞当作生物学重复。",
      "- delta_EXP_minus_CON > 0 表示 EXP 分数较高；组间 P 值为样本级 Wilcoxon 检验并做 BH 校正。",
      "- GSVA 是相对转录富集分数，不直接代表蛋白活性或代谢通量。"
    ),
    file.path(output_dir, "README_results_CN.md")
  )
}

input_rds <- arg_value("--input")
gmt_file <- arg_value("--gmt")
output_dir <- arg_value("--output")
if (any(vapply(list(input_rds, gmt_file, output_dir), is.null, logical(1)))) {
  stop("用法：Rscript gsva_workflow.R --input annotated.rds --gmt signatures.gmt --output results/gsva")
}
run_workflow(
  input_rds, gmt_file, output_dir,
  assay = arg_value("--assay", "RNA"), layer = arg_value("--layer", "data"),
  sample_col = arg_value("--sample-col", "sample"), group_col = arg_value("--group-col", "group"),
  celltype_col = arg_value("--celltype-col", "cell_type"),
  min_size = as.integer(arg_value("--min-size", "5")), max_size = as.integer(arg_value("--max-size", "500")),
  workers = as.integer(arg_value("--workers", "4")), top_n = as.integer(arg_value("--top-n", "12"))
)

}
if (sys.nframe() == 0L) main()
