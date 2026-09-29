# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# GO/KEGG ORA 与 preranked GSEA 独立入口
# 输入：主流程导出的差异表达 CSV；默认 logFC > 0 表示 EXP 高于 CON。

# 功能：读取 --name value 命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：把 ggplot 同时保存为矢量 PDF 和 600 dpi PNG。
save_plot <- function(plot, stem, width = 9, height = 7) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
}

# 功能：按物种选择 OrgDb 和 KEGG organism 代码。
species_resources <- function(species) {
  if (species == "human") {
    list(orgdb = org.Hs.eg.db::org.Hs.eg.db, kegg = "hsa")
  } else {
    list(orgdb = org.Mm.eg.db::org.Mm.eg.db, kegg = "mmu")
  }
}

# 功能：将基因 symbol 映射为唯一 Entrez ID，并保留可审计映射表。
map_symbols <- function(symbols, orgdb) {
  mapping <- AnnotationDbi::select(
    orgdb, keys = unique(as.character(symbols)),
    keytype = "SYMBOL", columns = "ENTREZID"
  )
  unique(mapping[!is.na(mapping$ENTREZID) & mapping$ENTREZID != "", c("SYMBOL", "ENTREZID")])
}

# 功能：分别对上调和下调基因运行 GO BP 与 KEGG 过度富集。
run_ora <- function(deg, mapping, resources, gene_col, logfc_col, padj_col, padj_cutoff, logfc_cutoff) {
  universe <- unique(mapping$ENTREZID[mapping$SYMBOL %in% deg[[gene_col]]])
  significant <- deg[!is.na(deg[[padj_col]]) & deg[[padj_col]] < padj_cutoff & abs(deg[[logfc_col]]) > logfc_cutoff, ]
  directions <- list(
    Up = significant[[gene_col]][significant[[logfc_col]] > 0],
    Down = significant[[gene_col]][significant[[logfc_col]] < 0]
  )
  lapply(directions, function(symbols) {
    genes <- unique(mapping$ENTREZID[mapping$SYMBOL %in% symbols])
    if (!length(genes)) return(list(GO = NULL, KEGG = NULL))
    list(
      GO = clusterProfiler::enrichGO(
        gene = genes, universe = universe, OrgDb = resources$orgdb,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      KEGG = clusterProfiler::enrichKEGG(
        gene = genes, universe = universe, organism = resources$kegg,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
      )
    )
  })
}

# 功能：保存每个方向的完整 ORA 表，不用非显著通路填满 Top N。
save_ora_tables <- function(results, output_dir) {
  for (direction in names(results)) {
    for (database in names(results[[direction]])) {
      result <- results[[direction]][[database]]
      table <- if (is.null(result)) data.frame() else as.data.frame(result)
      write.csv(table, file.path(output_dir, paste0(database, "_", direction, "_all.csv")), row.names = FALSE)
    }
  }
}

# 功能：绘制 Up 向右、Down 向左的 Top 10 富集条形图。
plot_diverging_ora <- function(results, database, output_dir, top_n = 10L) {
  tables <- lapply(c("Up", "Down"), function(direction) {
    result <- results[[direction]][[database]]
    if (is.null(result)) return(data.frame())
    table <- as.data.frame(result)
    table <- table[!is.na(table$p.adjust) & table$p.adjust < 0.05, ]
    table <- head(table[order(table$p.adjust), ], top_n)
    if (!nrow(table)) return(data.frame())
    table$direction <- direction
    table$signed_log10_padj <- -log10(pmax(table$p.adjust, .Machine$double.xmin)) * ifelse(direction == "Up", 1, -1)
    table
  })
  table <- do.call(rbind, tables)
  write.csv(table, file.path(output_dir, paste0(database, "_Up_Down_plot_data.csv")), row.names = FALSE)
  if (!nrow(table)) return(invisible(NULL))
  table$Description <- factor(table$Description, levels = rev(unique(table$Description)))
  plot <- ggplot2::ggplot(table, ggplot2::aes(Description, signed_log10_padj, fill = direction)) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(Up = "#E64B35", Down = "#3C5488")) +
    ggplot2::labs(x = NULL, y = "Signed -log10(adjusted P)", fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(color = "black"))
  save_plot(plot, file.path(output_dir, paste0(database, "_Up_Down_enrichment")), 10, 8)
}

# 功能：建立包含全部可检验基因的确定性 preranked 列表。
build_ranked_genes <- function(deg, mapping, gene_col, logfc_col, pvalue_col, stat_col) {
  table <- deg[!is.na(deg[[gene_col]]) & !is.na(deg[[logfc_col]]), ]
  pvalue <- pmax(table[[pvalue_col]], .Machine$double.xmin)
  table$rank_metric <- if (!is.null(stat_col) && stat_col %in% colnames(table)) {
    table[[stat_col]]
  } else {
    sign(table[[logfc_col]]) * stats::qnorm(pvalue / 2, lower.tail = FALSE)
  }
  table <- table[order(pvalue, -abs(table[[logfc_col]]), as.character(table[[gene_col]])), ]
  table <- table[!duplicated(table[[gene_col]]), ]
  ranked <- merge(table, mapping, by.x = gene_col, by.y = "SYMBOL")
  ranked <- ranked[order(-ranked$rank_metric, ranked$ENTREZID), ]
  ranked <- ranked[!duplicated(ranked$ENTREZID), ]
  # 极小的确定性偏移只打破完全相同的统计值，不改变实际排序方向。
  ranked$rank_metric <- ranked$rank_metric + seq_len(nrow(ranked)) * 1e-12
  values <- ranked$rank_metric
  names(values) <- ranked$ENTREZID
  list(values = sort(values, decreasing = TRUE), table = ranked)
}

# 功能：运行 GO BP、KEGG 及可选自定义 GMT 的 preranked GSEA。
run_gsea <- function(ranks, resources, gsea_gmt = NULL) {
  results <- list(
    GO = clusterProfiler::gseGO(
      geneList = ranks, OrgDb = resources$orgdb, keyType = "ENTREZID",
      ont = "BP", minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, by = "fgsea"
    ),
    KEGG = clusterProfiler::gseKEGG(
      geneList = ranks, organism = resources$kegg,
      minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, by = "fgsea"
    )
  )
  if (!is.null(gsea_gmt)) {
    term2gene <- clusterProfiler::read.gmt(gsea_gmt)
    results$GMT <- clusterProfiler::GSEA(
      ranks, TERM2GENE = term2gene,
      minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, by = "fgsea"
    )
  }
  results
}

# 功能：保存完整 GSEA 表，并为可选目标通路绘制标准 enrichment curve。
save_gsea_results <- function(results, output_dir, target_pathway = NULL) {
  for (database in names(results)) {
    table <- as.data.frame(results[[database]])
    write.csv(table, file.path(output_dir, paste0("GSEA_", database, "_all.csv")), row.names = FALSE)
    if (!is.null(target_pathway) && nrow(table)) {
      index <- which(table$ID == target_pathway | table$Description == target_pathway)
      if (length(index)) {
        plot <- enrichplot::gseaplot2(
          results[[database]], geneSetID = index[1],
          title = table$Description[index[1]], pvalue_table = TRUE
        )
        save_plot(plot, file.path(output_dir, paste0("GSEA_", database, "_target")), 8.2, 6.7)
      }
    }
  }
}

# 功能：串联映射、Up/Down ORA、全基因 GSEA、绘图和结果说明。
run_workflow <- function(
  input_csv, output_dir, species,
  gene_col, logfc_col, pvalue_col, padj_col, stat_col,
  padj_cutoff, logfc_cutoff, gsea_gmt, target_pathway
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  deg <- read.csv(input_csv, check.names = FALSE)
  resources <- species_resources(species)
  mapping <- map_symbols(deg[[gene_col]], resources$orgdb)
  write.csv(mapping, file.path(output_dir, "gene_symbol_to_entrez_mapping.csv"), row.names = FALSE)

  ora <- run_ora(deg, mapping, resources, gene_col, logfc_col, padj_col, padj_cutoff, logfc_cutoff)
  save_ora_tables(ora, output_dir)
  plot_diverging_ora(ora, "GO", output_dir)
  plot_diverging_ora(ora, "KEGG", output_dir)

  ranked <- build_ranked_genes(deg, mapping, gene_col, logfc_col, pvalue_col, stat_col)
  write.csv(ranked$table, file.path(output_dir, "GSEA_ranked_gene_list.csv"), row.names = FALSE)
  gsea <- run_gsea(ranked$values, resources, gsea_gmt)
  save_gsea_results(gsea, output_dir, target_pathway)
  if (!is.null(gsea_gmt)) file.copy(gsea_gmt, file.path(output_dir, "GSEA_gene_sets_snapshot.gmt"), overwrite = TRUE)

  writeLines(
    c(
      "# 富集分析结果说明", "",
      paste0("- 输入差异表：", normalizePath(input_csv, winslash = "/", mustWork = FALSE)),
      paste0("- 物种：", species, "；比较方向：正 logFC = EXP 高于 CON。"),
      paste0("- ORA 阈值：adjusted P < ", padj_cutoff, " 且 |logFC| > ", logfc_cutoff, "；Up/Down 分开富集。"),
      "- GO/KEGG ORA 使用全部可映射、可检验基因作为背景。",
      "- GSEA 使用全部可检验基因排序，不先按显著性筛选；正 NES 表示 EXP 端富集，负 NES 表示 CON 端富集。",
      "- 所有富集结果采用 BH 校正；转录富集不能直接解释为蛋白活性或代谢通量。"
    ),
    file.path(output_dir, "README_results_CN.md")
  )
}

input_csv <- arg_value("--input")
output_dir <- arg_value("--output")
if (is.null(input_csv) || is.null(output_dir)) {
  stop("用法：Rscript enrichment_workflow.R --input differential_expression.csv --output results/enrichment --species human")
}
run_workflow(
  input_csv, output_dir, species = arg_value("--species", "human"),
  gene_col = arg_value("--gene-col", "gene"), logfc_col = arg_value("--logfc-col", "avg_log2FC"),
  pvalue_col = arg_value("--pvalue-col", "p_val"), padj_col = arg_value("--padj-col", "p_val_adj"),
  stat_col = arg_value("--stat-col", NULL),
  padj_cutoff = as.numeric(arg_value("--padj-cutoff", "0.05")),
  logfc_cutoff = as.numeric(arg_value("--logfc-cutoff", "0.25")),
  gsea_gmt = arg_value("--gsea-gmt", NULL), target_pathway = arg_value("--target-pathway", NULL)
)

}
if (sys.nframe() == 0L) main()
