# 从现有原始表达计算给定基因集的 AUCell；不替换项目的基因集或统计单位。

# 将未检出基因的排名设为 NA；避免稀疏零值参与 AUCell 排名尾部。
mask_aucell_zero_ranks <- function(ranks, detected) {
  stopifnot(is.matrix(ranks), length(detected) == ncol(ranks),
            all(is.finite(detected)), all(detected >= 0 & detected <= nrow(ranks)))
  ranks[sweep(ranks, 2, detected, ">")] <- NA_integer_
  stopifnot(all(colSums(!is.na(ranks)) == detected))
  ranks
}

# 对一块细胞计算 AUCell，并排除未检出基因；依赖 AUCell 与 BiocParallel。
aucell_score_block <- function(expr, gene_sets, auc_max_rank) {
  rankings <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, splitByBlocks = TRUE,
    keepZeroesAsNA = FALSE, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  # AUCell 1.32's keepZeroesAsNA path mishandles sparse matrices/zero coordinates.
  # Positive counts always precede zero-count ties, so mask ranks beyond nDetected.
  ranks <- AUCell::getRanking(rankings)
  detected <- Matrix::colSums(expr > 0)
  ranks <- mask_aucell_zero_ranks(ranks, detected)
  SummarizedExperiment::assay(rankings, "ranking") <- ranks
  auc <- AUCell::AUCell_calcAUC(gene_sets, rankings, aucMaxRank = auc_max_rank,
                               normAUC = TRUE, nCores = 1, verbose = FALSE)
  AUCell::getAUC(auc)
}
