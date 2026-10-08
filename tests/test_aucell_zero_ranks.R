# 合成排名验证零表达掩码；本检查不需要运行 AUCell 推断。
source("src/R/aucell_scores.R")
ranks <- matrix(c(1L, 2L, 3L, 3L, 1L, 2L), nrow = 3)
masked <- mask_aucell_zero_ranks(ranks, c(1, 2))
stopifnot(identical(colSums(!is.na(masked)), c(1, 2)),
          masked[1, 1] == 1, is.na(masked[2, 1]), is.na(masked[1, 2]))
stopifnot(inherits(try(mask_aucell_zero_ranks(ranks, c(4, 1)), silent = TRUE), "try-error"))
cat("AUCELL_ZERO_RANKS_OK\n")
