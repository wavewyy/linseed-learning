# 数据预处理：输入为“基因 × 样本”表达矩阵，输出仍为表达矩阵。
# 使用 base R 实现，不调用官方 linseed 的算法函数。
# 当前处理基因级矩阵；探针到基因的注释合并、GEO 下载不在本模块内。

preprocess_expression <- function(X, config) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X) || any(dim(X) == 0L)) {
    stop("X 必须是非空数值矩阵：行是基因，列是样本。")
  }
  if (any(!is.finite(X)) || any(X < 0)) {
    stop("X 不能含 NA、无穷值或负数。")
  }
  if (is.null(rownames(X)) || is.null(colnames(X)) ||
      anyNA(rownames(X)) || anyNA(colnames(X)) ||
      any(!nzchar(rownames(X))) || any(!nzchar(colnames(X))) ||
      anyDuplicated(rownames(X)) || anyDuplicated(colnames(X))) {
    stop("X 必须有非空且唯一的基因名、样本名。")
  }
  top_genes <- config$top_genes
  if (!is.numeric(top_genes) || length(top_genes) != 1L ||
      !is.finite(top_genes) || top_genes < 1 || top_genes != floor(top_genes)) {
    stop("config$top_genes 必须是正整数。")
  }

  # 1. 选择样本。NULL 表示全部；10:42 表示选第 10 到 42 列。
  samples <- config$samples
  if (!is.null(samples)) {
    if (length(samples) == 0L || anyNA(samples) || anyDuplicated(samples)) {
      stop("samples 必须非空、不含缺失值且不重复。")
    }
    if (is.character(samples)) {
      if (any(!samples %in% colnames(X))) stop("samples 中有不存在的样本名。")
    } else if (is.numeric(samples)) {
      if (any(!is.finite(samples)) || any(samples != floor(samples)) ||
          any(samples < 1 | samples > ncol(X))) stop("samples 列编号超出范围。")
    } else stop("samples 应为列编号或样本名，或 NULL。")
    X <- X[, samples, drop = FALSE]
  }

  # 2. 去掉在所选样本中表达量全为 0 的基因。
  X <- X[rowSums(X) > 0, , drop = FALSE]
  if (nrow(X) == 0L) stop("去掉全零基因后，没有可用基因。")

  # 3. 确定输入尺度。auto 沿用官方的启发式规则：最大值>100
  #    或极差>100 时认为是线性数据，否则认为是 log2(x+1)。
  #    这只是猜测！真实数据应根据来源明确选择 linear 或 log2p1。
  scale <- config$input_scale
  if (is.null(scale)) scale <- "auto"
  if (length(scale) != 1L || is.na(scale) ||
      !scale %in% c("auto", "linear", "log2p1")) {
    stop("input_scale 只能是 auto、linear 或 log2p1。")
  }
  if (scale == "auto") {
    scale <- if (max(X) > 100 || diff(range(X)) > 100) "linear" else "log2p1"
    message("自动判断输入尺度：", scale, "（真实数据请明确指定）")
  }
  X_log <- if (scale == "linear") log2(X + 1) else X

  # 4. 在对数表达尺度上计算行总和，降序选取高表达基因。
  #    与直接在原始表达量上求总和排序不同。
  scores <- rowSums(X_log)
  n_keep <- min(top_genes, nrow(X_log))
  selected <- order(scores, decreasing = TRUE)[seq_len(n_keep)]
  X_log <- X_log[selected, , drop = FALSE]

  # 5. 返回线性表达量。线性输入直接取原值，避免往返转换的误差。
  X_out <- if (scale == "linear") X[selected, , drop = FALSE] else 2^X_log - 1
  if (any(!is.finite(X_out))) stop("转回线性尺度时溢出，请检查 input_scale。")

  # 6. 去除表达模式完全相同的重复行，保留排序靠前的那一行。
  #    去重后可能少于 top_genes，这是正常的。
  X_out <- X_out[!duplicated(X_out), , drop = FALSE]
  X_out
}
