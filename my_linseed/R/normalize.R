# 按基因归一化：每行是一个基因，每列是一个样本。
# 公式：X_norm[i, j] = X[i, j] / sum(X[i, ])
# 例如 c(10, 20, 30) -> c(1/6, 2/6, 3/6)。
# 输入必须是预处理后的非负线性表达量；输出每行之和约等于 1。
# 这一步不改变矩阵大小、行列顺序或名称，也不是按样本做 CPM/TPM。

normalize_genes <- function(X, config = list()) {
  # config 保留给 main.R 统一调用；本步骤没有需要设置的算法参数。
  if (!is.matrix(X) || !is.numeric(X) || any(dim(X) == 0L)) {
    stop("X 必须是非空数值矩阵，行是基因，列是样本。")
  }
  if (any(!is.finite(X)) || any(X < 0)) {
    stop("归一化输入不能含 NA、无穷值或负数。")
  }

  # 1. 求每个基因在所有样本中的表达总和。
  gene_totals <- rowSums(X)
  if (any(!is.finite(gene_totals))) stop("基因表达总和溢出，请检查输入数值。")
  if (any(gene_totals <= 0)) {
    stop("存在全零基因，无法除以 0；请先执行预处理去掉全零行。")
  }

  # 2. sweep(..., MARGIN = 1) 表示逐行计算。
  #    每一行都除以该行对应的 gene_totals。
  X_norm <- sweep(X, MARGIN = 1, STATS = gene_totals, FUN = "/")

  # 3. 浮点运算允许微小误差，不能用 == 1 做严格比较。
  if (any(abs(rowSums(X_norm) - 1) > 1e-10)) {
    stop("归一化校验失败：存在行总和未接近 1。")
  }
  X_norm
}
