# 第 7 步：利用角点做解卷积。
# 对照 LinSeed 的 pureDsa（MIT，见 THIRD_PARTY_LICENSES/LinSeed-MIT.txt）。
# NMF 只提供非负最小二乘求解器；下面自己组织两次求解。
deconvolve_expression <- function(X, X_norm, projection, corners, config) {
  if (!requireNamespace("NMF", quietly = TRUE)) stop("需要安装 NMF 包。")
  valid_matrix <- function(A) {
    is.matrix(A) && is.numeric(A) && length(A) > 0L &&
      all(is.finite(A)) && all(A >= 0)
  }
  if (!valid_matrix(X) || !valid_matrix(X_norm) ||
      !identical(dim(X), dim(X_norm)) ||
      !identical(dimnames(X), dimnames(X_norm))) {
    stop("X 与 X_norm 必须是同维度、同基因和样本顺序的非负有限矩阵。")
  }
  if (is.null(rownames(X)) || is.null(colnames(X)) ||
      anyDuplicated(rownames(X)) || anyDuplicated(colnames(X))) {
    stop("表达矩阵需要唯一的基因名和样本名。")
  }
  totals <- rowSums(X)
  if (any(!is.finite(totals)) || any(totals <= 0) ||
      max(abs(X_norm - sweep(X, 1L, totals, "/"))) > 1e-8) {
    stop("X_norm 应为 X 除以每个基因的行总和。")
  }
  E <- corners$endpoints # 已还原到样本空间：样本 × 类型
  k <- config$cell_types
  if (length(k) != 1L || !is.finite(k) || k < 2 || k != as.integer(k)) {
    stop("cell_types 必须是至少为 2 的整数。")
  }
  if (!valid_matrix(E) || !identical(dim(E), c(ncol(X), as.integer(k))) ||
      !identical(rownames(E), colnames(X))) {
    stop("角点必须为非负的样本×类型矩阵，样本顺序须与表达矩阵一致。")
  }
  if (qr(E)$rank < k) stop("角点不独立，无法区分指定的细胞类型。")
  # projection 留在接口中；角点已在第 6 步还原，因此这里无需再次投影。

  # 1. 找到每个角点的缩放系数 a，使 E %*% a 尽量接近全 1 向量。
  #    a >= 0；这等价于让每个样本的各类型比例之和尽量接近 1。
  a <- as.numeric(NMF::fcnnls(E, matrix(1, ncol(X), 1L))$x)
  P <- diag(a, nrow = k) %*% t(E) # 类型 × 样本

  # 2. 固定比例 P，为每个基因求非负表达特征 S，使 X_norm ≈ S %*% P。
  S <- t(NMF::.fcnnls(t(P), t(X_norm), pseudo = TRUE)$coef)
  components <- colnames(E)
  if (is.null(components)) components <- paste0("Component", seq_len(k))
  dimnames(P) <- list(components, colnames(X))
  dimnames(S) <- list(rownames(X), components)
  names(a) <- components
  if (!valid_matrix(P) || !valid_matrix(S)) stop("非负最小二乘返回无效结果。")
  if (any(a == 0)) warning("有类型的缩放系数为 0；请检查角点及类型数。")

  # 3. 恢复每个基因归一化前的尺度；P 保持相同。
  #    因为每个基因分别求最小二乘，行乘回总量即可恢复线性尺度的 S。
  S_linear <- sweep(S, 1L, totals, "*")
  reconstructed <- S %*% P
  reconstructed_linear <- S_linear %*% P
  if (any(!is.finite(S_linear)) || any(!is.finite(reconstructed_linear))) {
    stop("线性尺度结果溢出，请检查输入数值大小。")
  }
  diagnostics <- data.frame(
    scale = c("gene_normalized", "linear"),
    frobenius_error = c(norm(X_norm - reconstructed, "F"),
                        norm(X - reconstructed_linear, "F"))
  )
  diagnostics$relative_error <- diagnostics$frobenius_error /
    c(norm(X_norm, "F"), norm(X, "F"))
  list(signatures = S, proportions = P, signatures_linear = S_linear,
       reconstructed = reconstructed, reconstructed_linear = reconstructed_linear,
       scaling = a, proportion_sums = colSums(P), diagnostics = diagnostics)
}
