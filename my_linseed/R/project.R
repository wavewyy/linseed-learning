# 第 5 步：将基因的“跨样本表达模式”投影到低维空间。
# 输入 X_norm：基因×样本，使用第 4 步筛选后的归一化矩阵。
# 三类组分时：用 2 个中心化主方向描述平面，再加 1 个均值正交方向。
# 因此坐标矩阵是 3×基因数；每一列是一个基因，不是一个细胞。
# 当前实现对应官方 spherize=FALSE，使用 base R 的 svd()。

project_genes <- function(X_norm, config) {
  if (!is.matrix(X_norm) || !is.numeric(X_norm) ||
      any(dim(X_norm) == 0L) || any(!is.finite(X_norm)) || any(X_norm < 0)) {
    stop("X_norm 必须是非空、非负且有限的数值矩阵。")
  }
  if (is.null(rownames(X_norm)) || is.null(colnames(X_norm)) ||
      anyDuplicated(rownames(X_norm)) || anyDuplicated(colnames(X_norm))) {
    stop("投影输入需要唯一的基因名和样本名。")
  }
  if (any(abs(rowSums(X_norm) - 1) > 1e-8)) stop("请先按基因归一化。")
  k <- config$cell_types
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 2L || k != floor(k)) stop("cell_types 必须是至少为 2 的整数。")
  if (nrow(X_norm) < k || ncol(X_norm) < k) {
    stop("基因数和样本数都必须不少于类型数。")
  }

  # 1. 转置：让基因成为点（列），样本成为原始空间的坐标轴（行）。
  Y <- t(X_norm)

  # 2. 对所有基因点求平均位置，并将点云移到原点。
  #    此处的均值沿“基因”方向计算，不是再次对每个基因求平均。
  mean_vector <- rowMeans(Y)
  centered <- sweep(Y, 1, mean_vector, "-")

  # 3. SVD 找到点云变化最大的 k-1 个正交方向。
  decomposition <- svd(centered, nu = k - 1L, nv = 0L)
  singular_values <- decomposition$d
  rank_tolerance <- max(dim(centered)) * .Machine$double.eps * max(singular_values)
  numerical_rank <- sum(singular_values > rank_tolerance)
  if (numerical_rank < k - 1L) {
    stop("中心化数据秩不足，无法得到要求的 k-1 个方向；请检查数据或类型数。")
  }
  directions <- decomposition$u[, seq_len(k - 1L), drop = FALSE]

  # 4. 取均值中垂直于这些主方向的分量，构造第 k 个基向量。
  #    它保留平面相对原点的位置，供后续角点搜索和还原使用。
  mean_orthogonal <- mean_vector - as.vector(directions %*% crossprod(directions, mean_vector))
  orthogonal_length <- sqrt(sum(mean_orthogonal^2))
  if (!is.finite(orthogonal_length) ||
      orthogonal_length <= .Machine$double.eps * max(1, sqrt(sum(mean_vector^2)))) {
    stop("无法构造均值正交方向，请检查归一化输入。")
  }
  basis <- cbind(directions, mean_orthogonal / orthogonal_length)
  axis_names <- c(paste0("PC", seq_len(k - 1L)), "mean_orthogonal")
  dimnames(basis) <- list(colnames(X_norm), axis_names)

  # 5. 先将点投影到 k-1 维仿射平面，再用 basis 表示坐标。
  #    与官方的 proj %*% centered + mean 等价，但不创建大投影矩阵。
  plane_points <- sweep(directions %*% crossprod(directions, centered), 1, mean_vector, "+")
  coordinates <- crossprod(basis, plane_points)
  dimnames(coordinates) <- list(axis_names, rownames(X_norm))

  # 6. 保存质量统计：以后可以检查选择 k 是否损失了过多信息。
  variance <- singular_values^2
  variance_fraction <- variance / sum(variance)
  relative_error <- sqrt(sum((Y - plane_points)^2) / sum(Y^2))
  list(
    coordinates = coordinates,             # k×基因；前 k-1 行可用于画图
    basis = basis,                         # 样本×k；basis %*% coordinates 可还原投影点
    mean = setNames(mean_vector, colnames(X_norm)),
    singular_values = singular_values,
    variance_fraction = variance_fraction, # 中心化点云各方向解释的方差比例
    retained_variance = sum(variance_fraction[seq_len(k - 1L)]),
    relative_error = relative_error,
    cell_types = k
  )
}
