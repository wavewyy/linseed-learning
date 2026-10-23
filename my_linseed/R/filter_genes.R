# 第 4 步：基因关系 -> 网络权重 -> 随机检验 -> 基因筛选。
# 输入 X_norm：基因×样本，每行之和为 1。
# 返回 list：保留的基因名、逐基因统计、线性关系矩阵、Spearman 矩阵。
# 与官方一样：固定网络连接位置，对边权重进行有放回抽样。
# 这不是打乱表达矩阵，也没有做多重检验校正。

filter_genes <- function(X_norm, config) {
  if (!is.matrix(X_norm) || !is.numeric(X_norm) ||
      nrow(X_norm) < 2L || ncol(X_norm) < 2L ||
      any(!is.finite(X_norm)) || any(X_norm < 0)) {
    stop("X_norm 必须是至少 2×2 的非负有限数值矩阵。")
  }
  genes <- rownames(X_norm)
  if (is.null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("需要非空且唯一的基因名。")
  }
  if (any(abs(rowSums(X_norm) - 1) > 1e-8)) stop("请先进行按基因归一化。")
  B <- config$permutations
  cutoff <- config$pvalue_cutoff
  threshold <- config$spearman_threshold
  if (is.null(threshold)) threshold <- 0
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 1 || B != floor(B)) {
    stop("permutations 必须是正整数。")
  }
  if (length(cutoff) != 1L || !is.numeric(cutoff) || !is.finite(cutoff) ||
      cutoff <= 0 || cutoff >= 1) stop("pvalue_cutoff 必须在 0 与 1 之间。")
  if (length(threshold) != 1L || !is.numeric(threshold) || !is.finite(threshold) ||
      threshold < -1 || threshold >= 1) stop("spearman_threshold 必须在 [-1,1) 内。")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("请先安装 Matrix 包。")
  if (1 / (B + 1) >= cutoff) {
    warning("当前随机次数使最小可能 p 值 >= 筛选阈值，严格 p < 阈值无法保留基因。")
  }

  # 1. 计算官方 pairwiseR2 对应的对称线性关系。
  # distance[i,j] = sum((基因i - 基因j)^2)
  # TSS[i] = sum((基因i - mean(基因i))^2)
  # linearity[i,j] = 1 - distance[i,j]/2 * (1/TSS[i] + 1/TSS[j])
  n <- nrow(X_norm)
  squares <- rowSums(X_norm^2)
  distance <- outer(squares, squares, "+") - 2 * tcrossprod(X_norm)
  distance <- pmax(distance, 0)  # 消除浮点运算导致的微小负距离
  diag(distance) <- 0
  centered <- sweep(X_norm, 1, rowMeans(X_norm), "-")
  tss <- rowSums(centered^2)
  # 恒定基因没有跨样本变化，TSS 接近 0 时不参与网络。
  constant <- tss <= .Machine$double.eps * pmax(squares, .Machine$double.xmin)
  inverse_tss <- numeric(n)
  inverse_tss[!constant] <- 1 / tss[!constant]
  linearity <- 1 - distance * outer(inverse_tss, inverse_tss, "+") / 2
  linearity <- pmax(linearity, 0)  # 官方默认将负关系置为 0
  linearity[constant, ] <- 0
  linearity[, constant] <- 0
  dimnames(linearity) <- list(genes, genes)

  # 2. Spearman 相关性：比较基因在不同样本中的表达排序。
  spearman <- matrix(0, n, n, dimnames = list(genes, genes))
  active <- which(!constant)
  if (length(active)) {
    spearman[active, active] <- cor(t(X_norm[active, , drop = FALSE]), method = "spearman")
  }

  # 3. 只取上三角，避免重复计数；不包含自身连接。
  connected <- upper.tri(linearity) & linearity > 0 & spearman > threshold
  ij <- which(connected, arr.ind = TRUE)
  weights <- (linearity * spearman)[connected]
  # 默认 threshold=0，因此边权为正。阈值可调整以研究官方接口。
  degrees <- tabulate(c(ij[, 1], ij[, 2]), nbins = n)
  network_strength <- function(w) {
    network <- Matrix::sparseMatrix(i = ij[, 1], j = ij[, 2], x = w, dims = c(n, n))
    as.numeric(Matrix::rowSums(network + Matrix::t(network)))
  }
  observed <- network_strength(weights)

  # 4. 固定边的位置，从现有权重中有放回抽样。
  # p = (随机强度 >= 实际强度的次数 + 1) / (B + 1)。
  exceedances <- numeric(n)
  if (length(weights)) {
    for (b in seq_len(B)) {
      # sample.int 避免 sample(单个数值) 被误解成 1:该数值。
      random_weights <- weights[sample.int(length(weights), length(weights), replace = TRUE)]
      randomized <- network_strength(random_weights)
      exceedances <- exceedances + as.numeric(randomized >= observed)
    }
    pvalues <- (exceedances + 1) / (B + 1)
  } else {
    pvalues <- rep(1, n)
  }
  pvalues[constant] <- 1

  # 5. 与官方一致：严格使用 p < cutoff，不改成 <=。
  keep <- pvalues < cutoff & !constant
  statistics <- data.frame(
    gene = genes, degree = degrees, strength = observed,
    pvalue = pvalues, constant = constant, kept = keep,
    row.names = NULL
  )
  message(sprintf("关系网络：%d 个基因、%d 条边；筛选保留 %d 个基因。",
                  n, length(weights), sum(keep)))
  list(kept_genes = genes[keep], statistics = statistics,
       linearity = linearity, spearman = spearman)
}
