# 第 6 步：VCA 初始化 -> SISAL 优化 -> 对齐候选角点 -> 比较重建误差。
# 依据 ctlab/LinSeed 的 vca.R、sisal.R、smartSearchCorners 改写。
# Copyright (c) 2019 Konstantin Zaitsev and other contributors.
# 原项目为 MIT License，完整许可见本项目 THIRD_PARTY_LICENSES/LinSeed-MIT.txt。
# 本文件不调用 linseed 的算法；corpcor/NMF 仅提供伪逆与非负最小二乘。

# VCA：每次向一个与已有方向正交的随机方向投影，挑最远的点。
corners_vca <- function(Y, k) {
  n <- ncol(Y)
  U <- svd(tcrossprod(Y) / n, nu = k, nv = 0)$u
  xp <- crossprod(U, Y)
  reconstructed <- U %*% xp
  u <- rowMeans(xp)
  denominator <- colSums(xp * u)
  if (any(abs(denominator) < .Machine$double.eps)) stop("VCA 投影分母接近 0。")
  y <- sweep(xp, 2, denominator, "/")
  A <- matrix(0, k, k)
  A[k, 1] <- 1
  indices <- integer(k)
  for (i in seq_len(k)) {
    w <- matrix(runif(k), ncol = 1)
    f <- w - A %*% corpcor::pseudoinverse(A) %*% w
    length_f <- sqrt(sum(f^2))
    if (length_f < .Machine$double.eps) stop("VCA 初始化方向退化。")
    f <- f / length_f
    indices[i] <- which.max(abs(crossprod(f, y)))
    A[, i] <- y[, indices[i]]
  }
  if (anyDuplicated(indices)) stop("VCA 选择了重复的初始点。")
  reconstructed[, indices, drop = FALSE]
}

# SISAL：在低维空间中优化单纯形，平衡体积和落在单纯形外的数据惩罚。
# Q 是角点矩阵 M 的逆，Q %*% Y 给出每个点相对单纯形的系数。
corners_sisal <- function(projection, tau, iterations, initial = NULL) {
  Y <- projection$coordinates
  basis <- projection$basis
  k <- nrow(Y)
  mu <- k * 1000 / ncol(Y)
  hinge <- function(Z) pmax(-Z, 0)
  soft_negative <- function(Z, threshold) {
    z <- pmax(abs(Z + threshold / 2) - threshold / 2, 0)
    z / (z + threshold / 2) * (Z + threshold / 2)
  }
  objective <- function(Q) {
    log_volume <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
    -log_volume + tau * sum(hinge(Q %*% Y))
  }

  # 1. 初始单纯形放大，避免一开始过度切掉数据。
  if (is.null(initial)) {
    M <- corners_vca(Y, k)
    M <- M + k * sweep(M, 1, rowMeans(M), "-")
  } else M <- initial
  Q <- solve(M)

  # 2. 预先计算增广拉格朗日迭代所需矩阵。
  # R 按列展开矩阵，因此这里的 Kronecker 积次序与官方一致。
  H <- 1e-6 * diag(k^2)
  AAT <- kronecker(tcrossprod(Y), diag(k))
  B <- kronecker(diag(k), matrix(1, 1, k))
  qm <- matrix(rowSums(solve(tcrossprod(Y), Y)), ncol = 1)
  inverse_F <- solve(H + mu * AAT)
  G0 <- inverse_F %*% t(B) %*% solve(B %*% inverse_F %*% t(B))
  qm_aux <- G0 %*% qm
  G <- inverse_F - G0 %*% B %*% inverse_F
  Z <- Q %*% Y
  dual <- 0 * Z
  history <- numeric(iterations)

  # 3. 外层更新体积项，内层更新系数和负值惩罚。
  for (iteration in seq_len(iterations)) {
    gradient <- matrix(-t(solve(Q)), ncol = 1)
    q0 <- matrix(Q, ncol = 1)
    Q0 <- Q
    b_aux <- H %*% q0 - gradient
    inner_iterations <- if (iteration == iterations) 100L else 4L
    quadratic <- function(q, Q) {
      delta <- q - q0
      as.numeric(crossprod(delta, gradient) + 0.5 * crossprod(delta, H %*% delta)) +
        tau * sum(hinge(Q %*% Y))
    }
    accepted <- FALSE
    # 官方使用无界 while；这里加上限，数值退化时明确报错。
    for (attempt in seq_len(200L)) {
      q <- matrix(Q, ncol = 1)
      old_objective <- objective(Q)
      old_quadratic <- quadratic(q, Q)
      for (inner in seq.int(2L, inner_iterations)) {
        rhs <- b_aux + mu * matrix((Z + dual) %*% t(Y), ncol = 1)
        q <- G %*% rhs + qm_aux
        Q <- matrix(q, k, k)
        Z <- soft_negative(Q %*% Y - dual, tau / mu)
        dual <- dual - (Q %*% Y - Z)
      }
      if (old_quadratic >= quadratic(q, Q)) {
        for (backtrack in seq_len(100L)) {
          value <- objective(Q)
          if (is.finite(value) && value <= old_objective) break
          Q <- (Q + Q0) / 2
        }
        value <- objective(Q)
        if (!is.finite(value) || value > old_objective + 1e-10) {
          stop("SISAL 回溯未收敛。")
        }
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("SISAL 内层未收敛，请检查数据或参数。")
    history[iteration] <- objective(Q)
  }

  # 4. 将角点还原到样本空间；负表达角点向数据均值收缩。
  endpoints <- basis %*% solve(Q)
  if (any(endpoints < 0)) {
    difference <- projection$mean - endpoints
    shifts <- vapply(seq_len(k), function(i) {
      negative <- endpoints[, i] < 0
      if (any(negative)) -min(endpoints[negative, i] / difference[negative, i]) + 1e-5 else 1e-5
    }, numeric(1))
    endpoints <- endpoints + difference %*% diag(shifts, nrow = k)
  }
  if (any(!is.finite(endpoints)) || any(endpoints < -1e-10)) {
    stop("SISAL 产生无效角点。")
  }
  list(endpoints = endpoints, objective = history)
}

# 用角点临时估计系数并计算重建误差，用于挑选角点，而非第7步正式输出。
corners_reconstruction_error <- function(endpoints, X_norm) {
  multiplier <- NMF::fcnnls(endpoints, matrix(1, nrow(endpoints), 1))$x
  H <- diag(as.numeric(multiplier), nrow = ncol(endpoints)) %*% t(endpoints)
  coefficients <- NMF::.fcnnls(t(H), t(X_norm), pseudo = TRUE)$coef
  norm(X_norm - t(coefficients) %*% H, "F")
}

find_corners <- function(projection, X_norm, config) {
  if (!requireNamespace("corpcor", quietly = TRUE) || !requireNamespace("NMF", quietly = TRUE)) {
    stop("角点搜索需要已有的 corpcor 和 NMF 数值求解包。")
  }
  k <- config$cell_types
  Y <- projection$coordinates
  if (!is.matrix(Y) || !is.matrix(X_norm) || any(!is.finite(Y)) ||
      !identical(dim(Y), c(k, nrow(X_norm))) ||
      !identical(colnames(Y), rownames(X_norm)) ||
      !identical(rownames(projection$basis), colnames(X_norm))) {
    stop("投影坐标、基因名、样本名或维度不匹配。")
  }
  taus <- config$corner_taus
  iterations <- config$sisal_iterations
  if (!is.numeric(taus) || !length(taus) || any(!is.finite(taus)) || any(taus <= 0)) {
    stop("corner_taus 必须是正数向量。")
  }
  if (!is.numeric(iterations) || length(iterations) != 1L || !is.finite(iterations) ||
      iterations < 1L || iterations != floor(iterations)) stop("sisal_iterations 必须为正整数。")

  # 1. 每个 tau 独立初始化并优化，得到多个候选单纯形。
  candidates <- vector("list", length(taus))
  diagnostics <- data.frame(tau = taus, success = FALSE, error = NA_real_, note = "")
  for (i in seq_along(taus)) {
    message(sprintf("SISAL 候选 %d/%d：tau=%g，迭代=%d", i, length(taus), taus[i], iterations))
    candidates[i] <- list(tryCatch(corners_sisal(projection, taus[i], iterations), error = function(e) {
      diagnostics$note[i] <<- conditionMessage(e)
      NULL
    }))
  }
  valid <- which(!vapply(candidates, is.null, logical(1)))
  if (!length(valid)) stop("所有 SISAL 候选都失败：", paste(diagnostics$note, collapse = "; "))

  # 2. 同一组分在不同候选中可能位于不同列；按均值中心化相关性对齐。
  reference <- sweep(candidates[[valid[1]]]$endpoints, 1, projection$mean, "-")
  for (i in valid) {
    centered <- sweep(candidates[[i]]$endpoints, 1, projection$mean, "-")
    correlations <- cor(reference, centered)
    if (any(!is.finite(correlations))) {
      candidates[i] <- list(NULL)
      diagnostics$note[i] <- "角点相关性退化"
      next
    }
    order <- apply(correlations, 1, which.max)
    # 不允许将两个组分都映射到同一角点。
    if (anyDuplicated(order)) {
      candidates[i] <- list(NULL)
      diagnostics$note[i] <- "角点对齐出现重复，跳过此候选"
      next
    }
    candidates[[i]]$endpoints <- candidates[[i]]$endpoints[, order, drop = FALSE]
    diagnostics$success[i] <- TRUE
  }
  valid <- which(diagnostics$success)
  if (!length(valid)) stop("没有可用且可对齐的角点候选。")
  score <- function(E) tryCatch({
    value <- corners_reconstruction_error(E, X_norm)
    if (is.finite(value)) value else Inf
  }, error = function(e) Inf)
  for (i in valid) diagnostics$error[i] <- score(candidates[[i]]$endpoints)
  valid <- valid[is.finite(diagnostics$error[valid])]
  if (!length(valid)) stop("所有候选角点都无法完成重建评分。")

  # 3. 仿官方逐个角点比较不同 tau 候选，采用误差更小的组合。
  # 与官方的累计停止计数不同，这里只接受严格改善，并设置搜索轮数上限。
  choice <- rep(valid[1], k)
  assemble <- function(index) do.call(cbind, lapply(seq_len(k), function(j) candidates[[index[j]]]$endpoints[, j]))
  current_error <- score(assemble(choice))
  converged <- FALSE
  trace <- list()
  for (sweep in seq_len(100L)) {
    improved <- FALSE
    for (j in seq_len(k)) {
      errors <- vapply(valid, function(i) {
        proposal <- choice
        proposal[j] <- i
        score(assemble(proposal))
      }, numeric(1))
      best <- which.min(errors)
      if (errors[best] < current_error - 1e-12 * max(1, current_error)) {
        choice[j] <- valid[best]
        current_error <- errors[best]
        improved <- TRUE
      }
    }
    trace[[sweep]] <- data.frame(sweep = sweep, error = current_error)
    if (!improved) { converged <- TRUE; break }
  }
  if (!converged) warning("角点组合搜索达到轮数上限，请检查搜索记录。")
  endpoints <- assemble(choice)
  dimnames(endpoints) <- list(rownames(projection$basis), paste0("Component_", seq_len(k)))
  coordinates <- crossprod(projection$basis, endpoints)
  distances <- vapply(seq_len(k), function(j) {
    sqrt(colSums(sweep(Y, 1, coordinates[, j], "-")^2))
  }, numeric(ncol(Y)))
  dimnames(distances) <- list(colnames(Y), colnames(endpoints))
  message(sprintf("已找到 %d 个角点；重建 Frobenius 误差 %.6g。", k, current_error))
  list(endpoints = endpoints, coordinates = coordinates, distances = distances,
       selected_tau = setNames(taus[choice], colnames(endpoints)),
       reconstruction_error = current_error, candidates = diagnostics,
       search_trace = do.call(rbind, trace), converged = converged,
       objective_history = lapply(candidates, function(x) if (is.null(x)) NULL else x$objective))
}
