# my_linseed：项目入口
# 阅读顺序：参数配置 -> 算法接口 -> 主流程 -> 启动。
# 这里是算法框架，具体算法函数在 R文件夹中

# ---------- 1. 项目位置与运行环境 ----------
project_dir <- "D:/single cell项目/Linseed/my_linseed"
if (!dir.exists(project_dir)) stop("项目目录不存在，请修改 project_dir。")
setwd(project_dir)
.libPaths(c(file.path(dirname(project_dir), "library"), .libPaths()))

# ---------- 2. 参数配置：日常运行优先修改这里 ----------
config <- list(
  run_analysis = TRUE,      # FALSE：检查框架；TRUE：运行到 stop_after 指定步骤
  stop_after = "full",     # preprocess / normalize / filter / project / corners / full
  input_source = "demo",   # demo：quickstart 自带数据；csv：自己的 CSV
  samples = 10:42,          # demo 的混合样本；自己的数据一般改为 NULL
  input_scale = "auto",    # auto 仿官方猜测；真实数据明确设 linear 或 log2p1
  input_file = "data/bulk_expression.csv",
  cell_types = 3L,          # 第一版先手动指定类型数
  top_genes = 500L,         # 初学时使用较小的数据
  permutations = 100L,     # 随机检验次数
  pvalue_cutoff = 0.01,
  spearman_threshold = 0,  # 仅使用 Spearman 相关性 > 此值的正线性连接
  corner_taus = 2^seq(0, -8, -2), # 与 quickstart 相同：5 组候选
  sisal_iterations = 30L,         # 学习版；官方默认使用更多候选和迭代
  seed = 20261022L,
  results_dir = "results"
)

# ---------- 3. 算法模块：加载自己编写的函数 ----------
# 已实现全部 7 步；full 会运行到解卷积并保存结果。
# source() 负责读取你写的函数定义；不会安装 R 包。
modules <- c(
  preprocess_expression = "R/preprocess.R",
  normalize_genes = "R/normalize.R",
  filter_genes = "R/filter_genes.R",
  project_genes = "R/project.R",
  find_corners = "R/find_corners.R",
  deconvolve_expression = "R/deconvolve.R"
)

# 后续函数需要遵守的输入、输出约定：
# preprocess_expression(X, config) -> 处理后的表达矩阵，保留基因/样本名
# normalize_genes(X, config)      -> 同样行列顺序的归一化矩阵
# filter_genes(X_norm, config)    -> list(kept_genes, statistics, linearity, spearman)
# project_genes(X_norm, config)   -> list(coordinates, basis, mean, ...)
# find_corners(projection, X_norm, config) -> list(endpoints, coordinates, ...)
# deconvolve_expression(X, X_norm, projection, corners, config)
#   -> list(signatures = 基因×类型矩阵, proportions = 类型×样本矩阵)
# signatures 使用按基因归一化的尺度；signatures_linear 恢复线性表达尺度。

# ---------- 4. 通用辅助功能：不涉及 Linseed 算法 ----------
log_step <- function(text) {
  cat(format(Sys.time(), "[%H:%M:%S]"), text, "\n")
}

check_expression <- function(X) {
  if (!is.matrix(X) || !is.numeric(X)) stop("表达数据必须是数值矩阵。")
  if (nrow(X) < 3L || ncol(X) <= config$cell_types) {
    stop("需要至少 3 个基因，且样本数应大于指定的类型数。")
  }
  if (any(!is.finite(X)) || any(X < 0)) {
    stop("表达矩阵含缺失值、无穷值或负数，请检查输入。")
  }
  if (is.null(rownames(X)) || is.null(colnames(X)) ||
      anyDuplicated(rownames(X)) || anyDuplicated(colnames(X))) {
    stop("基因名与样本名必须存在且分别唯一。")
  }
  invisible(TRUE)
}

read_expression <- function(path) {
  if (!file.exists(path)) stop("找不到输入文件：", path)
  # CSV 格式：第一列是基因名；第一行是样本名；其余位置是表达量。
  # 例如：gene,sample_1,sample_2,...
  X <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  check_expression(X)
  X
}

load_modules <- function(required = names(modules)) {
  modules <- modules[required]
  missing_files <- modules[!file.exists(modules)]
  if (length(missing_files)) {
    stop("以下算法文件尚未实现：\n", paste(missing_files, collapse = "\n"),
         "\n请先完成这些文件，或将 run_analysis 设为 FALSE 检查框架。")
  }
  # 单独的环境避免误用当前交互会话里同名的旧函数。
  algorithm <- new.env(parent = globalenv())
  for (path in unique(modules)) sys.source(path, envir = algorithm)
  for (name in names(modules)) {
    if (!exists(name, envir = algorithm, inherits = FALSE) ||
        !is.function(algorithm[[name]])) stop("模块缺少函数：", name)
  }
  algorithm
}

# ---------- 5. 主流程：各算法完成后，从这里依次调用 ----------
run_pipeline <- function(config) {
  # 只加载本次运行需要的模块，不要求后续算法已经实现。
  required <- switch(config$stop_after,
    preprocess = "preprocess_expression",
    normalize = c("preprocess_expression", "normalize_genes"),
    filter = c("preprocess_expression", "normalize_genes", "filter_genes"),
    project = c("preprocess_expression", "normalize_genes", "filter_genes", "project_genes"),
    corners = c("preprocess_expression", "normalize_genes", "filter_genes", "project_genes", "find_corners"),
    full = names(modules)
  )
  algorithm <- load_modules(required)

  log_step("1/7 读取表达矩阵")
  if (config$input_source == "demo") {
    # 只从官方包取示例数据；预处理由我们自己的函数计算。
    demo <- new.env()
    data("datasetLiverBrainLung", package = "linseed", envir = demo)
    X <- demo$datasetLiverBrainLung
    check_expression(X)
  } else {
    X <- read_expression(config$input_file)
  }
  input_dim <- dim(X)

  log_step("2/7 预处理")
  X <- algorithm$preprocess_expression(X, config)
  check_expression(X)
  log_step(sprintf("预处理：%d 个基因 × %d 个样本 -> %d 个基因 × %d 个样本",
                   input_dim[1], input_dim[2], nrow(X), ncol(X)))
  if (config$stop_after == "preprocess") {
    return(list(preprocessed = X))
  }

  log_step("3/7 按基因归一化")
  X_norm <- algorithm$normalize_genes(X, config)
  check_expression(X_norm)
  if (!identical(dimnames(X), dimnames(X_norm))) {
    stop("归一化前后基因名、样本名及顺序必须保持一致。")
  }
  log_step(sprintf("归一化：%d 个基因 × %d 个样本；每行总和范围 %.12f–%.12f",
                   nrow(X_norm), ncol(X_norm),
                   min(rowSums(X_norm)), max(rowSums(X_norm))))
  if (config$stop_after == "normalize") {
    return(list(preprocessed = X, normalized = X_norm))
  }

  log_step("4/7 计算基因关系、显著性检验并筛选")
  filtering <- algorithm$filter_genes(X_norm, config)
  kept_genes <- filtering$kept_genes
  if (!is.character(kept_genes) ||
      anyDuplicated(kept_genes) || any(!kept_genes %in% rownames(X))) {
    stop("筛选模块应返回足够数量的、不重复且存在的基因名。")
  }
  X_selected <- X[kept_genes, , drop = FALSE]
  X_norm_selected <- X_norm[kept_genes, , drop = FALSE]
  if (config$stop_after == "filter") {
    return(list(preprocessed = X, normalized = X_norm,
                filtering = filtering, filtered = X_selected,
                filtered_normalized = X_norm_selected))
  }
  if (length(kept_genes) < config$cell_types) {
    stop("筛选后基因数少于类型数，无法继续投影。请检查数据和筛选参数。")
  }

  log_step("5/7 降维投影")
  projection <- algorithm$project_genes(X_norm_selected, config)
  log_step(sprintf("投影坐标：%d 个维度 × %d 个基因；中心化方差保留 %.2f%%",
                   nrow(projection$coordinates), ncol(projection$coordinates),
                   100 * projection$retained_variance))
  if (config$stop_after == "project") {
    return(list(preprocessed = X, normalized = X_norm,
                filtering = filtering, filtered = X_selected,
                filtered_normalized = X_norm_selected, projection = projection))
  }

  log_step("6/7 搜索角点")
  corners <- algorithm$find_corners(projection, X_norm_selected, config)
  if (config$stop_after == "corners") {
    return(list(preprocessed = X, normalized = X_norm,
                filtering = filtering, filtered = X_selected,
                filtered_normalized = X_norm_selected,
                projection = projection, corners = corners))
  }

  log_step("7/7 解卷积")
  fit <- algorithm$deconvolve_expression(
    X_selected, X_norm_selected, projection, corners, config
  )
  if (!is.list(fit) || !is.matrix(fit$proportions) ||
      !is.matrix(fit$signatures)) stop("解卷积结果必须包含两个矩阵。")
  P <- fit$proportions
  S <- fit$signatures
  if (!identical(dim(P), c(config$cell_types, ncol(X_selected))) ||
      !identical(dim(S), c(nrow(X_selected), config$cell_types))) {
    stop("解卷积矩阵维度不符合接口约定。")
  }
  if (any(!is.finite(P)) || any(!is.finite(S))) stop("输出包含非有限值。")
  if (any(P < 0) || any(S < 0)) stop("输出包含负值。")
  log_step(sprintf("解卷积：%d 个基因 × %d 个类型；归一化重建误差 %.6g",
                   nrow(S), ncol(S), fit$diagnostics$frobenius_error[1]))
  log_step(sprintf("各样本比例总和范围：%.6f–%.6f",
                   min(colSums(P)), max(colSums(P))))
  # 此处保留原始估计值，不自动强制比例之和等于 1。
  list(fit = fit, kept_genes = kept_genes, filtering = filtering,
       projection = projection, corners = corners)
}

# ---------- 6. 启动与保存 ----------
main <- function() {
  stopifnot(is.logical(config$run_analysis), length(config$run_analysis) == 1L,
            !is.na(config$run_analysis), config$cell_types >= 2L,
            config$cell_types == as.integer(config$cell_types),
            config$top_genes >= config$cell_types,
            config$permutations >= 1L,
            config$pvalue_cutoff > 0, config$pvalue_cutoff < 1)
  stopifnot(config$stop_after %in% c("preprocess", "normalize", "filter", "project", "corners", "full"),
            config$input_source %in% c("demo", "csv"))
  set.seed(config$seed)
  for (path in c("R", "data", config$results_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  log_step(paste("项目：", getwd()))

  if (!config$run_analysis) {
    log_step("框架检查模式：跳过数据读取和算法计算。")
    print(data.frame(
      function_name = names(modules), file = unname(modules),
      implemented_file = file.exists(modules), row.names = NULL
    ))
    cat("\n框架启动成功。FALSE 表示对应算法文件还未创建。\n")
    cat("全部 7 步已实现；stop_after = full 时运行完整流程。\n")
    cat("全部准备好后，将 config 中的 run_analysis 改为 TRUE。\n")
    return(invisible(NULL))
  }

  result <- run_pipeline(config)
  # 每次实际分析使用独立目录，避免覆盖旧结果。
  dir.create(config$results_dir, showWarnings = FALSE)
  run_dir <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-"),
                      tmpdir = config$results_dir)
  dir.create(run_dir)
  if (config$stop_after %in% c("preprocess", "normalize", "filter", "project", "corners")) {
    write.csv(result$preprocessed, file.path(run_dir, "preprocessed_expression.csv"))
    if (config$stop_after %in% c("normalize", "filter", "project", "corners")) {
      write.csv(result$normalized, file.path(run_dir, "normalized_expression.csv"))
    }
  } else {
    write.csv(result$fit$proportions, file.path(run_dir, "proportions.csv"))
    write.csv(result$fit$signatures, file.path(run_dir, "signatures.csv"))
    write.csv(result$fit$signatures_linear, file.path(run_dir, "signatures_linear.csv"))
    write.csv(result$fit$reconstructed, file.path(run_dir, "reconstructed_normalized.csv"))
    write.csv(result$fit$reconstructed_linear, file.path(run_dir, "reconstructed_linear.csv"))
    write.csv(result$fit$diagnostics, file.path(run_dir, "deconvolution_errors.csv"), row.names = FALSE)
    write.csv(data.frame(sample = names(result$fit$proportion_sums),
                         proportion_sum = unname(result$fit$proportion_sums)),
              file.path(run_dir, "proportion_sums.csv"), row.names = FALSE)
  }
  if (!is.null(result$filtering)) {
    write.csv(result$filtering$statistics, file.path(run_dir, "gene_statistics.csv"), row.names = FALSE)
    writeLines(result$filtering$kept_genes, file.path(run_dir, "kept_genes.txt"))
  }
  if (config$stop_after %in% c("filter", "project", "corners")) {
    write.csv(result$filtered, file.path(run_dir, "filtered_expression.csv"))
    write.csv(result$filtered_normalized, file.path(run_dir, "filtered_normalized_expression.csv"))
  }
  if (!is.null(result$projection)) {
    # CSV 以每行一个基因的方式导出；R 内存中坐标仍为“维度×基因”。
    write.csv(t(result$projection$coordinates), file.path(run_dir, "projection_coordinates.csv"))
    write.csv(result$projection$basis, file.path(run_dir, "projection_basis.csv"))
    write.csv(data.frame(
      component = seq_along(result$projection$singular_values),
      singular_value = result$projection$singular_values,
      variance_fraction = result$projection$variance_fraction
    ), file.path(run_dir, "projection_variance.csv"), row.names = FALSE)
  }
  if (!is.null(result$corners)) {
    write.csv(result$corners$endpoints, file.path(run_dir, "corner_endpoints.csv"))
    write.csv(t(result$corners$coordinates), file.path(run_dir, "corner_coordinates.csv"))
    write.csv(result$corners$distances, file.path(run_dir, "gene_corner_distances.csv"))
    write.csv(result$corners$candidates, file.path(run_dir, "corner_candidates.csv"), row.names = FALSE)
    write.csv(result$corners$search_trace, file.path(run_dir, "corner_search_trace.csv"), row.names = FALSE)
    write.csv(data.frame(component = names(result$corners$selected_tau),
                         tau = unname(result$corners$selected_tau)),
              file.path(run_dir, "corner_selected_tau.csv"), row.names = FALSE)
  }
  saveRDS(list(config = config, result = result), file.path(run_dir, "analysis.rds"))
  dput(config, file = file.path(run_dir, "config.R"))
  # 绘图单独加载，不改变算法结果；每次运行自动导出 PNG。
  plot_env <- new.env(parent = globalenv())
  sys.source("R/plot_results.R", envir = plot_env)
  plot_env$plot_results(result, config, file.path(run_dir, "plots"))
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  log_step(paste("分析完成，结果保存在：", run_dir))
  invisible(result)
}

analysis <- main()

