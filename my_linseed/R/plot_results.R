# 学习版结果图：只读取计算结果，不改变比例或重新运行算法。
# 输出 PNG，可直接在 VS Code 打开。使用 R 自带绘图，无需额外安装包。
plot_results <- function(result, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  palette <- grDevices::hcl.colors(config$cell_types, "Dark 3")
  draw <- function(filename, code) {
    grDevices::png(file.path(output_dir, filename), width = 1600, height = 1000,
                   res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (.Platform$OS.type == "windows") {
      grDevices::windowsFonts(study = grDevices::windowsFont("Microsoft YaHei"))
      par(family = "study")
    }
    par(mar = c(5, 5, 4, 2), las = 1, col.axis = "#394553",
        col.lab = "#243447", col.main = "#243447")
    force(code)
  }
  if (!is.null(result$fit)) {
    P <- result$fit$proportions
    draw("01_proportions.png", {
      par(mar = c(11, 5, 5, 2))
      barplot(P, col = palette, border = NA, las = 2, cex.names = .65,
              ylim = c(0, max(1, colSums(P)) * 1.24),
              main = "各样本的类型比例估计", ylab = "估计比例（未强制归一化）")
      abline(h = 1, lty = 2, col = "#697586")
      legend("top", legend = rownames(P), fill = palette, bty = "n", horiz = TRUE,
             cex = .85)
      mtext("虚线表示总和为 1；Component 编号尚未对应具体细胞类型", side = 3,
            line = .3, cex = .75)
    })
  }
  if (!is.null(result$projection)) {
    pr <- result$projection
    draw("02_projection_corners.png", {
      # k=2 时仅有一个中心化方向，另一个轴为均值正交方向。
      dims <- 1:2
      xy <- t(pr$coordinates[dims, , drop = FALSE])
      E <- if (!is.null(result$corners)) t(result$corners$coordinates[dims, , drop = FALSE]) else NULL
      limits <- rbind(xy, E)
      plot(xy, xlim = extendrange(limits[, 1], f = .18),
           ylim = extendrange(limits[, 2], f = .18), asp = 1,
           pch = 16, cex = .65, col = "#4285A480",
           xlab = rownames(pr$coordinates)[1], ylab = rownames(pr$coordinates)[2],
           main = "基因投影与搜索到的角点")
      if (!is.null(E)) {
        if (nrow(E) >= 3) {
          hull <- chull(E)
          polygon(E[hull, , drop = FALSE], border = "#697586", lty = 2)
        }
        points(E, pch = 17, col = palette, cex = 1.6)
        text(E, labels = rownames(E), pos = 3, cex = .8, col = palette)
      }
      mtext("圆点为基因，三角形为角点；类型数大于 3 时这里只显示二维切片", side = 3,
            line = .3, cex = .75)
    })
    draw("03_variance.png", {
      n <- min(15L, length(pr$variance_fraction))
      barplot(100 * pr$variance_fraction[seq_len(n)], names.arg = seq_len(n),
              col = ifelse(seq_len(n) <= config$cell_types - 1L, "#258C91", "#C9D4DF"),
              border = NA, xlab = "中心化矩阵的奇异方向", ylab = "方差占比（%）",
              main = "降维保留了多少信息？")
      mtext(sprintf("绿色为保留方向，合计保留 %.2f%% 的中心化方差",
                    100 * pr$retained_variance), side = 3, line = .3, cex = .8)
    })
  }
  if (!is.null(result$filtering)) {
    st <- result$filtering$statistics
    draw("04_gene_selection.png", {
      plot(st$degree, st$strength, pch = 16, cex = .8,
           col = ifelse(st$kept, "#D56B35", "#AAB7C480"),
           xlab = "连接数（degree）", ylab = "连接权重总和（strength）",
           main = "显著性检验与基因筛选")
      legend("topleft", c(sprintf("保留 %d 个基因", sum(st$kept)),
                           sprintf("未保留 %d 个基因", sum(!st$kept))),
             col = c("#D56B35", "#AAB7C4"), pch = 16, bty = "n")
      mtext(sprintf("筛选条件：p < %g；颜色来自置换检验结果", config$pvalue_cutoff),
            side = 3, line = .3, cex = .8)
    })
  }
  invisible(output_dir)
}
