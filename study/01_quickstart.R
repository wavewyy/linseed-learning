# Linseed 原教程参数：对照 source/README.md。
# 与原文一样从 GEO 读取 GSE19830，首次运行可能需要联网下载。
# 在 VS Code 中用 Ctrl+Shift+P -> R: Create R terminal，选中代码按 Ctrl+Enter。
.libPaths(c("D:/single cell项目/Linseed/library", .libPaths()))
setwd("D:/single cell项目/Linseed")
library(linseed)
# 行是基因，列是样本；10:42 是混合样本。
# 原文没有指定随机种子，因此不承诺每次结果完全一致。
lo <- LinseedObject$new("GSE19830", samples=10:42, topGenes=10000)
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100)
lo$filterDatasetByPval(0.01)
lo$setCellTypeNumber(3)
lo$project("filtered")
# 按原文调用默认搜索：taus = 2^seq(0, -20, -1)，sisalIter = 100。
lo$smartSearchCorners(dataset="filtered", error="norm")
lo$deconvolveByEndpoints()
stopifnot(all(is.finite(lo$proportions)))
# 输出到单独目录，与以前 500 基因的 quickstart 结果区分。
dir.create("results/quickstart_original", recursive=TRUE, showWarnings=FALSE)
write.csv(lo$proportions, "results/quickstart_original/proportions.csv")
saveRDS(lo, "results/quickstart_original/linseed_object.rds")
png("results/quickstart_original/proportions.png", width=1200, height=700)
print(plotProportions(lo$proportions))
dev.off()
writeLines(capture.output(sessionInfo()), "results/quickstart_original/sessionInfo.txt")
cat("QUICKSTART_OK: results/quickstart_original\n")

