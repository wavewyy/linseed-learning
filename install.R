options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 120)
.libPaths(c("D:/single cell项目/Linseed/library", .libPaths()))
Sys.setenv(R_LIBS_USER="D:/single cell项目/Linseed/library", RTOOLS45_HOME="D:/single cell项目/Linseed/rtools45")
Sys.setenv(PATH=paste("D:/single cell项目/Linseed/rtools45/usr/bin", "D:/single cell项目/Linseed/rtools45/x86_64-w64-mingw32.static.posix/bin", Sys.getenv("PATH"), sep=";"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", type="binary")
options(repos=c(BioCsoft="https://bioconductor.org/packages/3.22/bioc", CRAN="https://cloud.r-project.org"))
deps <- c("Rcpp", "RcppArmadillo", "BH", "GEOquery", "Biobase", "corpcor", "combinat", "NMF", "R6", "dplyr", "ggplot2", "reshape2", "preprocessCore", "progress", "Rtsne", "languageserver")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, type="binary")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
install.packages("D:/single cell项目/Linseed/source", repos=NULL, type="source", INSTALL_opts="--no-multiarch")
stopifnot(requireNamespace("linseed", quietly=TRUE))
writeLines(capture.output(sessionInfo()), "D:/single cell项目/Linseed/logs/sessionInfo.txt")
cat("LINSEED_INSTALL_OK\n")


