list(run_analysis = TRUE, stop_after = "corners", input_source = "demo", 
    samples = 10:42, input_scale = "auto", input_file = "data/bulk_expression.csv", 
    cell_types = 3L, top_genes = 500L, permutations = 100L, pvalue_cutoff = 0.01, 
    spearman_threshold = 0, corner_taus = c(1, 0.25, 0.0625, 
    0.015625, 0.00390625), sisal_iterations = 30L, seed = 20261022L, 
    results_dir = "results")
