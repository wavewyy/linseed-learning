list(run_analysis = TRUE, stop_after = "filter", input_source = "demo", 
    samples = 10:42, input_scale = "auto", input_file = "data/bulk_expression.csv", 
    cell_types = 3L, top_genes = 500L, permutations = 100L, pvalue_cutoff = 0.01, 
    spearman_threshold = 0, seed = 20261022L, results_dir = "results")
