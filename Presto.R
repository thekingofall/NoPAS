library(Rcpp)

cppFunction('
NumericVector cpp_normalize_robust(NumericVector x) {
    int n = x.length();
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());
    
    double q10 = sorted[int(n * 0.1)];
    double q90 = sorted[int(n * 0.9)];
    double med = (n % 2 == 0) ? 
        (sorted[n/2 - 1] + sorted[n/2])/2.0 : 
        sorted[n/2];
    
    double range = q90 - q10;
    if(range == 0) range = 1;
    
    NumericVector result(n);
    for(int i = 0; i < n; i++) {
        result[i] = (x[i] - med) / range;
    }
    return result;
}
')

cppFunction('
NumericVector cpp_wilcox(NumericMatrix ranks_matrix,
                        IntegerVector pathway_idx,
                        int n_genes) {
    int n_samples = ranks_matrix.ncol();
    NumericVector results(n_samples);
    
    int n1 = pathway_idx.length();
    int n2 = n_genes - n1;
    
    for(int j = 0; j < n_samples; j++) {
        NumericVector col = ranks_matrix.column(j);
        double rank_sum = 0;
        for(int i = 0; i < n1; i++) {
            rank_sum += col[pathway_idx[i]];
        }
        double U = rank_sum - (n1 * (n1 + 1))/2;
        results[j] = U/(n1 * n2);
    }
    
    return results;
}
')


calculate_presto_scores <- function(expr_matrix, pathways, scaler_method = "robust", min_genes = 3) {

  if(scaler_method == "robust") {
    norm_expr <- t(apply(expr_matrix, 1, cpp_normalize_robust))
    rownames(norm_expr) <- rownames(expr_matrix)
    colnames(norm_expr) <- colnames(expr_matrix)
  } else {
    stop("Currently only robust scaling is implemented in Rcpp version")
  }
  

  ranks_matrix <- apply(norm_expr, 2, rank, ties.method = "average")
  rownames(ranks_matrix) <- rownames(norm_expr)
  

  scores <- matrix(NA, 
                  nrow = length(pathways), 
                  ncol = ncol(norm_expr),
                  dimnames = list(names(pathways), colnames(norm_expr)))
  

  for(pathway_name in names(pathways)) {
    genes_in_pathway <- intersect(pathways[[pathway_name]], rownames(ranks_matrix))
    if(length(genes_in_pathway) < min_genes) next
    
    pathway_idx <- match(genes_in_pathway, rownames(ranks_matrix)) - 1
  

    scores[pathway_name, ] <- cpp_wilcox(ranks_matrix, pathway_idx, nrow(ranks_matrix))
  }
  
  return(scores)
}


calculate_significance <- function(expr_matrix, 
                                 pathways, 
                                 n_permutations = 1000, 
                                 seed = 42,
                                 scaler_method = "robust",
                                 min_genes = 3) {
  
  set.seed(seed)
  real_scores <- calculate_presto_scores(expr_matrix, pathways, scaler_method, min_genes)
  

  n_pathways <- length(pathways)
  n_samples <- ncol(expr_matrix)

  library(parallel)
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)

  clusterExport(cl, c("calculate_presto_scores", "cpp_normalize_robust", "cpp_wilcox"))
  
  random_scores <- parLapply(cl, 1:n_permutations, function(i) {
    shuffled_expr <- expr_matrix
    rownames(shuffled_expr) <- sample(rownames(shuffled_expr))
    as.vector(calculate_presto_scores(shuffled_expr, pathways, scaler_method, min_genes))
  })
  
  stopCluster(cl)

  random_scores <- do.call(rbind, random_scores)
  real_scores_vec <- as.vector(real_scores)
  
  pvalues <- matrix(
    colMeans(random_scores >= matrix(
      real_scores_vec,
      nrow = n_permutations,
      ncol = length(real_scores_vec),
      byrow = TRUE
    )),
    nrow = n_pathways,
    ncol = n_samples,
    dimnames = dimnames(real_scores)
  )
  
  return(list(
    scores = real_scores,
    pvalues = pvalues
  ))
}
