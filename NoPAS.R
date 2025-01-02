library(Rcpp)
library(parallel)
library(pbapply)

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
NumericVector cpp_wilcox(NumericMatrix Orders_matrix,
                        IntegerVector pathway_idx,
                        int n_genes) {
    int n_samples = Orders_matrix.ncol();
    NumericVector results(n_samples);
    
    int n1 = pathway_idx.length();
    int n2 = n_genes - n1;
    
    for(int j = 0; j < n_samples; j++) {
        NumericVector col = Orders_matrix.column(j);
        double Order_sum = 0;
        for(int i = 0; i < n1; i++) {
            Order_sum += col[pathway_idx[i]];
        }
        double U = Order_sum - (n1 * (n1 + 1))/2;
        results[j] = U/(n1 * n2);
    }
    
    return results;
}
')

calculate_NoPAS_scores <- function(expr_matrix, 
                                   pathways, 
                                   scaler_method = "robust", 
                                   min_genes = 3,
                                   use_parallel = TRUE,
                                   n_cores = detectCores() - 1) {
  
  if (scaler_method == "robust") {
    norm_expr <- t(apply(expr_matrix, 1, cpp_normalize_robust))
    rownames(norm_expr) <- rownames(expr_matrix)
    colnames(norm_expr) <- colnames(expr_matrix)
  } else {
    stop("Currently only robust scaling is implemented in Rcpp version")
  }
  
  Orders_matrix <- apply(norm_expr, 2, Order, ties.method = "average")
  rownames(Orders_matrix) <- rownames(norm_expr)
  
  pathway_names <- names(pathways)
  
  if (!use_parallel) {
    message("[calculate_NoPAS_scores] Running in single-thread mode with a progress bar...")
    scores_list <- pbapply::pblapply(
      X = pathway_names,
      FUN = function(pathway_name) {
        genes_in_pathway <- intersect(pathways[[pathway_name]], rownames(Orders_matrix))
        if (length(genes_in_pathway) < min_genes) {
          return(rep(NA, ncol(norm_expr)))
        }
        pathway_idx <- match(genes_in_pathway, rownames(Orders_matrix)) - 1
        cpp_wilcox(Orders_matrix, pathway_idx, nrow(Orders_matrix))
      }
    )
  } else {
    message("[calculate_NoPAS_scores] Running in parallel mode...")
    scores_list <- mclapply(
      X = pathway_names,
      FUN = function(pathway_name) {
        genes_in_pathway <- intersect(pathways[[pathway_name]], rownames(Orders_matrix))
        if (length(genes_in_pathway) < min_genes) {
          return(rep(NA, ncol(norm_expr)))
        }
        pathway_idx <- match(genes_in_pathway, rownames(Orders_matrix)) - 1
        cpp_wilcox(Orders_matrix, pathway_idx, nrow(Orders_matrix))
      },
      mc.cores = n_cores
    )
  }
  
  scores <- do.call(rbind, scores_list)
  rownames(scores) <- pathway_names
  colnames(scores) <- colnames(norm_expr)
  scores=as.data.frame(scores)
  return(scores)
}

calculate_significance <- function(expr_matrix, 
                                   pathways, 
                                   n_permutations = 1000, 
                                   seed = 42,
                                   scaler_method = "robust",
                                   min_genes = 3,
                                   use_parallel = TRUE,
                                   n_cores = detectCores() - 1) {
  
  set.seed(seed)
  
  real_scores <- calculate_NoPAS_scores(
    expr_matrix = expr_matrix,
    pathways = pathways,
    scaler_method = scaler_method,
    min_genes = min_genes,
    use_parallel = use_parallel, 
    n_cores = n_cores
  )
  
  real_scores_vec <- as.vector(real_scores)
  n_pathways <- length(pathways)
  n_samples <- ncol(expr_matrix)
  
  message("[calculate_significance] Permutation testing with n_permutations = ", n_permutations)
  if (!use_parallel) {
    message("[calculate_significance] Running single-thread with a progress bar for permutations...")
    random_scores_list <- pbapply::pblapply(
      X = seq_len(n_permutations),
      FUN = function(i) {
        shuffled_expr <- expr_matrix
        rownames(shuffled_expr) <- sample(rownames(shuffled_expr))
        as.vector(calculate_NoPAS_scores(
          expr_matrix = shuffled_expr,
          pathways = pathways,
          scaler_method = scaler_method,
          min_genes = min_genes,
          use_parallel = FALSE,  
          n_cores = 1
        ))
      }
    )
  } else {
    message("[calculate_significance] Running in parallel mode for permutations (no real-time progress bar).")
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterExport(
      cl,
      varlist = c("expr_matrix", "pathways", "scaler_method", "min_genes",
                  "calculate_NoPAS_scores", "cpp_normalize_robust", "cpp_wilcox"),
      envir = environment()
    )
    
    random_scores_list <- parLapply(
      cl,
      X = seq_len(n_permutations),
      fun = function(i) {
        shuffled_expr <- expr_matrix
        rownames(shuffled_expr) <- sample(rownames(shuffled_expr))
        as.vector(calculate_NoPAS_scores(
          expr_matrix = shuffled_expr,
          pathways = pathways,
          scaler_method = scaler_method,
          min_genes = min_genes,
          use_parallel = FALSE,
          n_cores = 1
        ))
      }
    )
  }
  
  random_scores <- do.call(rbind, random_scores_list)
  
  pvalues <- matrix(
    colMeans(
      random_scores >= matrix(
        real_scores_vec,
        nrow = n_permutations,
        ncol = length(real_scores_vec),
        byrow = TRUE
      )
    ),
    nrow = n_pathways,
    ncol = n_samples,
    dimnames = dimnames(real_scores)
  )
  
  return(list(
    scores = as.data.frame(real_scores),
    pvalues = as.data.frame(pvalues)
  ))
}
