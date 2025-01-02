library(Rcpp)
library(parallel)
library(pbapply)

cppFunction('
NumericVector cpp_NoPAS_normalize_robust(NumericVector x) {
    int n = x.length();
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());
    
    double q10 = sorted[int(n * 0.1)];
    double q90 = sorted[int(n * 0.9)];
    double median = (n % 2 == 0) ? 
        (sorted[n/2 - 1] + sorted[n/2])/2.0 : 
        sorted[n/2];
    
    double range = q90 - q10;
    if(range == 0) range = 1;
    
    NumericVector normalized(n);
    for(int i = 0; i < n; i++) {
        normalized[i] = (x[i] - median) / range;
    }
    return normalized;
}
')

cppFunction('
NumericVector cpp_NoPAS_wilcox(NumericMatrix order_matrix,
                        IntegerVector pathway_indices,
                        int total_genes) {
    int n_samples = order_matrix.ncol();
    NumericVector wilcox_results(n_samples);
    
    int n_pathway_genes = pathway_indices.length();
    int n_other_genes = total_genes - n_pathway_genes;
    
    for(int j = 0; j < n_samples; j++) {
        NumericVector sample_order = order_matrix.column(j);
        double rank_sum = 0;
        for(int i = 0; i < n_pathway_genes; i++) {
            rank_sum += sample_order[pathway_indices[i]];
        }
        double U = rank_sum - (n_pathway_genes * (n_pathway_genes + 1))/2;
        wilcox_results[j] = U / (n_pathway_genes * n_other_genes);
    }
    
    return wilcox_results;
}
')

calculate_NoPAS_scores <- function(expr_matrix, 
                                   pathways, 
                                   scaler_method = "robust", 
                                   min_genes = 3,
                                   use_parallel = TRUE,
                                   n_cores = detectCores() - 1) {
  
  if (scaler_method == "robust") {
    NoPAS_normalized <- t(apply(expr_matrix, 1, cpp_NoPAS_normalize_robust))
    rownames(NoPAS_normalized) <- rownames(expr_matrix)
    colnames(NoPAS_normalized) <- colnames(expr_matrix)
  } else {
    stop("目前仅支持 Rcpp 版本的稳健归一化方法")
  }
  
  NoPAS_order_matrix <- apply(NoPAS_normalized, 2, rank, ties.method = "average")
  rownames(NoPAS_order_matrix) <- rownames(NoPAS_normalized)
  
  pathway_names <- names(pathways)
  
  if (!use_parallel) {
    message("[calculate_NoPAS_scores] 以单线程模式运行，并显示进度条...")
    scores_list <- pbapply::pblapply(
      X = pathway_names,
      FUN = function(pathway_name) {
        genes_in_pathway <- intersect(pathways[[pathway_name]], rownames(NoPAS_order_matrix))
        
        if (length(genes_in_pathway) < min_genes) {
          return(rep(NA, ncol(NoPAS_normalized)))
        }
        
        pathway_indices <- match(genes_in_pathway, rownames(NoPAS_order_matrix)) - 1
        
        cpp_NoPAS_wilcox(NoPAS_order_matrix, pathway_indices, nrow(NoPAS_order_matrix))
      }
    )
  } else {
    message("[calculate_NoPAS_scores] 以并行模式运行...")
    scores_list <- mclapply(
      X = pathway_names,
      FUN = function(pathway_name) {
        genes_in_pathway <- intersect(pathways[[pathway_name]], rownames(NoPAS_order_matrix))
        
        if (length(genes_in_pathway) < min_genes) {
          return(rep(NA, ncol(NoPAS_normalized)))
        }
        
        pathway_indices <- match(genes_in_pathway, rownames(NoPAS_order_matrix)) - 1
        
        cpp_NoPAS_wilcox(NoPAS_order_matrix, pathway_indices, nrow(NoPAS_order_matrix))
      },
      mc.cores = n_cores
    )
  }
  
  NoPAS_scores <- do.call(rbind, scores_list)
  rownames(NoPAS_scores) <- pathway_names
  colnames(NoPAS_scores) <- colnames(NoPAS_normalized)
  NoPAS_scores <- as.data.frame(NoPAS_scores)
  
  return(NoPAS_scores)
}

calculate_NoPAS_significance <- function(expr_matrix, 
                                         pathways, 
                                         n_permutations = 1000, 
                                         seed = 42,
                                         scaler_method = "robust",
                                         min_genes = 3,
                                         use_parallel = TRUE,
                                         n_cores = detectCores() - 1) {
  
  set.seed(seed)
  
  real_NoPAS_scores <- calculate_NoPAS_scores(
    expr_matrix = expr_matrix,
    pathways = pathways,
    scaler_method = scaler_method,
    min_genes = min_genes,
    use_parallel = use_parallel, 
    n_cores = n_cores
  )
  
  real_scores_vec <- as.vector(real_NoPAS_scores)
  
  n_pathways <- length(pathways)
  n_samples <- ncol(expr_matrix)
  
  message("[calculate_NoPAS_significance] 进行置换测试，置换次数 = ", n_permutations)
  
  if (!use_parallel) {
    message("[calculate_NoPAS_significance] 以单线程模式运行置换，并显示进度条...")
    random_scores_list <- pbapply::pblapply(
      X = seq_len(n_permutations),
      FUN = function(i) {
        shuffled_expr <- expr_matrix
        rownames(shuffled_expr) <- sample(rownames(shuffled_expr))
        
        shuffled_scores <- calculate_NoPAS_scores(
          expr_matrix = shuffled_expr,
          pathways = pathways,
          scaler_method = scaler_method,
          min_genes = min_genes,
          use_parallel = FALSE,  
          n_cores = 1
        )
        
        return(as.vector(shuffled_scores))
      }
    )
  } else {
    message("[calculate_NoPAS_significance] 以并行模式运行置换测试（无实时进度条）。")
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterExport(
      cl,
      varlist = c("expr_matrix", "pathways", "scaler_method", "min_genes",
                  "calculate_NoPAS_scores", "cpp_NoPAS_normalize_robust", "cpp_NoPAS_wilcox"),
      envir = environment()
    )
    
    random_scores_list <- parLapply(
      cl,
      X = seq_len(n_permutations),
      fun = function(i) {
        shuffled_expr <- expr_matrix
        rownames(shuffled_expr) <- sample(rownames(shuffled_expr))
        
        shuffled_scores <- calculate_NoPAS_scores(
          expr_matrix = shuffled_expr,
          pathways = pathways,
          scaler_method = scaler_method,
          min_genes = min_genes,
          use_parallel = FALSE,
          n_cores = 1
        )
        
        return(as.vector(shuffled_scores))
      }
    )
  }
  
  random_NoPAS_scores <- do.call(rbind, random_scores_list)
  
  pvalues <- matrix(
    colMeans(
      random_NoPAS_scores >= matrix(
        real_scores_vec,
        nrow = n_permutations,
        ncol = length(real_scores_vec),
        byrow = TRUE
      )
    ),
    nrow = n_pathways,
    ncol = n_samples,
    dimnames = dimnames(real_NoPAS_scores)
  )
  
  pvalues <- as.data.frame(pvalues)
  
  return(list(
    scores = real_NoPAS_scores,
    pvalues = pvalues
  ))
}
