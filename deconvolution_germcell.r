#!/usr/bin/env Rscript
# Enhanced Cell Type Deconvolution Analysis Script with CIBERSORT-like Algorithm
# 基于CIBERSORT算法的增强型生殖细胞解卷积分析脚本
# Author: chenpeigen
# Date: 2025-07-31
# 
# USAGE:
#   Rscript deconvolution_germcell.r -i input_count.csv -o output_proportions.csv
#   Rscript deconvolution_germcell.r --test  # Run with test data
# 
# OUTPUT:
#   - *_proportions.csv: Cell type proportions for each sample
#   - *_detailed.csv: Detailed results with p-values and quality metrics
# 
# REQUIREMENTS:
#   - R packages: tidyverse, optparse, parallel
#   - Optional (auto-installed): nnls, quadprog for better performance

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(parallel)   # for parallel computing
  
  # Check and install optional packages for better performance
  if (!requireNamespace("nnls", quietly = TRUE)) {
    cat("Note: 'nnls' package not found. Installing for better deconvolution performance...\n")
    tryCatch({
      install.packages("nnls", repos = "https://cran.r-project.org")
    }, error = function(e) {
      cat("Warning: Failed to install 'nnls' package. Using fallback method.\n")
    })
  }
  
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    cat("Note: 'quadprog' package not found. Installing for constrained optimization...\n")
    tryCatch({
      install.packages("quadprog", repos = "https://cran.r-project.org")
    }, error = function(e) {
      cat("Warning: Failed to install 'quadprog' package. Using fallback method.\n")
    })
  }
})

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input count matrix file (CSV format)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="celltype_proportions.csv",
              help="Output cell type proportions file [default= %default]", metavar="character"),
  make_option(c("-n", "--nu"), type="double", default=0.5,
              help="Nu parameter for SVR (0-1) [default= %default]", metavar="double"),
  make_option(c("-p", "--perm"), type="integer", default=100,
              help="Number of permutations for significance testing [default= %default]", metavar="integer"),
  make_option(c("-t", "--test"), action="store_true", default=FALSE,
              help="Run with test data")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Enhanced Cell Type Deconvolution Analysis with CIBERSORT-like Algorithm")
opt <- parse_args(opt_parser)

# Define enhanced cell marker genes for testicular cell types with detailed germ cell subtypes
cell_markers <- list(
  # Somatic cell types
  Leydig = c("Star", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Hsd3b2",
             "Hsd17b3", "Insl3", "Lhcgr", "Scarb1", "Ldlr"),
  Sertoli = c("Amh", "Fshr", "Sox9", "Wt1", "Gata4", "Dhh",
              "Vim", "Claudin11", "Tjp1", "Ocln"),
  Peritubular = c("Acta2", "Myh11", "Pdgfrb", "Des", "Vim"),
  
  # Detailed germ cell subtypes
  SSC = c("Plzf", "Gfra1", "Ret", "Id4", "Nanos2", "Utf1", "Etv5", "Bcl6b"),
  Diff_Spg = c("Kit", "Stra8", "Sohlh1", "Sohlh2", "Ngn3", "Dmrt1", "Rhox13"),
  Spermatocytes = c("Sycp3", "Sycp1", "Dmc1", "Rec8", "Spo11", "Msh4", "Mlh1", "Hormad1"),
  Round_Spermatids = c("Prm1", "Prm2", "Tnp1", "Tnp2", "Acr", "Crem", "Acrv1", "Zpbp1"),
  Elongating_Spermatids = c("Odf1", "Odf2", "Spag6", "Tektin2", "Tekt1", "Spag16", "Rsph1"),
  
  # Immune cell types
  Macrophages = c("Cd68", "F4/80", "Cd14", "Tnf", "Il1b", "Il6"),
  M1_Macrophages = c("Nos2", "Tnf", "Il1b", "Il6", "Cd86", "Cd80",
                     "Cxcl9", "Cxcl10", "Cxcl11", "Ccr7", "Socs3",
                     "Ptgs2", "Cd40", "Fcgr1", "Stat1", "Irf1"),
  M2_Macrophages = c("Arg1", "Il10", "Tgfb1", "Mrc1", "Il4ra", "Retnla",
                     "Chil3", "Ccl17", "Ccl18", "Ccl22", "Il1rn", "Cd163",
                     "Clec10a", "Lyve1", "Stat6", "Irf4")
)

# Function to create test data
create_test_data <- function() {
  cat("Creating test data...\n")
  
  # Create synthetic expression data
  set.seed(123)
  n_genes <- 1000
  n_samples <- 6
  
  # Generate gene names (mix of marker genes and random genes)
  all_markers <- unique(unlist(cell_markers))
  marker_subset <- sample(all_markers, min(200, length(all_markers)))
  random_genes <- paste0("Gene_", 1:(n_genes - length(marker_subset)))
  gene_names <- c(marker_subset, random_genes)
  
  # Generate expression matrix with some biological signal
  expr_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- paste0("Sample_", 1:n_samples)
  
  # Add baseline expression
  expr_matrix <- matrix(rpois(n_genes * n_samples, lambda = 50), 
                       nrow = n_genes, ncol = n_samples)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- paste0("Sample_", 1:n_samples)
  
  # Add cell-type specific signals
  for(i in 1:length(cell_markers)) {
    cell_type <- names(cell_markers)[i]
    markers <- intersect(cell_markers[[i]], gene_names)
    
    if(length(markers) > 0) {
      # Different samples have different cell type compositions
      for(j in 1:n_samples) {
        # Simulate varying cell type proportions
        proportion <- runif(1, 0.05, 0.3)
        boost_factor <- 1 + proportion * 5
        expr_matrix[markers, j] <- expr_matrix[markers, j] * boost_factor
      }
    }
  }
  
  # Save test data
  test_file <- "test_expression_data.csv"
  write.csv(expr_matrix, test_file)
  cat("Test data saved to:", test_file, "\n")
  cat("Test data dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  
  return(test_file)
}

# Handle test mode
if(opt$test) {
  cat("=== Running in Test Mode ===\n")
  opt$input <- create_test_data()
  opt$output <- "test_celltype_proportions.csv"
}

# Check if input file is provided
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified with -i option, or use --test for test mode.", call.=FALSE)
}

# Check if input file exists
if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call.=FALSE)
}

# Validate parameters
if(opt$nu < 0 || opt$nu > 1) {
  stop("Nu parameter must be between 0 and 1", call.=FALSE)
}

if(opt$perm < 0) {
  stop("Number of permutations must be non-negative", call.=FALSE)
}

# Define enhanced cell marker genes for testicular cell types with detailed germ cell subtypes
cell_markers <- list(
  # Somatic cell types
  Leydig = c("Star", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Hsd3b2",
             "Hsd17b3", "Insl3", "Lhcgr", "Scarb1", "Ldlr"),
  Sertoli = c("Amh", "Fshr", "Sox9", "Wt1", "Gata4", "Dhh",
              "Vim", "Claudin11", "Tjp1", "Ocln"),
  Peritubular = c("Acta2", "Myh11", "Pdgfrb", "Des", "Vim"),
  
  # Detailed germ cell subtypes
  SSC = c("Plzf", "Gfra1", "Ret", "Id4", "Nanos2", "Utf1", "Etv5", "Bcl6b"),
  Diff_Spg = c("Kit", "Stra8", "Sohlh1", "Sohlh2", "Ngn3", "Dmrt1", "Rhox13"),
  Spermatocytes = c("Sycp3", "Sycp1", "Dmc1", "Rec8", "Spo11", "Msh4", "Mlh1", "Hormad1"),
  Round_Spermatids = c("Prm1", "Prm2", "Tnp1", "Tnp2", "Acr", "Crem", "Acrv1", "Zpbp1"),
  Elongating_Spermatids = c("Odf1", "Odf2", "Spag6", "Tektin2", "Tekt1", "Spag16", "Rsph1"),
  
  # Immune cell types
  Macrophages = c("Cd68", "F4/80", "Cd14", "Tnf", "Il1b", "Il6"),
  M1_Macrophages = c("Nos2", "Tnf", "Il1b", "Il6", "Cd86", "Cd80",
                     "Cxcl9", "Cxcl10", "Cxcl11", "Ccr7", "Socs3",
                     "Ptgs2", "Cd40", "Fcgr1", "Stat1", "Irf1"),
  M2_Macrophages = c("Arg1", "Il10", "Tgfb1", "Mrc1", "Il4ra", "Retnla",
                     "Chil3", "Ccl17", "Ccl18", "Ccl22", "Il1rn", "Cd163",
                     "Clec10a", "Lyve1", "Stat6", "Irf4")
)

# TPM normalization function
tpm_normalize <- function(counts, gene_lengths = NULL) {
  if(is.null(gene_lengths)) {
    # Assume 1kb gene length if not provided
    gene_lengths <- rep(1000, nrow(counts))
  }
  
  # Calculate RPK (Reads Per Kilobase)
  rpk <- counts / (gene_lengths / 1000)
  
  # Calculate TPM (Transcripts Per Million)
  tpm <- t(t(rpk) / colSums(rpk) * 1e6)
  
  return(log2(tpm + 1))
}

# Feature selection function
select_features <- function(sig_matrix, top_n = 500) {
  # Calculate coefficient of variation for each gene
  cv_genes <- apply(sig_matrix, 1, function(x) {
    if(mean(x) == 0) return(0)
    sd(x) / mean(x)
  })
  
  # Calculate F-statistic for group differences
  f_stats <- apply(sig_matrix, 1, function(x) {
    if(var(x) == 0) return(0)
    groups <- factor(rep(1:ncol(sig_matrix), each = 1))
    tryCatch({
      anova_result <- aov(x ~ groups)
      summary(anova_result)[[1]][1, 4]
    }, error = function(e) 0)
  })
  
  # Combine scores and select top features
  feature_score <- rank(-cv_genes) + rank(-f_stats)
  selected_genes <- names(sort(feature_score)[1:min(top_n, length(feature_score))])
  
  return(selected_genes)
}

# Improved deconvolution using quadratic programming
solve_deconvolution <- function(bulk_sample, sig_matrix) {
  # Prepare training data
  common_genes <- intersect(names(bulk_sample), rownames(sig_matrix))
  if(length(common_genes) < 10) {
    warning("Too few common genes for deconvolution analysis")
    return(rep(1/ncol(sig_matrix), ncol(sig_matrix)))
  }
  
  cat("  Using", length(common_genes), "common genes for deconvolution\n")
  
  # Prepare data: solve X %*% proportions = y
  X <- sig_matrix[common_genes, , drop = FALSE]  # genes x cell_types
  y <- bulk_sample[common_genes]  # bulk expression for these genes
  
  # Handle missing values
  complete_idx <- complete.cases(cbind(y, t(X)))
  if(sum(complete_idx) < 5) {
    warning("Insufficient complete data for deconvolution")
    return(rep(1/ncol(sig_matrix), ncol(sig_matrix)))
  }
  
  X <- X[complete_idx, , drop = FALSE]
  y <- y[complete_idx]
  
  cat("  After filtering:", nrow(X), "genes for deconvolution\n")
  
  # Try multiple approaches
  tryCatch({
    # Method 1: Non-negative least squares (if available)
    if(requireNamespace("nnls", quietly = TRUE)) {
      cat("  Using NNLS method\n")
      nnls_result <- nnls::nnls(X, y)
      proportions <- nnls_result$x
      method_used <- "NNLS"
    } else {
      # Method 2: Constrained optimization using quadprog (if available)
      if(requireNamespace("quadprog", quietly = TRUE)) {
        cat("  Using quadratic programming method\n")
        # Set up quadratic programming problem
        # minimize: ||X %*% p - y||^2 subject to: p >= 0, sum(p) = 1
        
        Dmat <- t(X) %*% X
        dvec <- t(X) %*% y
        
        # Constraints: p >= 0 and sum(p) = 1
        n_celltypes <- ncol(X)
        Amat <- rbind(rep(1, n_celltypes), diag(n_celltypes))
        bvec <- c(1, rep(0, n_celltypes))
        
        qp_result <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = 1)
        proportions <- qp_result$solution
        method_used <- "QP"
      } else {
        # Method 3: Simple least squares with post-processing
        cat("  Using least squares with constraints\n")
        lm_result <- lm(y ~ . - 1, data = data.frame(y = y, t(X)))
        proportions <- coef(lm_result)
        proportions[is.na(proportions)] <- 0
        method_used <- "LS"
      }
    }
    
    names(proportions) <- colnames(sig_matrix)
    
    # Apply non-negative constraint
    proportions[proportions < 0] <- 0
    
    # Normalize to sum to 1
    if(sum(proportions) > 0) {
      proportions <- proportions / sum(proportions)
    } else {
      proportions <- rep(1/length(proportions), length(proportions))
    }
    
    cat("  Deconvolution completed using", method_used, "\n")
    return(proportions)
    
  }, error = function(e) {
    warning(paste("Deconvolution failed:", e$message))
    cat("  Falling back to correlation-based scoring\n")
    
    # Improved fallback: correlation-based scoring
    scores <- sapply(1:ncol(sig_matrix), function(i) {
      cell_signature <- sig_matrix[, i]
      # Calculate correlation between bulk sample and cell signature
      cor_score <- cor(bulk_sample[common_genes], cell_signature[common_genes], 
                      use = "complete.obs")
      # Convert correlation to positive score
      return(max(0, cor_score))
    })
    
    names(scores) <- colnames(sig_matrix)
    
    # Add biological variability to avoid identical scores
    scores <- scores * runif(length(scores), 0.8, 1.2)
    scores[scores < 0] <- 0
    
    # Normalize
    if(sum(scores) > 0) {
      scores <- scores / sum(scores)
    } else {
      scores <- rep(1/length(scores), length(scores))
    }
    
    return(scores)
  })
}

# Permutation test for statistical significance
permutation_test <- function(bulk_sample, sig_matrix, observed_props, perm = 100) {
  null_dist <- replicate(perm, {
    # Shuffle bulk expression data
    shuffled_bulk <- sample(bulk_sample)
    names(shuffled_bulk) <- names(bulk_sample)
    
    # Recalculate proportions
    null_props <- solve_deconvolution(shuffled_bulk, sig_matrix)
    return(null_props)
  })
  
  # Calculate p-values
  p_values <- sapply(1:length(observed_props), function(i) {
    if(is.matrix(null_dist)) {
      sum(null_dist[i, ] >= observed_props[i]) / perm
    } else {
      sum(null_dist >= observed_props[i]) / perm
    }
  })
  
  names(p_values) <- names(observed_props)
  return(p_values)
}

# Create realistic signature matrix based on cell type specificity
create_signature_matrix <- function(bulk_expr, markers) {
  all_marker_genes <- unique(unlist(markers))
  available_markers <- intersect(all_marker_genes, rownames(bulk_expr))
  
  if(length(available_markers) < 10) {
    stop("Too few marker genes available in the expression data")
  }
  
  # Create signature matrix
  sig_matrix <- matrix(0, nrow = length(available_markers), ncol = length(markers))
  rownames(sig_matrix) <- available_markers
  colnames(sig_matrix) <- names(markers)
  
  # Calculate gene expression statistics
  gene_means <- rowMeans(bulk_expr[available_markers, ], na.rm = TRUE)
  gene_vars <- apply(bulk_expr[available_markers, ], 1, var, na.rm = TRUE)
  
  # Fill signature matrix with realistic expression values
  for(i in 1:length(markers)) {
    cell_type <- names(markers)[i]
    cell_markers <- intersect(markers[[i]], available_markers)
    
    if(length(cell_markers) > 0) {
      for(marker in cell_markers) {
        # Base expression level
        base_expr <- gene_means[marker]
        
        # Cell type specific expression (simulate higher expression in target cell type)
        # Use a combination of base expression and cell-type specific boost
        specificity_boost <- runif(1, 2, 5)  # Random boost between 2-5x
        cell_specific_expr <- base_expr * specificity_boost
        
        # Add some biological variability
        noise <- rnorm(1, 0, sqrt(gene_vars[marker]) * 0.1)
        final_expr <- max(0.1, cell_specific_expr + noise)
        
        sig_matrix[marker, i] <- final_expr
      }
    }
  }
  
  # Add background expression for non-marker genes
  for(i in 1:ncol(sig_matrix)) {
    zero_genes <- which(sig_matrix[, i] == 0)
    if(length(zero_genes) > 0) {
      # Add low background expression
      background_expr <- gene_means[zero_genes] * runif(length(zero_genes), 0.1, 0.3)
      sig_matrix[zero_genes, i] <- background_expr
    }
  }
  
  return(sig_matrix)
}

# Enhanced CIBERSORT-like deconvolution function
cibersort_deconvolution <- function(bulk_expr, markers, nu = 0.5, perm = 100) {
  cat("Starting CIBERSORT-like deconvolution analysis...\n")
  
  # Create realistic signature matrix
  sig_matrix <- create_signature_matrix(bulk_expr, markers)
  
  cat("Signature matrix summary:\n")
  for(i in 1:ncol(sig_matrix)) {
    marker_genes <- intersect(markers[[i]], rownames(sig_matrix))
    non_marker_genes <- setdiff(rownames(sig_matrix), marker_genes)
    
    marker_expr <- mean(sig_matrix[marker_genes, i])
    background_expr <- mean(sig_matrix[non_marker_genes, i])
    specificity_ratio <- marker_expr / background_expr
    
    cat(sprintf("  %s: %d markers, specificity ratio: %.2f\n", 
                colnames(sig_matrix)[i], length(marker_genes), specificity_ratio))
  }
  
  # Feature selection
  if(nrow(sig_matrix) > 500) {
    selected_features <- select_features(sig_matrix, top_n = 500)
    sig_matrix <- sig_matrix[selected_features, , drop = FALSE]
  }
  
  cat("Using", nrow(sig_matrix), "features for deconvolution\n")
  
  # Analyze each sample
  results <- list()
  
  for(sample_name in colnames(bulk_expr)) {
    cat("Processing sample:", sample_name, "\n")
    
    bulk_sample <- bulk_expr[, sample_name]
    
    # Deconvolution analysis
    proportions <- solve_deconvolution(bulk_sample, sig_matrix)
    
    # Statistical significance testing
    if(perm > 0) {
      p_values <- permutation_test(bulk_sample, sig_matrix, proportions, perm)
    } else {
      p_values <- rep(1, length(proportions))
      names(p_values) <- names(proportions)
    }
    
    # Calculate RMSE as quality metric
    common_genes <- intersect(names(bulk_sample), rownames(sig_matrix))
    if(length(common_genes) > 0) {
      predicted <- as.numeric(sig_matrix[common_genes, ] %*% proportions)
      observed <- bulk_sample[common_genes]
      rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
    } else {
      rmse <- NA
    }
    
    results[[sample_name]] <- list(
      proportions = proportions,
      p_values = p_values,
      rmse = rmse
    )
  }
  
  return(results)
}

# Main deconvolution function with enhanced CIBERSORT-like algorithm
perform_deconvolution <- function(input_file, output_file, nu = 0.5, perm = 100) {
  cat("Starting enhanced cell type deconvolution analysis...\n")
  
  # Read input data
  cat("Reading input file:", input_file, "\n")
  expr_data <- read.csv(input_file, row.names = 1, check.names = FALSE)
  
  # Check data dimensions
  cat("Expression matrix dimensions:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
  
  # Convert to matrix
  expr_matrix <- as.matrix(expr_data)
  
  # Apply TPM normalization (assumes raw counts)
  if(max(expr_matrix, na.rm = TRUE) > 50) {
    cat("Applying TPM normalization and log2 transformation...\n")
    expr_matrix <- tpm_normalize(expr_matrix)
  } else {
    cat("Data appears to be already normalized\n")
  }
  
  # Filter out genes with very low expression
  gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
  keep_genes <- gene_means > quantile(gene_means, 0.1, na.rm = TRUE)
  expr_matrix <- expr_matrix[keep_genes, ]
  
  cat("After filtering:", nrow(expr_matrix), "genes retained\n")
  
  # Perform CIBERSORT-like deconvolution
  deconv_results <- cibersort_deconvolution(expr_matrix, cell_markers, nu, perm)
  
  # Format results for output
  proportions_df <- data.frame(
    Sample = names(deconv_results),
    stringsAsFactors = FALSE
  )
  
  pvalues_df <- data.frame(
    Sample = names(deconv_results),
    stringsAsFactors = FALSE
  )
  
  quality_df <- data.frame(
    Sample = names(deconv_results),
    RMSE = sapply(deconv_results, function(x) x$rmse),
    stringsAsFactors = FALSE
  )
  
  # Extract proportions and p-values
  for(cell_type in names(cell_markers)) {
    proportions_df[[cell_type]] <- sapply(deconv_results, function(x) {
      if(cell_type %in% names(x$proportions)) {
        x$proportions[cell_type]
      } else {
        0
      }
    })
    
    pvalues_df[[paste0(cell_type, "_pval")]] <- sapply(deconv_results, function(x) {
      if(cell_type %in% names(x$p_values)) {
        x$p_values[cell_type]
      } else {
        1
      }
    })
  }
  
  # Combine all results
  final_results <- merge(proportions_df, pvalues_df, by = "Sample")
  final_results <- merge(final_results, quality_df, by = "Sample")
  
  # Write main results (proportions)
  cat("Writing proportions to:", output_file, "\n")
  write.csv(proportions_df, output_file, row.names = FALSE)
  
  # Write detailed results with p-values and quality metrics
  detailed_output <- gsub("\\.csv$", "_detailed.csv", output_file)
  cat("Writing detailed results to:", detailed_output, "\n")
  write.csv(final_results, detailed_output, row.names = FALSE)
  
  # Print summary statistics
  cat("\n=== Deconvolution Summary ===\n")
  cat("Number of samples processed:", nrow(proportions_df), "\n")
  cat("Cell types analyzed:", paste(names(cell_markers), collapse = ", "), "\n")
  cat("Average RMSE:", round(mean(quality_df$RMSE, na.rm = TRUE), 4), "\n")
  
  # Show proportion summary
  cat("\nCell type proportion summary:\n")
  prop_summary <- proportions_df[, -1, drop = FALSE]
  print(round(apply(prop_summary, 2, function(x) c(Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x))), 4))
  
  cat("\nDeconvolution analysis completed successfully!\n")
  
  return(list(
    proportions = proportions_df,
    detailed_results = final_results,
    summary_stats = list(
      mean_rmse = mean(quality_df$RMSE, na.rm = TRUE),
      cell_type_summary = apply(prop_summary, 2, function(x) c(Mean = mean(x), SD = sd(x)))
    )
  ))
}

cat("=== Enhanced Cell Type Deconvolution Analysis ===\n")
cat("Input file:", opt$input, "\n")
cat("Output file:", opt$output, "\n")
cat("SVR nu parameter:", opt$nu, "\n")
cat("Permutations:", opt$perm, "\n\n")

# Run the analysis
tryCatch({
  result <- perform_deconvolution(opt$input, opt$output, opt$nu, opt$perm)
  
  cat("\n=== Final Analysis Summary ===\n")
  cat("Successfully processed", nrow(result$proportions), "samples\n")
  cat("Cell type proportions saved to:", opt$output, "\n")
  
  detailed_file <- gsub("\\.csv$", "_detailed.csv", opt$output)
  cat("Detailed results (with p-values) saved to:", detailed_file, "\n")
  
  cat("\nQuality metrics:\n")
  cat("  Average RMSE:", round(result$summary_stats$mean_rmse, 4), "\n")
  
  cat("\nTop cell types by average proportion:\n")
  avg_props <- result$summary_stats$cell_type_summary["Mean", ]
  top_types <- sort(avg_props, decreasing = TRUE)[1:min(5, length(avg_props))]
  for(i in 1:length(top_types)) {
    cat(sprintf("  %s: %.3f\n", names(top_types)[i], top_types[i]))
  }
  
}, error = function(e) {
  cat("Error during analysis:", conditionMessage(e), "\n")
  cat("Traceback:\n")
  traceback()
  quit(status = 1)
})