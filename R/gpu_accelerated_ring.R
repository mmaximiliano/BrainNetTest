# GPU-Accelerated Ring Graph Generation
# Requires: torch package for GPU acceleration

#' Check if GPU is available and return device
#' @return torch device (cuda or cpu)
get_torch_device <- function() {
  if (requireNamespace("torch", quietly = TRUE)) {
    if (torch::cuda_is_available()) {
      return(torch::torch_device("cuda"))
    }
  }
  return(NULL)
}

#' GPU-accelerated ring distance computation
#' @param N Number of nodes
#' @param device Torch device
#' @return Distance matrix as torch tensor
compute_ring_distances_gpu <- function(N, device) {
  # Create index tensors
  i_idx <- torch::torch_arange(0, N-1, device = device)$unsqueeze(2)
  j_idx <- torch::torch_arange(0, N-1, device = device)$unsqueeze(1)
  
  # Compute differences
  diff <- torch::torch_abs(i_idx - j_idx)
  
  # Ring distance: min(diff, N - diff)
  ring_dist <- torch::torch_min(diff, N - diff)
  
  return(ring_dist)
}

#' GPU-accelerated probability matrix computation
#' @param D Distance matrix (torch tensor)
#' @param lambda Decay parameter
#' @param perturbed_mask Boolean mask for perturbed edges
#' @param lambda_alt Alternative lambda for perturbed edges
#' @param p_const Constant probability for perturbed edges
#' @param perturbation_type Type of perturbation
#' @return Probability matrix as torch tensor
compute_probabilities_gpu <- function(D, lambda, perturbed_mask = NULL, 
                                    lambda_alt = NULL, p_const = NULL,
                                    perturbation_type = "none") {
  # Base probabilities
  P <- torch::torch_exp(-lambda * D)
  
  # Apply perturbations if needed
  if (!is.null(perturbed_mask) && perturbation_type != "none") {
    if (perturbation_type %in% c("lambda_half", "lambda_double")) {
      P_perturbed <- torch::torch_exp(-lambda_alt * D)
      P <- torch::torch_where(perturbed_mask, P_perturbed, P)
    } else if (perturbation_type %in% c("const_high", "const_low")) {
      P <- torch::torch_where(perturbed_mask, p_const, P)
    }
  }
  
  return(P)
}

#' Batch generate ring graphs on GPU
#' @param n_graphs Number of graphs to generate
#' @param N Number of nodes
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation
#' @param lambda_alt Alternative lambda
#' @param p_const Constant probability
#' @param device Torch device
#' @return List of adjacency matrices
batch_generate_ring_graphs_gpu <- function(n_graphs, N, lambda, 
                                         perturbed_nodes = c(),
                                         perturbation_type = "none",
                                         lambda_alt = NULL, 
                                         p_const = NULL,
                                         device) {
  
  # Compute distance matrix once
  D <- compute_ring_distances_gpu(N, device)
  
  # Create perturbation mask if needed
  perturbed_mask <- NULL
  if (length(perturbed_nodes) > 0 && perturbation_type != "none") {
    mask <- torch::torch_zeros(N, N, dtype = torch::torch_bool(), device = device)
    for (node in perturbed_nodes - 1) {  # Convert to 0-based indexing
      mask[node, ] <- TRUE
      mask[, node] <- TRUE
    }
    perturbed_mask <- mask
  }
  
  # Compute probability matrix
  P <- compute_probabilities_gpu(D, lambda, perturbed_mask, 
                               lambda_alt, p_const, perturbation_type)
  
  # Get upper triangular indices
  upper_tri <- torch::torch_triu(torch::torch_ones(N, N, device = device), diagonal = 1)
  upper_tri_mask <- upper_tri$to(dtype = torch::torch_bool())
  
  # Generate all graphs in batch
  graphs <- list()
  
  # Process in mini-batches to avoid memory issues
  batch_size <- min(100, n_graphs)
  n_batches <- ceiling(n_graphs / batch_size)
  
  for (b in 1:n_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, n_graphs)
    current_batch_size <- end_idx - start_idx + 1
    
    # Generate random values for entire batch
    rand_vals <- torch::torch_rand(current_batch_size, N, N, device = device)
    
    # Create edges based on probabilities
    edge_batch <- (rand_vals < P$unsqueeze(1))$to(dtype = torch::torch_float32())
    
    # Apply upper triangular mask and symmetrize
    edge_batch <- edge_batch * upper_tri_mask$unsqueeze(1)
    edge_batch <- edge_batch + edge_batch$transpose(-1, -2)
    
    # Convert to R matrices
    edge_batch_cpu <- edge_batch$cpu()$numpy()
    
    for (i in 1:current_batch_size) {
      graphs[[start_idx + i - 1]] <- edge_batch_cpu[i,,]
    }
  }
  
  return(graphs)
}
