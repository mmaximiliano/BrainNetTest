# GPU-Accelerated Ring Graph Generation
# Requires: torch package for GPU acceleration

#' Check if GPU is available and return device
#' @return torch device (cuda or cpu)
get_torch_device <- function() {
  # First check if .gpu_env exists and is enabled
  if (exists(".gpu_env", mode = "environment") && 
      !is.null(.gpu_env$enabled) && 
      .gpu_env$enabled &&
      !is.null(.gpu_env$device)) {
    return(.gpu_env$device)
  }
  
  # Otherwise check if torch is available and CUDA is available
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

#' Batch generate ring graphs on GPU with optimized memory usage
#' @param n_graphs Number of graphs to generate
#' @param N Number of nodes per graph
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation
#' @param lambda_alt Alternative lambda for perturbed edges
#' @param p_const Constant probability for perturbed edges
#' @param device Torch device
#' @param batch_size Number of graphs to process in each batch
#' @return List of adjacency matrices
batch_generate_ring_graphs_gpu <- function(n_graphs, N, lambda, perturbed_nodes = c(),
                                         perturbation_type = "none",
                                         lambda_alt = NULL, p_const = NULL, 
                                         device = NULL, batch_size = NULL) {
  
  if (is.null(device)) {
    device <- get_torch_device()
  }
  
  # Determine optimal batch size based on available GPU memory
  if (is.null(batch_size)) {
    # Estimate memory requirement per graph
    mem_per_graph <- N * N * 4 * 3  # float32, 3 matrices (D, P, adj)
    
    # Try to get available GPU memory
    available_memory <- tryCatch({
      # Get total memory and currently allocated
      props <- torch::cuda_get_device_properties(0)
      total_mem <- props$total_memory
      
      # Reserve 20% for overhead
      usable_mem <- total_mem * 0.8
      
      # Calculate batch size
      max_batch <- floor(usable_mem / mem_per_graph)
      
      # Cap at reasonable limits
      min(max_batch, n_graphs, 1000)
    }, error = function(e) {
      # Fallback batch sizes based on N
      if (N <= 100) 100
      else if (N <= 1000) 50
      else 10
    })
    
    batch_size <- max(1, available_memory)
  }
  
  # Process in batches to manage memory
  all_graphs <- list()
  n_batches <- ceiling(n_graphs / batch_size)
  
  for (batch_idx in 1:n_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, n_graphs)
    current_batch_size <- end_idx - start_idx + 1
    
    # Generate batch of graphs
    batch_graphs <- generate_ring_graph_batch_gpu(
      current_batch_size, N, lambda, perturbed_nodes,
      perturbation_type, lambda_alt, p_const, device
    )
    
    # Convert to list of matrices
    for (i in 1:current_batch_size) {
      all_graphs[[start_idx + i - 1]] <- batch_graphs[i,,]$cpu()$to(dtype = torch::torch_float())$numpy()
    }
    
    # Clear GPU cache periodically
    if (batch_idx %% 5 == 0) {
      torch::cuda_empty_cache()
    }
  }
  
  return(all_graphs)
}

#' Batch generate ring graphs on GPU with CUDA stream support
#' @param n_graphs Number of graphs to generate
#' @param N Number of nodes per graph
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation
#' @param lambda_alt Alternative lambda for perturbed edges
#' @param p_const Constant probability for perturbed edges
#' @param device Torch device to use
#' @param batch_size Batch size for generation (NULL for auto)
#' @param cuda_stream Optional CUDA stream for concurrent execution
#' @return List of adjacency matrices
batch_generate_ring_graphs_gpu <- function(n_graphs, N, lambda, perturbed_nodes = c(),
                                         perturbation_type = "none",
                                         lambda_alt = NULL, p_const = NULL,
                                         device = NULL, batch_size = NULL,
                                         cuda_stream = NULL) {
  
  if (is.null(device)) {
    device <- get_torch_device()
  }
  
  # Set CUDA stream context if provided
  if (!is.null(cuda_stream) && !is.null(torch$cuda$stream)) {
    stream_context <- torch$cuda$stream(cuda_stream)
  }
  
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
  
  # Ensure operations are synchronized if using streams
  if (!is.null(cuda_stream)) {
    torch$cuda$synchronize()
  }
  
  return(graphs)
}

#' Create multiple CUDA streams for concurrent GPU operations
#' @param n_streams Number of CUDA streams to create
#' @return List of CUDA streams or NULL if not available
create_cuda_streams <- function(n_streams = 4) {
  if (!exists("torch") || is.null(torch$cuda$Stream)) {
    return(NULL)
  }
  
  streams <- list()
  for (i in 1:n_streams) {
    tryCatch({
      streams[[i]] <- torch$cuda$Stream()
    }, error = function(e) {
      warning(sprintf("Failed to create CUDA stream %d: %s", i, e$message))
      return(NULL)
    })
  }
  
  return(streams)
}

#' Get optimal batch size based on GPU memory and network size
#' @param N Number of nodes
#' @param available_memory Available GPU memory in bytes
#' @param n_concurrent Number of concurrent operations
#' @return Optimal batch size
get_optimal_batch_size_concurrent <- function(N, available_memory = NULL, n_concurrent = 1) {
  if (is.null(available_memory) && exists("torch") && torch$cuda$is_available()) {
    # Get available memory accounting for concurrent operations
    total_memory <- torch$cuda$get_device_properties(0)$total_memory
    allocated_memory <- torch$cuda$memory_allocated()
    available_memory <- (total_memory - allocated_memory) / n_concurrent
  }
  
  # Adjust batch sizes for concurrent execution
  bytes_per_element <- 4  # float32
  elements_per_graph <- N * N
  bytes_per_graph <- elements_per_graph * bytes_per_element
  
  # Reserve memory for operations and concurrent execution
  memory_overhead_factor <- 4 * n_concurrent
  usable_memory <- available_memory / memory_overhead_factor
  
  # Calculate batch size with safety margin
  batch_size <- floor(usable_memory / bytes_per_graph)
  
  # Apply network-size specific limits adjusted for concurrency
  if (N <= 100) {
    batch_size <- min(batch_size, floor(1000 / n_concurrent))
  } else if (N <= 1000) {
    batch_size <- min(batch_size, floor(100 / n_concurrent))
  } else if (N <= 10000) {
    batch_size <- min(batch_size, floor(10 / n_concurrent))
  } else {
    batch_size <- min(batch_size, floor(5 / n_concurrent))
  }
  
  return(max(1, batch_size))
}

#' GPU-accelerated bootstrap sampling for edge statistics
#' @param edge_stats Matrix of edge statistics (n_edges x n_populations)
#' @param n_bootstrap Number of bootstrap samples
#' @param device Torch device
#' @return List with bootstrap results
gpu_bootstrap_edge_stats <- function(edge_stats, n_bootstrap, device = NULL) {
  
  if (is.null(device)) {
    device <- get_torch_device()
  }
  
  # Convert to tensor
  stats_tensor <- torch::torch_tensor(edge_stats, device = device)
  n_samples <- nrow(edge_stats)
  n_pops <- ncol(edge_stats)
  
  # Generate all bootstrap indices at once
  bootstrap_indices <- torch::torch_randint(
    0, n_samples, 
    size = c(n_bootstrap, n_samples),
    dtype = torch::torch_long(),
    device = device
  )
  
  # Gather bootstrap samples
  bootstrap_samples <- list()
  for (pop in 1:n_pops) {
    pop_stats <- stats_tensor[, pop]
    # Index into population statistics
    bootstrap_samples[[pop]] <- pop_stats$index_select(0, bootstrap_indices$view(-1))$view(c(n_bootstrap, n_samples))
  }
  
  # Compute bootstrap means
  bootstrap_means <- list()
  for (pop in 1:n_pops) {
    bootstrap_means[[pop]] <- bootstrap_samples[[pop]]$mean(dim = 2)
  }
  
  # Compute test statistics
  if (n_pops == 2) {
    # Two-sample test
    diff_means <- bootstrap_means[[1]] - bootstrap_means[[2]]
    test_stats <- diff_means$cpu()$numpy()
  } else {
    # Multi-sample test (simplified)
    # Would need to implement proper multi-sample test
    test_stats <- NULL
  }
  
  return(list(
    bootstrap_means = bootstrap_means,
    test_statistics = test_stats
  ))
}

#' GPU-accelerated computation of edge presence across populations
#' @param populations List of population graphs
#' @param device Torch device
#' @return Matrix of edge presence indicators
gpu_compute_edge_presence <- function(populations, device = NULL) {
  
  if (is.null(device)) {
    device <- get_torch_device()
  }
  
  # Get dimensions
  n_pops <- length(populations)
  n_graphs_per_pop <- length(populations[[1]])
  N <- nrow(populations[[1]][[1]])
  
  # Stack all graphs into tensors per population
  pop_tensors <- list()
  for (p in 1:n_pops) {
    # Convert list of matrices to 3D tensor
    graphs_array <- array(unlist(populations[[p]]), dim = c(N, N, n_graphs_per_pop))
    pop_tensors[[p]] <- torch::torch_tensor(graphs_array, device = device)$permute(c(3, 1, 2))
  }
  
  # Get upper triangular indices
  upper_indices <- torch::torch_triu_indices(N, N, offset = 1L, device = device)
  n_edges <- upper_indices$size(2)
  
  # Extract edge values for all graphs
  edge_presence <- matrix(0, nrow = n_edges, ncol = n_pops)
  
  for (p in 1:n_pops) {
    # Extract upper triangular values for all graphs in population
    edge_vals <- pop_tensors[[p]]$index_select(
      2, upper_indices[1,]
    )$index_select(
      3, upper_indices[2,]
    )
    
    # Sum across graphs to get presence count
    edge_counts <- edge_vals$sum(dim = 1)$cpu()$numpy()
    edge_presence[, p] <- edge_counts
  }
  
  return(edge_presence)
}

#' GPU-accelerated Fisher's exact test for edge differences
#' @param edge_presence Matrix of edge presence counts
#' @param device Torch device
#' @return Vector of p-values
gpu_fisher_test_edges <- function(edge_presence, device = NULL) {
  
  if (is.null(device)) {
    device <- get_torch_device()
  }
  
  # For now, fall back to CPU for Fisher's exact test
  # as it requires special functions not easily implemented in torch
  
  # However, we can use GPU to quickly identify candidate edges
  # based on large differences in presence
  
  presence_tensor <- torch::torch_tensor(edge_presence, device = device)
  n_graphs_per_pop <- max(edge_presence)  # Assuming same size populations
  
  # Compute proportions
  props <- presence_tensor / n_graphs_per_pop
  
  # For two populations, compute absolute difference
  if (ncol(edge_presence) == 2) {
    prop_diff <- torch::torch_abs(props[, 1] - props[, 2])
    
    # Identify edges with large differences (candidates for testing)
    threshold <- 0.1  # 10% difference threshold
    candidate_edges <- torch::torch_where(prop_diff > threshold)[1]$cpu()$numpy() + 1
    
    # Only test candidate edges on CPU
    p_values <- rep(1, nrow(edge_presence))
    
    if (length(candidate_edges) > 0) {
      for (edge_idx in candidate_edges) {
        # Perform Fisher's exact test on CPU
        test_result <- fisher.test(matrix(c(
          edge_presence[edge_idx, 1],
          n_graphs_per_pop - edge_presence[edge_idx, 1],
          edge_presence[edge_idx, 2],
          n_graphs_per_pop - edge_presence[edge_idx, 2]
        ), nrow = 2))
        
        p_values[edge_idx] <- test_result$p.value
      }
    }
    
    return(p_values)
  }
  
  # For multiple populations, would need to implement appropriate test
  return(rep(1, nrow(edge_presence)))
}
