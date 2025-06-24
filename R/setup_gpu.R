#' Setup GPU acceleration for BrainNetTest
#' 
#' This function helps install and configure GPU support
#' @export
setup_gpu_acceleration <- function() {
  # Check if torch is installed
  if (!requireNamespace("torch", quietly = TRUE)) {
    message("Installing torch package for GPU acceleration...")
    install.packages("torch")
  }
  
  # Install torch with CUDA support
  if (requireNamespace("torch", quietly = TRUE)) {
    # Check if torch is already installed with CUDA
    if (!torch::cuda_is_available()) {
      message("Installing torch with CUDA support...")
      torch::install_torch(reinstall = FALSE)
    }
    
    # Verify GPU is available
    if (torch::cuda_is_available()) {
      gpu_info <- torch::cuda_device_properties(0)
      message(sprintf("GPU detected: %s", gpu_info$name))
      message(sprintf("GPU memory: %.2f GB", gpu_info$total_memory / 1e9))
      return(TRUE)
    } else {
      warning("No GPU detected. This could be due to:\n",
              "1. No compatible NVIDIA GPU on this system\n",
              "2. Missing CUDA drivers\n",
              "3. Incompatible CUDA version\n",
              "Falling back to CPU.")
      return(FALSE)
    }
  }
  
  warning("Failed to load torch package. Falling back to CPU.")
  return(FALSE)
}

#' Benchmark GPU vs CPU performance
#' @param N Number of nodes to test
#' @param n_graphs Number of graphs to generate
#' @export
benchmark_gpu_performance <- function(N = 1000, n_graphs = 100) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("torch package is required for GPU benchmarking")
  }
  
  if (!torch::cuda_is_available()) {
    stop("GPU not available for benchmarking")
  }
  
  lambda <- find_lambda(N)
  
  # CPU timing
  cpu_start <- Sys.time()
  cpu_graphs <- generate_population(n_graphs, N, lambda)
  cpu_time <- as.numeric(Sys.time() - cpu_start, units = "secs")
  
  # GPU timing
  device <- torch::torch_device("cuda")
  gpu_start <- Sys.time()
  gpu_graphs <- batch_generate_ring_graphs_gpu(n_graphs, N, lambda, device = device)
  gpu_time <- as.numeric(Sys.time() - gpu_start, units = "secs")
  
  message(sprintf("CPU time: %.2f seconds", cpu_time))
  message(sprintf("GPU time: %.2f seconds", gpu_time))
  message(sprintf("Speedup: %.2fx", cpu_time / gpu_time))
  
  return(list(cpu_time = cpu_time, gpu_time = gpu_time, speedup = cpu_time / gpu_time))
}

#' Check CUDA compatibility
#' @return List with CUDA status information
#' @export
check_cuda_compatibility <- function() {
  results <- list(
    torch_installed = requireNamespace("torch", quietly = TRUE),
    cuda_available = FALSE,
    cuda_version = NA,
    gpu_name = NA,
    gpu_memory = NA,
    driver_version = NA
  )
  
  if (results$torch_installed) {
    results$cuda_available <- torch::cuda_is_available()
    
    if (results$cuda_available) {
      # Get GPU properties
      gpu_props <- torch::cuda_device_properties(0)
      results$gpu_name <- gpu_props$name
      results$gpu_memory <- gpu_props$total_memory / 1e9  # Convert to GB
      
      # Try to get CUDA version
      tryCatch({
        results$cuda_version <- torch::cuda_version()
        results$driver_version <- torch::cuda_driver_version()
      }, error = function(e) {
        # Some torch versions may not have these functions
      })
    }
  }
  
  # Print summary
  cat("CUDA Compatibility Check:\n")
  cat("------------------------\n")
  cat(sprintf("torch installed: %s\n", ifelse(results$torch_installed, "Yes", "No")))
  cat(sprintf("CUDA available: %s\n", ifelse(results$cuda_available, "Yes", "No")))
  
  if (results$cuda_available) {
    cat(sprintf("GPU: %s\n", results$gpu_name))
    cat(sprintf("GPU memory: %.2f GB\n", results$gpu_memory))
    if (!is.na(results$cuda_version)) {
      cat(sprintf("CUDA version: %s\n", results$cuda_version))
    }
    if (!is.na(results$driver_version)) {
      cat(sprintf("Driver version: %s\n", results$driver_version))
    }
  } else if (results$torch_installed) {
    cat("\nPossible reasons for CUDA unavailability:\n")
    cat("1. No compatible NVIDIA GPU present\n")
    cat("2. CUDA drivers not installed\n")
    cat("3. torch not built with CUDA support\n")
    cat("\nTry running: torch::install_torch(reinstall = TRUE)\n")
  }
  
  return(invisible(results))
}
