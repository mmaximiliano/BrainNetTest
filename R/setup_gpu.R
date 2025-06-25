#' Check CUDA compatibility
#' @return List with compatibility information
check_cuda_compatibility <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    return(list(
      compatible = FALSE,
      cuda_available = FALSE,
      device_count = 0,
      devices = list(),
      message = "torch package not installed"
    ))
  }
  
  # Check if CUDA is available
  cuda_available <- tryCatch({
    torch::cuda_is_available()
  }, error = function(e) {
    FALSE
  })
  
  if (!cuda_available) {
    return(list(
      compatible = FALSE,
      cuda_available = FALSE,
      device_count = 0,
      devices = list(),
      message = "CUDA not available"
    ))
  }
  
  # Get device count
  device_count <- tryCatch({
    torch::cuda_device_count()
  }, error = function(e) {
    0
  })
  
  if (device_count == 0) {
    return(list(
      compatible = FALSE,
      cuda_available = TRUE,
      device_count = 0,
      devices = list(),
      message = "No CUDA devices found"
    ))
  }
  
  # Get device properties using available functions
  devices <- list()
  for (i in seq_len(device_count)) {
    device_info <- tryCatch({
      # Get compute capability using available function
      capability <- NULL
      if (exists("cuda_get_device_capability", where = asNamespace("torch"))) {
        cap <- torch::cuda_get_device_capability(i - 1)  # 0-indexed
        if (!is.null(cap) && length(cap) >= 2) {
          capability <- paste0(cap[1], ".", cap[2])
        }
      }
      
      # Get runtime version
      runtime_version <- NULL
      if (exists("cuda_runtime_version", where = asNamespace("torch"))) {
        runtime_version <- torch::cuda_runtime_version()
      }
      
      # Try to get memory info
      memory_gb <- NA
      if (exists("cuda_memory_stats", where = asNamespace("torch"))) {
        tryCatch({
          # Set current device
          old_device <- torch::cuda_current_device()
          torch::cuda_set_device(i)  # 1-indexed in torch
          
          mem_stats <- torch::cuda_memory_stats(device = i - 1)
          if (!is.null(mem_stats) && "allocated_bytes.all.current" %in% names(mem_stats)) {
            # Get total memory from reserved or use a reasonable estimate
            if ("reserved_bytes.all.peak" %in% names(mem_stats)) {
              total_memory <- mem_stats$reserved_bytes.all.peak
            } else {
              # Estimate based on typical GPU memory sizes
              total_memory <- 8 * 1024^3  # Default to 8GB
            }
            memory_gb <- round(total_memory / 1024^3, 2)
          }
          
          # Restore device
          torch::cuda_set_device(old_device + 1)
        }, error = function(e) {
          # Memory stats failed, continue
        })
      }
      
      list(
        index = i - 1,
        name = paste("CUDA Device", i - 1),
        compute_capability = ifelse(!is.null(capability), capability, "Unknown"),
        memory_gb = memory_gb,
        runtime_version = runtime_version
      )
    }, error = function(e) {
      # Minimal device info if properties can't be retrieved
      list(
        index = i - 1,
        name = paste("CUDA Device", i - 1),
        compute_capability = "Unknown",
        memory_gb = NA,
        runtime_version = NA
      )
    })
    
    devices[[i]] <- device_info
  }
  
  return(list(
    compatible = TRUE,
    cuda_available = TRUE,
    device_count = device_count,
    devices = devices,
    message = paste("Found", device_count, "CUDA device(s)")
  ))
}

# GPU Setup and Acceleration for BrainNetTest
# ============================================

# Initialize GPU environment
.gpu_env <- new.env(parent = emptyenv())
.gpu_env$enabled <- FALSE
.gpu_env$device <- NULL

#' Setup GPU acceleration for BrainNetTest
#' @param device_index GPU device index to use (default: 0)
#' @param verbose Print setup information
#' @return TRUE if GPU setup successful, FALSE otherwise
#' @export
setup_gpu_acceleration <- function(device_index = 0, verbose = TRUE) {
  # Check compatibility
  compat <- check_cuda_compatibility()
  
  if (!compat$compatible) {
    if (verbose) {
      message("GPU acceleration not available: ", compat$message)
    }
    .gpu_env$device <- NULL
    .gpu_env$enabled <- FALSE
    return(FALSE)
  }
  
  # Set device
  tryCatch({
    if (device_index >= compat$device_count) {
      device_index <- 0
      if (verbose) {
        message("Requested device index ", device_index, " not available. Using device 0.")
      }
    }
    
    # Create device string
    device_str <- paste0("cuda:", device_index)
    
    # Try to create torch device directly without cuda_set_device
    .gpu_env$device <- torch::torch_device(device_str)
    .gpu_env$enabled <- TRUE
    
    if (verbose) {
      device_info <- compat$devices[[device_index + 1]]
      message("GPU acceleration enabled")
      message("Using device: ", device_info$name)
      if (!is.na(device_info$compute_capability) && device_info$compute_capability != "Unknown") {
        message("Compute capability: ", device_info$compute_capability)
      }
      if (!is.na(device_info$memory_gb)) {
        message("Memory: ", device_info$memory_gb, " GB")
      }
      if (!is.null(device_info$runtime_version)) {
        message("CUDA Runtime version: ", device_info$runtime_version)
      }
    }
    
    # Warm up GPU
    tryCatch({
      x <- torch::torch_randn(c(100, 100), device = .gpu_env$device)
      y <- torch::torch_matmul(x, x)
      rm(x, y)
      torch::cuda_empty_cache()
      
      if (verbose) {
        message("GPU warmup completed")
      }
    }, error = function(e) {
      # GPU warmup failed, but don't fail the whole setup
      if (verbose) {
        message("GPU warmup skipped")
      }
    })
    
    return(TRUE)
    
  }, error = function(e) {
    if (verbose) {
      message("Failed to setup GPU: ", e$message)
    }
    .gpu_env$device <- NULL
    .gpu_env$enabled <- FALSE
    return(FALSE)
  })
}

#' Check if GPU is enabled
#' @return TRUE if GPU is enabled, FALSE otherwise
#' @export
is_gpu_enabled <- function() {
  exists(".gpu_env", mode = "environment") && 
    !is.null(.gpu_env$enabled) && 
    .gpu_env$enabled
}

#' Get current GPU device
#' @return torch device object or NULL
#' @export
get_gpu_device <- function() {
  if (is_gpu_enabled()) {
    return(.gpu_env$device)
  }
  return(NULL)
}

#' Benchmark GPU performance
#' @param sizes Vector of matrix sizes to test
#' @param verbose Print benchmark results
#' @return Data frame with benchmark results
#' @export
benchmark_gpu_performance <- function(sizes = c(100, 500, 1000, 2500, 5000), 
                                    verbose = TRUE) {
  if (!is_gpu_enabled()) {
    if (verbose) message("GPU not enabled. Run setup_gpu_acceleration() first.")
    return(NULL)
  }
  
  results <- data.frame(
    size = integer(),
    operation = character(),
    device = character(),
    time_seconds = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size in sizes) {
    if (verbose) cat(sprintf("Benchmarking size %dx%d...\n", size, size))
    
    # Test matrix multiplication
    for (device_type in c("cpu", "cuda")) {
      device <- if (device_type == "cuda") .gpu_env$device else "cpu"
      
      tryCatch({
        # Create random matrices
        x <- torch::torch_randn(c(size, size), device = device)
        y <- torch::torch_randn(c(size, size), device = device)
        
        # Ensure tensors are on device
        if (device_type == "cuda") {
          torch::cuda_synchronize()
        }
        
        # Time the operation
        start_time <- Sys.time()
        z <- torch::torch_matmul(x, y)
        if (device_type == "cuda") {
          torch::cuda_synchronize()
        }
        end_time <- Sys.time()
        
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        
        results <- rbind(results, data.frame(
          size = size,
          operation = "matmul",
          device = device_type,
          time_seconds = time_taken,
          stringsAsFactors = FALSE
        ))
        
        # Clean up
        rm(x, y, z)
        if (device_type == "cuda") {
          torch::cuda_empty_cache()
        }
        
      }, error = function(e) {
        if (verbose) cat(sprintf("  Failed for %s: %s\n", device_type, e$message))
      })
    }
  }
  
  if (verbose && nrow(results) > 0) {
    cat("\nBenchmark Results:\n")
    cat("==================\n")
    
    # Calculate speedups
    for (s in unique(results$size)) {
      size_results <- results[results$size == s & results$operation == "matmul", ]
      if (nrow(size_results) == 2) {
        cpu_time <- size_results$time_seconds[size_results$device == "cpu"]
        gpu_time <- size_results$time_seconds[size_results$device == "cuda"]
        speedup <- cpu_time / gpu_time
        
        cat(sprintf("Size %dx%d: CPU=%.4fs, GPU=%.4fs, Speedup=%.1fx\n",
                    s, s, cpu_time, gpu_time, speedup))
      }
    }
  }
  
  return(results)
}
