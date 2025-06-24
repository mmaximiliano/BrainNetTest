# GPU Troubleshooting Script
# ──────────────────────────────────────────────────────────────────────────────

# Load required functions
tryCatch({
  source("R/utils.R")
  source("R/setup_gpu.R")
  source("R/gpu_accelerated_ring.R")
}, error = function(e) {
  message("Error loading required files: ", e$message)
})

# Check system information
cat("=== System Information ===\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
cat(sprintf("OS: %s\n", Sys.info()["sysname"]))
cat(sprintf("Windows version: %s\n", if(Sys.info()["sysname"] == "Windows") 
                                      paste(Sys.info()["release"], Sys.info()["version"]) else "N/A"))
cat("\n")

# Check for torch package
cat("=== Torch Package ===\n")
torch_available <- requireNamespace("torch", quietly = TRUE)
cat(sprintf("torch package installed: %s\n", if(torch_available) "Yes" else "No"))

if (torch_available) {
  # Get torch version
  tryCatch({
    torch_version <- utils::packageVersion("torch")
    cat(sprintf("torch version: %s\n", torch_version))
  }, error = function(e) {
    cat("Could not determine torch version\n")
  })
  
  # Check CUDA availability
  cuda_available <- FALSE
  tryCatch({
    cuda_available <- torch::cuda_is_available()
    cat(sprintf("CUDA available: %s\n", if(cuda_available) "Yes" else "No"))
    
    if (cuda_available) {
      device_count <- torch::cuda_device_count()
      cat(sprintf("CUDA device count: %d\n", device_count))
      
      # Get details of the first device
      device_props <- torch::cuda_device_properties(0)
      cat(sprintf("GPU name: %s\n", device_props$name))
      cat(sprintf("GPU memory: %.2f GB\n", device_props$total_memory / 1e9))
      cat(sprintf("CUDA capability: %d.%d\n", device_props$major, device_props$minor))
    }
  }, error = function(e) {
    cat("Error checking CUDA availability: ", e$message, "\n")
  })
  
  if (!cuda_available) {
    cat("\n=== NVIDIA GPU Detection ===\n")
    # On Windows, try to detect NVIDIA GPUs using the system
    if (Sys.info()["sysname"] == "Windows") {
      tryCatch({
        nvidia_smi <- system("nvidia-smi", intern = TRUE, ignore.stderr = TRUE)
        cat("nvidia-smi output:\n")
        cat(paste(nvidia_smi[1:10], collapse = "\n"), "\n")
        cat("...\n")
      }, error = function(e) {
        cat("nvidia-smi not found or failed. This suggests NVIDIA drivers may not be installed.\n")
      })
    }
    
    cat("\n=== Troubleshooting Steps ===\n")
    cat("1. Make sure you have an NVIDIA GPU that supports CUDA\n")
    cat("2. Install the latest NVIDIA drivers from https://www.nvidia.com/Download/index.aspx\n")
    cat("3. Reinstall torch with CUDA support by running 'torch::install_torch(reinstall = TRUE)'\n")
    cat("4. If on Windows, make sure your GPU is not being exclusively used by another application\n")
  }
}

# Try to fix the issues
cat("\n=== Attempting to Fix Issues ===\n")

if (!torch_available) {
  cat("Installing torch package...\n")
  install.packages("torch")
} else if (!cuda_available) {
  cat("Reinstalling torch with CUDA support...\n")
  torch::install_torch(reinstall = TRUE)
}

# Try GPU setup again
cat("\n=== Testing GPU Setup ===\n")
tryCatch({
  gpu_available <- setup_gpu_acceleration()
  cat(sprintf("GPU setup successful: %s\n", if(gpu_available) "Yes" else "No"))
  
  if (gpu_available) {
    cat("\n=== GPU Test ===\n")
    cat("Running small GPU benchmark...\n")
    test_device <- torch::torch_device("cuda")
    # Create and manipulate a small tensor to verify operation
    test_tensor <- torch::torch_rand(c(1000, 1000), device = test_device)
    result_tensor <- test_tensor %*% test_tensor
    cat(sprintf("GPU operation successful. Tensor shape: [%d, %d]\n", 
                result_tensor$shape[1], result_tensor$shape[2]))
  } else {
    cat("GPU setup failed. Please review the troubleshooting steps above.\n")
  }
}, error = function(e) {
  cat("Error during GPU setup: ", e$message, "\n")
})

cat("\n=== Verification Complete ===\n")
cat("If GPU issues persist, please ensure your system meets these requirements:\n")
cat("1. Compatible NVIDIA GPU\n")
cat("2. Latest NVIDIA drivers installed\n")
cat("3. CUDA Toolkit compatible with your GPU\n")
cat("4. For Windows users: check if Windows is using the integrated GPU instead of NVIDIA GPU\n")
cat("5. For Windows users: check Windows Graphics Settings to prefer High Performance NVIDIA GPU\n")
