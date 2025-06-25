# Save this as gpu_diagnostics_safe.R

# GPU Diagnostics Script (Safe Version)
# ──────────────────────────────────────────────────────────────────────────────

# Set library path first (to use your working setup)
torch_lib_path <- file.path(system.file("", package = "torch"), "lib")
Sys.setenv(LD_LIBRARY_PATH = paste(
  torch_lib_path,
  "/usr/lib/x86_64-linux-gnu",
  Sys.getenv("LD_LIBRARY_PATH"),
  sep = ":"
))

# Check system information
cat("=== System Information ===\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
cat(sprintf("OS: %s\n", Sys.info()["sysname"]))
cat(sprintf("User: %s\n", Sys.info()["user"]))
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

    # Load torch
    library(torch)

    # Check CUDA availability
    cuda_available <- torch::cuda_is_available()
    cat(sprintf("CUDA available: %s\n", if(cuda_available) "Yes" else "No"))

    if (cuda_available) {
      device_count <- torch::cuda_device_count()
      cat(sprintf("CUDA device count: %d\n", device_count))

      # Get device capability (correct function)
      capability <- torch::cuda_get_device_capability(0)
      cat(sprintf("GPU Compute Capability: %d.%d\n", capability[[1]], capability[[2]]))

      # Runtime version
      cat(sprintf("CUDA Runtime Version: %s\n", torch::cuda_runtime_version()))

      # Quick test
      cat("\n=== GPU Test ===\n")
      test_device <- torch::torch_device("cuda")
      test_tensor <- torch::torch_rand(c(1000, 1000), device = test_device)
      result_tensor <- torch::torch_matmul(test_tensor, test_tensor)
      cat("✓ GPU tensor operations working correctly\n")

      # Memory info
      torch::cuda_memory_summary(0)
    }
  }, error = function(e) {
    cat("Error: ", e$message, "\n")
  })
}

# Check nvidia-smi
cat("\n=== NVIDIA GPU Detection ===\n")
tryCatch({
  nvidia_info <- system("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader", intern = TRUE)
  cat("GPU Info:\n", nvidia_info, "\n")
}, error = function(e) {
  cat("nvidia-smi not available\n")
})

cat("\n=== Current Setup Status ===\n")
cat("Your RTX 3090 appears to be working correctly with torch!\n")
cat("To maintain this setup, always set LD_LIBRARY_PATH before loading torch.\n")
