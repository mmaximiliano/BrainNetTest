# Initialize GPU support when the package is loaded

.onLoad <- function(libname, pkgname) {
  # Quietly check for torch
  has_torch <- suppressWarnings(suppressMessages(
    requireNamespace("torch", quietly = TRUE)
  ))
  
  if (has_torch) {
    # Only print a message if torch is loaded
    packageStartupMessage("torch package detected. Use setup_gpu_acceleration() to enable GPU support.")
    
    # Check if CUDA is available
    has_cuda <- FALSE
    tryCatch({
      has_cuda <- torch::cuda_is_available()
      if (has_cuda) {
        device_info <- torch::cuda_device_properties(0)
        packageStartupMessage(sprintf("CUDA-capable GPU detected: %s", device_info$name))
      }
    }, error = function(e) {
      # Silently fail, don't show error messages on startup
    })
  }
  
  invisible()
}
