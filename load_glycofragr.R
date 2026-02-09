# Simple glycofragR loader - avoids devtools installation issues
# Just loads glycofragR directly from source

message("Loading glycofragR from source...")

# Get the current working directory
pkg_root <- getwd()
glycofragr_path <- file.path(pkg_root, "glycofragR")

if (!dir.exists(glycofragr_path)) {
  stop("glycofragR directory not found at: ", glycofragr_path)
}

message(paste("Package path:", glycofragr_path))

# Try using devtools::load_all first
tryCatch({
  # Try to load devtools without requiring it (avoid namespace issues)
  if (requireNamespace("devtools", quietly = TRUE)) {
    message("Using devtools::load_all...")
    devtools::load_all(glycofragr_path, quiet = TRUE)
    message("✓ glycofragR loaded via devtools")
  } else {
    stop("devtools not available")
  }
}, error = function(e) {
  message("Note: devtools not available, using manual loading...")
  
  # Manual fallback: Source all R files and set up namespace
  r_dir <- file.path(glycofragr_path, "R")
  
  if (!dir.exists(r_dir)) {
    stop("R directory not found in glycofragR package")
  }
  
  # Source all .R files in order
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  r_files <- sort(r_files)  # Ensure consistent order
  
  message(paste("Sourcing", length(r_files), "R files..."))
  
  for (f in r_files) {
    message(paste("  -", basename(f)))
    source(f)
  }
  
  message("✓ glycofragR loaded manually")
})

# Verify glycofragR functions are available
if (exists("glycan_analysis_new")) {
  message("✓ glycofragR functions available!")
} else {
  stop("glycofragR functions not loaded - something went wrong")
}

message("\nYou can now use glycofragR:")
message("  analysis <- glycan_analysis_new('4502', 'permethylated_reduced')")
message("  df <- glycan_analysis_theoretical_df(analysis)")
