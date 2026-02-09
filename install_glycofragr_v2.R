# Install glycofragR package from local source

message("Step 1: Installing devtools (if needed)...")
if (!require("devtools", character.only = TRUE)) {
  message("  Installing devtools...")
  install.packages("devtools", repos = "https://cran.r-project.org", quiet = TRUE)
}

message("✓ devtools ready")

message("\nStep 2: Building and installing glycofragR...")
message("This may take a moment...")

# Get the current working directory (should be the package root)
pkg_root <- getwd()
glycofragr_path <- file.path(pkg_root, "glycofragR")

# Check if glycofragR directory exists
if (!dir.exists(glycofragr_path)) {
  stop("glycofragR directory not found at: ", glycofragr_path)
}

message(paste("  Package location:", glycofragr_path))

# Install the package
tryCatch({
  devtools::install(
    pkg = glycofragr_path, 
    force = TRUE, 
    upgrade = "never",
    dependencies = FALSE,
    quiet = TRUE
  )
  message("✓ glycofragR built and installed successfully!")
}, error = function(e) {
  message("Warning: Installation encountered an issue:")
  message(paste("  ", e$message))
  message("\nTrying alternative method...")
  
  # Try load_all as fallback
  devtools::load_all(glycofragr_path)
  message("✓ glycofragR loaded (development mode)")
})

message("\n✓ glycofragR is ready!")
message("\nYou can now use: library(glycofragR)")
message("\nNext steps:")
message("1. Run: source('test_basic.R')")
message("2. Or: source('GlycofragR_test_UPDATED.R')")
