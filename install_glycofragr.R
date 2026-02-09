# Install glycofragR package from local source

message("Step 1: Installing devtools (if needed)...")
if (!require("devtools", character.only = TRUE)) {
  install.packages("devtools", repos = "https://cran.r-project.org")
}

message("✓ devtools ready")

message("\nStep 2: Installing glycofragR from local source...")
message("This may take a moment...")

# Get the current working directory (should be the package root)
pkg_root <- getwd()
glycofragr_path <- file.path(pkg_root, "glycofragR")

message(paste("Package path:", glycofragr_path))

# Install the package
devtools::install(pkg = glycofragr_path, force = TRUE, upgrade = "never")

message("\n✓ glycofragR installed successfully!")
message("\nYou can now use: library(glycofragR)")
message("\nNext steps:")
message("1. Run: source('test_basic.R')")
message("2. Or: source('GlycofragR_test_UPDATED.R')")
