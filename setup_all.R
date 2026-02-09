# Complete Setup: Install everything in one script
# Run this once to set up your R environment completely

message("================================================================================")
message("GlycofragR Complete Setup")
message("================================================================================\n")

# Step 1: Install CRAN packages
message("Step 1/4: Installing CRAN packages (reticulate, writexl, dplyr)...")
source('install_packages.R')

message("\n================================================================================\n")

# Step 2: Install glycofragR
message("Step 2/4: Installing glycofragR package...")
source('install_glycofragr.R')

message("\n================================================================================\n")

# Step 3: Run basic test
message("Step 3/4: Running basic functionality test...")
source('test_basic.R')

message("\n================================================================================")
message("Setup Complete!")
message("================================================================================\n")

message("You can now use glycofragR in your R scripts:")
message("  library(glycofragR)")
message("\nNext steps:")
message("  1. Run full test: source('GlycofragR_test_UPDATED.R')")
message("  2. Or start using glycofragR in your own code")
message("\nFor more information, see:")
message("  - QUICK_REFERENCE.md")
message("  - R_WRAPPER_USAGE_UPDATE.md")
message("  - SETUP_INSTRUCTIONS.md")
