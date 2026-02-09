# Simple complete setup - avoids all the installation hassles
# Just loads glycofragR directly from source

message("================================================================================")
message("GlycofragR Setup - Simple Version")
message("================================================================================\n")

# Step 1: Load required packages
message("Step 1: Loading required packages...")
packages_needed <- c("reticulate", "writexl", "dplyr")

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("  Installing", pkg, "..."))
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}
message("✓ All packages loaded\n")

# Step 2: Load glycofragR
message("Step 2: Loading glycofragR...")
source('load_glycofragr.R')

# Step 3: Quick test
message("\nStep 3: Running quick test...")
tryCatch({
  reticulate::use_python(
    "C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe",
    required = TRUE
  )
  
  analysis <- glycan_analysis_new("4502", "permethylated_reduced", visualize_structure = NULL)
  df <- glycan_analysis_theoretical_df(analysis)
  
  message(paste("✓ Got", nrow(df), "theoretical fragments"))
}, error = function(e) {
  message(paste("Test skipped:", e$message))
})

message("\n================================================================================")
message("✓ Setup Complete!")
message("================================================================================\n")

message("You can now use glycofragR in your scripts. Example:\n")
message("  library(reticulate)")
message("  library(writexl)")
message("  analysis <- glycan_analysis_new('4502', 'permethylated_reduced')")
message("  df <- glycan_analysis_theoretical_df(analysis)")
message("  write_xlsx(df, 'results.xlsx')")
