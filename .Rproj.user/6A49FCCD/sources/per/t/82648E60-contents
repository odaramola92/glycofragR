# Glycan Structure Prediction and Fragmentation Example (R version)
# Equivalent to GLycan_fragments_example.py
# Compare R and Python results to verify R wrapper works correctly

library(glycofragR)
library(reticulate)
library(writexl)

# Configure Python environment
use_python("C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe",
           required = TRUE)

# Force non-GUI backend so PNGs can be saved without Tcl/Tk
reticulate::py_run_string("import matplotlib; matplotlib.use('Agg')")

# ========== CONFIGURATION ==========

GLYCAN_TYPE <- 'N'                      # 'N' for N-glycans, 'O' for O-glycans
GLYCAN_CODE <- '6632'                   # 5-digit composition code
MODIFICATION_TYPE <- 'permethylated_reduced'  # Reducing end modification
ISOMER_SENSITIVE <- FALSE               # Treat mirror images as distinct
MAX_STRUCTURES <- 100                   # Maximum structures to predict
FRAGMENT_TYPES <- c('BY')               # c('BY'), c('CZ'), or c('BY', 'CZ')
CHARGES <- c(1, 2, 3)                   # Charge states
VISUALIZE_STRUCTURE <- 'all'            # 'all', 'best', int, list, or NULL

OUTPUT_FILE <- paste0("glycan_analysis_", GLYCAN_CODE, "_", GLYCAN_TYPE, ".xlsx")

# =================================================

cat("================================================================================\n")
cat("GLYCAN FRAGMENTATION ANALYSIS (R VERSION)\n")
cat("================================================================================\n")
cat(sprintf("Glycan Type: %s, Code: %s\n", GLYCAN_TYPE, GLYCAN_CODE))
cat(sprintf("Modification: %s\n", MODIFICATION_TYPE))
cat(sprintf("Configuration: isomer_sensitive=%s, fragment_types=%s\n", 
            ISOMER_SENSITIVE, paste(FRAGMENT_TYPES, collapse=", ")))
cat("================================================================================\n\n")

# Create GlycanAnalysis
tryCatch({
  analysis <- glycan_analysis_new(
    glycan_code = GLYCAN_CODE,
    glycan_type = GLYCAN_TYPE,
    modification_type = MODIFICATION_TYPE,
    max_structures = MAX_STRUCTURES,
    isomer_sensitive = ISOMER_SENSITIVE,
    fragment_types = FRAGMENT_TYPES,
    charges = CHARGES,
    visualize_structure = VISUALIZE_STRUCTURE
  )
  cat("[OK] GlycanAnalysis created successfully\n")
}, error = function(e) {
  cat(sprintf("[ERROR] Failed to create analysis: %s\n", e$message))
  stop(e)
})

# Get data as native R data.frames
theoretical_df <- glycan_analysis_theoretical_df(analysis)
clean_df <- glycan_analysis_clean_df(analysis)
summary_df <- glycan_analysis_summary(analysis)
structures_df <- glycan_analysis_structures(analysis)

cat(sprintf("     %d structures predicted\n", nrow(structures_df)))
cat(sprintf("     %d theoretical fragments generated\n", nrow(theoretical_df)))
cat(sprintf("     %d unique fragments (deduplicated)\n\n", nrow(clean_df)))

# Display summary
cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
print(summary_df)

# Display structures
cat("\n================================================================================\n")
cat("PREDICTED STRUCTURES\n")
cat("================================================================================\n")
if (nrow(structures_df) > 0) {
  for (i in seq_len(nrow(structures_df))) {
    cat(sprintf("\nStructure %d:\n", structures_df$structure_id[i]))
    cat(sprintf("  %s\n", structures_df$structure_str[i]))
  }
}

# Display fragment summary
cat("\n================================================================================\n")
cat("FRAGMENT SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Total theoretical fragments: %d\n", nrow(theoretical_df)))
cat(sprintf("Unique fragments (deduplicated): %d\n", nrow(clean_df)))

if (nrow(clean_df) > 0) {
  if ("Fragment Type" %in% colnames(clean_df)) {
    cat("\nFragments by type:\n")
    frag_counts <- table(clean_df$`Fragment Type`)
    for (frag_type in names(frag_counts)) {
      cat(sprintf("  %s: %d\n", frag_type, frag_counts[frag_type]))
    }
  }
}

# Show sample fragments
cat("\n================================================================================\n")
cat("SAMPLE FRAGMENTS (First 10 rows)\n")
cat("================================================================================\n")
if (nrow(clean_df) > 0) {
  print(head(clean_df, 10))
}

# Export to Excel (no structures sheet; PNGs are saved via VISUALIZE_STRUCTURE)
cat("\n================================================================================\n")
cat("EXPORTING TO EXCEL\n")
cat("================================================================================\n")

tryCatch({
  write_xlsx(
    list(
      Summary = summary_df,
      Theoretical_Fragments = theoretical_df,
      Clean_Fragments = clean_df
    ),
    path = OUTPUT_FILE
  )
  cat(sprintf("[OK] Exported to: %s\n", OUTPUT_FILE))
}, error = function(e) {
  cat(sprintf("[WARNING] Could not export to Excel: %s\n", e$message))
})

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")

cat("\nComparison with Python:\n")
cat("- Check if structure count matches Python version\n")
cat("- Check if fragment counts match Python version\n")
cat("- Verify data.frame values are identical\n")
