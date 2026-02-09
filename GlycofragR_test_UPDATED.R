reticulate::use_python("C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe", required = TRUE)
library(glycofragR)
library(writexl)
library(dplyr)

# ===== EXAMPLE 1: Simple Data Export =====
# This now works smoothly without py_to_r() conversion!

analysis <- glycan_analysis_new("4502", "permethylated_reduced", 
                                visualize_structure = NULL)  # Set to NULL to skip visualization

# Get fragments as native R data.frames
theoretical_df <- glycan_analysis_theoretical_df(analysis)
clean_df <- glycan_analysis_clean_df(analysis)
summary_df <- glycan_analysis_summary(analysis)

# Directly export to Excel - no more conversion needed!
write_xlsx(
  list(
    theoretical_fragments = theoretical_df,
    clean_fragments = clean_df,
    summary = summary_df
  ),
  path = "glycan_analysis_results.xlsx"
)

cat("✓ Example 1 complete: Data exported to Excel\n")

# ===== EXAMPLE 2: Data Manipulation with dplyr =====

# Filter based on mass range
filtered <- clean_df %>%
  filter(mass > 500 & mass < 2000) %>%
  arrange(mass) %>%
  head(10)

cat("✓ Example 2 complete: Filtered and sorted fragments\n")
print(filtered)

# ===== EXAMPLE 3: Summary Statistics =====

print("Summary Information:")
print(summary_df)

cat("\nTheoretical fragments count:", nrow(theoretical_df), "\n")
cat("Clean fragments count:", nrow(clean_df), "\n")

# ===== EXAMPLE 4: Access Structures =====

structures_df <- glycan_analysis_structures(analysis)
cat("\nStructures predicted:", nrow(structures_df), "\n")
print(structures_df)

# ===== EXAMPLE 5: Custom Analysis =====

if (nrow(clean_df) > 0) {
  # Get unique compositions
  unique_compositions <- clean_df %>%
    distinct(composition) %>%
    nrow()
  
  cat("\nUnique compositions:", unique_compositions, "\n")
  
  # Get mass range
  mass_range <- clean_df %>%
    summarise(
      min_mass = min(mass, na.rm = TRUE),
      max_mass = max(mass, na.rm = TRUE),
      mean_mass = mean(mass, na.rm = TRUE),
      n_fragments = n()
    )
  
  print(mass_range)
}

cat("\n✓ All examples completed successfully!\n")
cat("✓ Notice: No py_to_r() needed - data.frames work natively in R!\n")
