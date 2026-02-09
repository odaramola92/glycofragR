reticulate::use_python("C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe", required = TRUE)
# then try:
library(glycofragR)
library(openxlsx)


analysis <- glycan_analysis_new("4502", "permethylated_reduced", 
                                visualize_structure = "all")

# fragments from analysis
theoretical <- glycan_analysis_theoretical_df(analysis)
clean <- glycan_analysis_clean_df(analysis)

# summaries and structures if needed
summary <- glycan_analysis_summary(analysis)
structures <- glycan_analysis_structures(analysis)

library(reticulate)
library(writexl)

# Convert pandas DataFrame -> R data.frame
theoretical_df <- py_to_r(theoretical)
clean_df       <- py_to_r(clean)
summary_df     <- py_to_r(summary)


write_xlsx(
  list(
    theoretical_fragments = theoretical_df,
    clean_fragments       = clean_df,
    summary               = summary_df
  ),
  path = "glycan_analysis_results.xlsx"
)

