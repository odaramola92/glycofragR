# R Wrapper - Updated Usage with Native Data.Frame Support

## Problem Solved
Previously, the R wrapper was returning Python pandas DataFrames, which didn't work smoothly in R. Now all DataFrames are automatically converted to native R `data.frame` objects!

## Changes Made

### Python Side (glycofrag)
Added new methods to both `GlycanAnalysis` and `GlycoPeptideAnalysis` classes:
- `get_theoretical_df_as_list()` - Returns theoretical fragments as list of dictionaries
- `get_clean_df_as_list()` - Returns clean/deduplicated fragments as list of dictionaries  
- `get_summary_df_as_list()` - Returns summary info as list of dictionaries
- `get_structures_as_list()` - Returns structure information as list of dictionaries

### R Side (glycofragR)
Updated wrapper functions to use the new conversion methods:
- `glycan_analysis_theoretical_df()` - Now returns native R data.frame
- `glycan_analysis_clean_df()` - Now returns native R data.frame
- `glycan_analysis_summary()` - Now returns native R data.frame
- `glycan_analysis_structures()` - Now returns native R data.frame

Added helper function:
- `.list_dict_to_df()` - Utility to convert Python list-of-dicts to R data.frame

## Usage - Before vs After

### BEFORE (Old Way - NO LONGER NEEDED)
```r
reticulate::use_python("C:/path/to/python.exe", required = TRUE)
library(glycofragR)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Had to use py_to_r() to convert Python dataframes
theoretical <- glycan_analysis_theoretical_df(analysis)
theoretical_df <- py_to_r(theoretical)  # ← Conversion step needed

clean <- glycan_analysis_clean_df(analysis)
clean_df <- py_to_r(clean)  # ← Conversion step needed

summary <- glycan_analysis_summary(analysis)
summary_df <- py_to_r(summary)  # ← Conversion step needed
```

### AFTER (New Way - MUCH SIMPLER!)
```r
reticulate::use_python("C:/path/to/python.exe", required = TRUE)
library(glycofragR)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Functions now return native R data.frames directly!
theoretical_df <- glycan_analysis_theoretical_df(analysis)  # Already a data.frame!
clean_df <- glycan_analysis_clean_df(analysis)              # Already a data.frame!
summary_df <- glycan_analysis_summary(analysis)             # Already a data.frame!
structures_df <- glycan_analysis_structures(analysis)       # Already a data.frame!

# Use them directly with normal R operations
head(theoretical_df)
nrow(clean_df)
write.csv(summary_df, "summary.csv")
```

## Example: Export to Excel

```r
library(writexl)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Now you can use the data.frames directly!
write_xlsx(
  list(
    theoretical_fragments = glycan_analysis_theoretical_df(analysis),
    clean_fragments = glycan_analysis_clean_df(analysis),
    summary = glycan_analysis_summary(analysis)
  ),
  path = "glycan_analysis_results.xlsx"
)
```

## Example: Work with Data in R

```r
library(dplyr)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Get data as native R data.frames
clean_frags <- glycan_analysis_clean_df(analysis)

# Use dplyr for data manipulation
filtered <- clean_frags %>%
  filter(mass > 500 & mass < 2000) %>%
  arrange(mass) %>%
  select(mass, composition, m_z)

# View summary statistics
summary_data <- glycan_analysis_summary(analysis)
print(summary_data)

# Plot if you have visualization needs
library(ggplot2)
ggplot(clean_frags, aes(x = mass, y = m_z)) +
  geom_point() +
  theme_minimal()
```

## Benefits

✅ **No more conversion headaches** - Native R data.frames from the start  
✅ **Full R compatibility** - Use with tidyverse, data.table, base R functions  
✅ **Simpler code** - No `py_to_r()` calls needed  
✅ **Better performance** - Direct conversion avoids intermediate steps  
✅ **Consistent interface** - Works the same for both `GlycanAnalysis` and `GlycoPeptideAnalysis`  

## Notes

- Empty results return empty `data.frame()` objects (not NULL)
- Missing/NaN values are converted to empty strings for R compatibility
- All methods work the same for both glycan-only and glycopeptide analysis
- For structure objects, the data.frame includes structure_id, structure_str, and structure_repr columns
