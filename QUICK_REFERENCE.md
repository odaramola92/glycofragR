# Quick Reference - R-Friendly Data Format

## The Problem (Now Solved)
Python pandas DataFrames were being returned to R, but they're not easy to work with in R. Now you get native R data.frames!

## What Changed - 30 Second Summary

### Before ❌
```r
analysis <- glycan_analysis_new("4502", "permethylated_reduced")
theoretical <- glycan_analysis_theoretical_df(analysis)
# Returns: Python pandas DataFrame object (not R-friendly)
# Required: py_to_r(theoretical) to convert
```

### After ✅
```r
analysis <- glycan_analysis_new("4502", "permethylated_reduced")
theoretical_df <- glycan_analysis_theoretical_df(analysis)
# Returns: Native R data.frame (ready to use immediately!)
```

## Methods Available (All Return Native R data.frames)

```r
# Get analysis results as native R data.frames
theoretical_df <- glycan_analysis_theoretical_df(analysis)
clean_df <- glycan_analysis_clean_df(analysis)
summary_df <- glycan_analysis_summary(analysis)
structures_df <- glycan_analysis_structures(analysis)

# All are native R data.frames - use them directly!
nrow(theoretical_df)           # Works!
head(clean_df)                 # Works!
write.csv(summary_df, "out.csv") # Works!
```

## Common Use Cases

### 1. Export to Excel
```r
write_xlsx(
  list(Theoretical = theoretical_df, Clean = clean_df, Summary = summary_df),
  path = "results.xlsx"
)
```

### 2. Filter & Sort
```r
significant <- clean_df %>%
  filter(mass > 500 & mass < 2000) %>%
  arrange(mass)
```

### 3. Get Statistics
```r
summary(clean_df$mass)
cat("Total fragments:", nrow(clean_df), "\n")
```

### 4. View Data
```r
head(theoretical_df)
str(clean_df)
View(summary_df)
```

## No More Needed
```r
# DON'T do this anymore:
theoretical <- glycan_analysis_theoretical_df(analysis)
theoretical_df <- py_to_r(theoretical)  # ← No longer necessary!

# JUST do this:
theoretical_df <- glycan_analysis_theoretical_df(analysis)  # ← Ready to use!
```

## What the Code Does (Behind the Scenes)

Python side:
```python
# New method in GlycanAnalysis class
def get_theoretical_df_as_list(self):
    return self.theoretical_df.fillna("").to_dict('records')
    # Converts pandas → list of dictionaries
```

R side:
```r
# New helper function in utils.R
.list_dict_to_df <- function(list_of_dicts) {
  do.call(rbind, lapply(list_of_dicts, as.data.frame))
  # Converts list of dicts → R data.frame
}
```

Wrapper function:
```r
glycan_analysis_theoretical_df <- function(analysis) {
  frag_list <- analysis$get_theoretical_df_as_list()  # Get list from Python
  .list_dict_to_df(frag_list)  # Convert to R data.frame
}
```

## Files Changed

**Python:**
- `glycofrag/glycan_analysis.py` - Added 4 methods
- `glycofrag/analysis.py` - Added 4 methods

**R:**
- `glycofragR/R/glycan_analysis.R` - Updated 4 functions
- `glycofragR/R/utils.R` - Added 1 helper function

**Documentation:**
- `R_WRAPPER_USAGE_UPDATE.md` - Complete guide
- `GlycofragR_test_UPDATED.R` - Working example
- `PYTHON_TO_R_CONVERSION_SUMMARY.md` - Technical details

## Status

✅ Python code: Ready (new methods added)
✅ R code: Ready (updated all wrappers)  
✅ Documentation: Complete
✅ Testing: See `GlycofragR_test_UPDATED.R`

Just use it! No more conversion steps needed.
