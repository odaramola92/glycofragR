# GlypPRM R Wrapper - Python DataFrame to R-Compatible Format Conversion

## Summary

You can now receive R-compatible data formats (native data.frames) from the Python analysis classes instead of Python pandas DataFrames. This makes working with the data in R much simpler and more intuitive.

## What Was Changed

### 1. **Python Code Changes** (2 files)

#### `glycofrag/glycan_analysis.py` - GlycanAnalysis class
Added 4 new methods:
- `get_theoretical_df_as_list()` → Returns list of dictionaries
- `get_clean_df_as_list()` → Returns list of dictionaries  
- `get_summary_df_as_list()` → Returns list of dictionaries
- `get_structures_as_list()` → Returns list of dictionaries

These methods convert pandas DataFrames into Python list-of-dictionaries format, which reticulate (Python/R bridge) handles much better.

#### `glycofrag/analysis.py` - GlycoPeptideAnalysis class
Added the same 4 methods (identical implementation) for glycopeptide analysis.

---

### 2. **R Code Changes** (2 files)

#### `glycofragR/R/glycan_analysis.R` - Wrapper functions
Updated functions:
- `glycan_analysis_theoretical_df()` - Now calls `get_theoretical_df_as_list()` and converts to data.frame
- `glycan_analysis_clean_df()` - Now calls `get_clean_df_as_list()` and converts to data.frame
- `glycan_analysis_summary()` - Now calls `get_summary_df_as_list()` and converts to data.frame
- `glycan_analysis_structures()` - Now calls `get_structures_as_list()` and converts to data.frame

#### `glycofragR/R/utils.R` - Helper function
Added:
- `.list_dict_to_df()` - Internal helper that converts Python list-of-dicts to R data.frame

---

## Technical Details

### Why This Works Better

**Problem:**
```
Python pandas DataFrame 
    ↓ (reticulate conversion)
R PyObject (doesn't behave like a data.frame)
    ↓ (user has to manually convert)
py_to_r() 
    ↓
R data.frame (finally!)
```

**Solution:**
```
Python list of dictionaries
    ↓ (reticulate handles this well)
R list (easily converted)
    ↓ (our helper function)
R data.frame (native, immediately usable!)
```

### Conversion Process

1. Python methods convert pandas DataFrames to list-of-dicts format
   ```python
   df.fillna("").to_dict('records')
   # [{col1: val1, col2: val2}, {col1: val3, col2: val4}, ...]
   ```

2. R receives the list and converts to data.frame
   ```r
   .list_dict_to_df(list_of_dicts)
   # Returns native R data.frame
   ```

---

## Usage Examples

### Simple Export (Most Common Use Case)
```r
library(glycofragR)
library(writexl)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Get data - now returns native R data.frames!
theoretical <- glycan_analysis_theoretical_df(analysis)
clean <- glycan_analysis_clean_df(analysis)
summary <- glycan_analysis_summary(analysis)

# Export directly - no conversion needed
write_xlsx(
  list(
    Theoretical = theoretical,
    Clean = clean,
    Summary = summary
  ),
  path = "results.xlsx"
)
```

### Data Manipulation with dplyr
```r
library(dplyr)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")
clean_df <- glycan_analysis_clean_df(analysis)

# Use native R operations
filtered <- clean_df %>%
  filter(mass > 500, mass < 2000) %>%
  arrange(mass) %>%
  select(mass, composition, m_z)
```

### Statistical Analysis
```r
summary_df <- glycan_analysis_summary(analysis)

# R's summary functions work directly
summary(summary_df)
str(summary_df)
```

---

## File Locations

**Files Modified:**
1. `glycofrag/glycan_analysis.py` - Added 4 conversion methods
2. `glycofrag/analysis.py` - Added 4 conversion methods
3. `glycofragR/R/glycan_analysis.R` - Updated 4 wrapper functions
4. `glycofragR/R/utils.R` - Added helper function

**New Documentation:**
1. `R_WRAPPER_USAGE_UPDATE.md` - Complete usage guide with before/after examples
2. `GlycofragR_test_UPDATED.R` - Example R script showing new usage

---

## Backward Compatibility

⚠️ **Breaking Change**: The old approach using `py_to_r()` on the results is no longer necessary.

**Before:**
```r
theoretical <- glycan_analysis_theoretical_df(analysis)
theoretical_df <- py_to_r(theoretical)  # ← Old line no longer needed
```

**After:**
```r
theoretical_df <- glycan_analysis_theoretical_df(analysis)  # ← Direct native data.frame
```

If you still use the old `py_to_r()` approach, it will still work (R will just convert a data.frame to a data.frame), but it's unnecessary.

---

## Testing the Changes

Run the updated test script to verify everything works:
```r
source("GlycofragR_test_UPDATED.R")
```

This should:
- ✓ Create analysis
- ✓ Get data.frames natively (no conversion needed)
- ✓ Export to Excel
- ✓ Manipulate data with dplyr
- ✓ Calculate statistics

---

## Benefits Gained

| Aspect | Before | After |
|--------|--------|-------|
| **Data Type** | Python DataFrame | R data.frame |
| **R Compatibility** | Requires conversion | Native, immediate |
| **Required Packages** | reticulate | None additional |
| **Extra Steps** | py_to_r() needed | Data ready to use |
| **tidyverse Support** | Limited | Full support |
| **Learning Curve** | Steep | Gentle |

---

## Future Enhancement Ideas

If you want even more R integration in the future, consider:

1. **Create S3 classes** - Wrap analysis results in custom R classes
2. **Add S3 methods** - Custom `print()`, `plot()`, `summary()` for results
3. **Create tibble versions** - Return tibbles instead of data.frames for modern R
4. **Add more utility functions** - Mass conversion helpers, composition parsing, etc.

---

## Questions or Issues?

If you encounter any issues with the new format:
1. Check the documentation in `R_WRAPPER_USAGE_UPDATE.md`
2. Review the example in `GlycofragR_test_UPDATED.R`
3. Ensure Python code is using the new methods (not old `theoretical_df` attribute directly)
4. Verify R package is installed with latest changes

---

**Summary**: Your R wrapper now returns universal, R-compatible data.frames directly, eliminating the need for manual Python-to-R conversion. Much simpler! 🎉
