# Migration Guide: From py_to_r() to Native R data.frames

## Overview

Your R wrapper has been updated to return native R `data.frame` objects directly from all analysis functions. This eliminates the need for manual `py_to_r()` conversion.

## Step-by-Step Migration

### Step 1: Identify Old Code

Look for patterns like this in your R scripts:

```r
# OLD PATTERN 1: Using py_to_r() on results
theoretical <- glycan_analysis_theoretical_df(analysis)
theoretical_df <- py_to_r(theoretical)

# OLD PATTERN 2: Chaining py_to_r()
clean_df <- py_to_r(glycan_analysis_clean_df(analysis))

# OLD PATTERN 3: Multiple conversions
theoretical <- py_to_r(glycan_analysis_theoretical_df(analysis))
clean <- py_to_r(glycan_analysis_clean_df(analysis))
summary <- py_to_r(glycan_analysis_summary(analysis))
```

### Step 2: Replace With New Code

Simply remove the `py_to_r()` wrapper:

```r
# NEW PATTERN 1: Direct use
theoretical_df <- glycan_analysis_theoretical_df(analysis)

# NEW PATTERN 2: Direct assignment
clean_df <- glycan_analysis_clean_df(analysis)

# NEW PATTERN 3: All direct
theoretical_df <- glycan_analysis_theoretical_df(analysis)
clean_df <- glycan_analysis_clean_df(analysis)
summary_df <- glycan_analysis_summary(analysis)
```

---

## Before and After Examples

### Example 1: Basic Workflow

**Before (Old Way):**
```r
library(glycofragR)
library(openxlsx)
library(reticulate)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Step 1: Get results
theoretical <- glycan_analysis_theoretical_df(analysis)
clean <- glycan_analysis_clean_df(analysis)
summary <- glycan_analysis_summary(analysis)

# Step 2: Convert to R data.frames
theoretical_df <- py_to_r(theoretical)
clean_df <- py_to_r(clean)
summary_df <- py_to_r(summary)

# Step 3: Export
write_xlsx(
  list(theoretical = theoretical_df, clean = clean_df, summary = summary_df),
  path = "results.xlsx"
)
```

**After (New Way):**
```r
library(glycofragR)
library(writexl)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Direct assignment - no conversion step!
theoretical_df <- glycan_analysis_theoretical_df(analysis)
clean_df <- glycan_analysis_clean_df(analysis)
summary_df <- glycan_analysis_summary(analysis)

# Export directly
write_xlsx(
  list(theoretical = theoretical_df, clean = clean_df, summary = summary_df),
  path = "results.xlsx"
)
```

---

### Example 2: Data Analysis

**Before (Old Way):**
```r
library(dplyr)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Get and convert
clean <- glycan_analysis_clean_df(analysis)
clean_df <- py_to_r(clean)

# Now do analysis
filtered <- clean_df %>%
  filter(mass > 500) %>%
  arrange(mass)
```

**After (New Way):**
```r
library(dplyr)

analysis <- glycan_analysis_new("4502", "permethylated_reduced")

# Direct analysis - no conversion!
filtered <- glycan_analysis_clean_df(analysis) %>%
  filter(mass > 500) %>%
  arrange(mass)
```

---

### Example 3: Batch Processing

**Before (Old Way):**
```r
glycan_codes <- c("4501", "4502", "4503")

for (code in glycan_codes) {
  analysis <- glycan_analysis_new(code, "permethylated_reduced")
  
  # Each one needs conversion
  clean <- glycan_analysis_clean_df(analysis)
  clean_df <- py_to_r(clean)
  
  summary <- glycan_analysis_summary(analysis)
  summary_df <- py_to_r(summary)
  
  # Write results
  write.csv(clean_df, paste0("clean_", code, ".csv"))
  write.csv(summary_df, paste0("summary_", code, ".csv"))
}
```

**After (New Way):**
```r
glycan_codes <- c("4501", "4502", "4503")

for (code in glycan_codes) {
  analysis <- glycan_analysis_new(code, "permethylated_reduced")
  
  # No conversion needed!
  clean_df <- glycan_analysis_clean_df(analysis)
  summary_df <- glycan_analysis_summary(analysis)
  
  # Write directly
  write.csv(clean_df, paste0("clean_", code, ".csv"))
  write.csv(summary_df, paste0("summary_", code, ".csv"))
}
```

---

## Function Reference

### Old Functions → New Equivalents

| Old Code | New Code | Return Type | Notes |
|----------|----------|-------------|-------|
| `glycan_analysis_theoretical_df(analysis)` + `py_to_r()` | `glycan_analysis_theoretical_df(analysis)` | data.frame | No conversion needed |
| `glycan_analysis_clean_df(analysis)` + `py_to_r()` | `glycan_analysis_clean_df(analysis)` | data.frame | No conversion needed |
| `glycan_analysis_summary(analysis)` + `py_to_r()` | `glycan_analysis_summary(analysis)` | data.frame | No conversion needed |
| `glycan_analysis_structures(analysis)` + `py_to_r()` | `glycan_analysis_structures(analysis)` | data.frame | Now returns data.frame |

---

## Common Issues & Solutions

### Issue: "No such variable" error

**Problem:**
```r
# You removed the py_to_r() line but still use the old variable name
theoretical <- glycan_analysis_theoretical_df(analysis)  # Forgot to rename or convert
head(theoretical)  # Error: theoretical is a data.frame, not a list-of-dicts
```

**Solution:**
```r
# Rename or assign correctly
theoretical_df <- glycan_analysis_theoretical_df(analysis)
head(theoretical_df)  # Works!
```

---

### Issue: "Object is not a data.frame" error

**Problem:**
```r
# Still trying to convert with py_to_r()
clean_df <- py_to_r(glycan_analysis_clean_df(analysis))
# py_to_r() attempts to convert a data.frame (no-op), might cause issues
```

**Solution:**
```r
# Remove the py_to_r() wrapper
clean_df <- glycan_analysis_clean_df(analysis)
```

---

### Issue: Column names not recognized

**Problem:**
```r
# After migration, column names might appear different
clean_df <- glycan_analysis_clean_df(analysis)
clean_df$fragment_id  # Returns NULL or error
```

**Solution:**
```r
# Check actual column names
colnames(clean_df)
head(clean_df)

# Use the correct column name
clean_df$m_z  # Or whatever the actual column is
```

---

## Backward Compatibility

The new functions are **compatible** with the old approach:

```r
# This still works (data.frame → data.frame is a no-op)
theoretical_df <- py_to_r(glycan_analysis_theoretical_df(analysis))
# But it's unnecessary - just use:
theoretical_df <- glycan_analysis_theoretical_df(analysis)
```

So you can migrate gradually - there's no breaking change, just unnecessary steps.

---

## Testing Your Migration

Use the provided test script to verify your setup:

```r
# Load the updated test script
source("GlycofragR_test_UPDATED.R")

# Should run without errors and produce:
# ✓ Example 1 complete: Data exported to Excel
# ✓ Example 2 complete: Filtered and sorted fragments
# ✓ All examples completed successfully!
```

---

## Optional Cleanup

If you had helper scripts that only handled `py_to_r()` conversion, you can now:

1. **Remove unused imports:**
   ```r
   # No longer needed explicitly
   # library(reticulate)  # Still needed for use_python(), but not for conversion
   ```

2. **Delete conversion helper functions** (if you had any)
   ```r
   # Remove any custom py_to_r wrappers you created
   ```

3. **Simplify existing functions** that were handling conversion

---

## Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Extra conversion step** | Required | Not needed |
| **Code complexity** | Higher | Lower |
| **Performance** | Slightly slower (extra step) | Faster (direct) |
| **R compatibility** | Limited | Full |
| **Lines of code** | More | Less |

---

## Questions?

Refer to these files for more info:
- `QUICK_REFERENCE.md` - 30-second summary
- `R_WRAPPER_USAGE_UPDATE.md` - Complete usage guide
- `PYTHON_TO_R_CONVERSION_SUMMARY.md` - Technical details
- `GlycofragR_test_UPDATED.R` - Working examples

Happy coding! 🎉
