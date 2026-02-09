# Quick Setup & Testing Instructions

## Step 1: Install Required CRAN Packages

In your R Interactive console, run:

```r
source('install_packages.R')
```

This will install:
- `reticulate` - Python/R bridge
- `writexl` - Excel export  
- `dplyr` - Data manipulation (optional)

**Expected output:** Package installation complete!

---

## Step 2: Install glycofragR Package

In your R Interactive console, run:

```r
source('install_glycofragr.R')
```

This will:
1. Install `devtools` (if needed)
2. Build and install the `glycofragR` package from the local source

**Expected output:** 
```
✓ glycofragR installed successfully!
```

---

## Step 3: Test Basic Functionality

In your R Interactive console, run:

```r
source('test_basic.R')
```

This test **does NOT require** dplyr or writexl - just the core glycofragR functionality.

**What it tests:**
- Python environment setup
- Creating GlycanAnalysis
- Getting theoretical fragments (native R data.frame)
- Getting clean fragments (native R data.frame)
- Getting summary (native R data.frame)
- Getting structures (native R data.frame)

**Expected output:**
```
✓ Got 45 theoretical fragments
✓ Got 43 clean fragments
✓ Got summary with 1 rows
✓ Got 1 structures
✓ ALL TESTS PASSED!
```

---

## Step 4: Run Full Test Suite (Optional)

After completing Steps 1-3, run:

```r
source('GlycofragR_test_UPDATED.R')
```

This runs comprehensive tests including:
- Data export to Excel
- Data manipulation with dplyr
- Statistical analysis
- Advanced examples

---

## If You Get Errors

### Error: "package 'reticulate' not found"

**Solution:** Run Step 1 first to install packages.

```r
source('install_packages.R')
```

---

### Error: "there is no package called 'glycofragR'"

**Solution:** Run Step 2 to install the package.

```r
source('install_glycofragr.R')
```

This will build and install glycofragR from the local source code.

---

### Error: "Python executable not found"

The test script uses:
```
C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe
```

Make sure this path is correct. If your Python is in a different location, edit the path in `test_basic.R` (line 7).

---

### Error: "glycofrag module not found"

Make sure the Python glycofrag package is installed in your test_env. Run in PowerShell:

```powershell
& "C:\Users\oluwa\OneDrive - Texas Tech University\GlypPRM\Package\test_env\Scripts\python.exe" -m pip install -e .
```

(Run from the package root where `pyproject.toml` is located)

---

## File Summary

| File | Purpose | Dependencies |
|------|---------|--------------|
| `install_packages.R` | Install CRAN packages (reticulate, writexl, dplyr) | Internet |
| `install_glycofragr.R` | Build & install glycofragR from source | devtools |
| `test_basic.R` | Quick test of core functionality | reticulate, glycofragR |
| `GlycofragR_test_UPDATED.R` | Full comprehensive test | reticulate, writexl, dplyr, glycofragR |

---

## Quick Copy-Paste Commands

### Install everything in order:
```r
source('install_packages.R')    # Step 1
source('install_glycofragr.R')  # Step 2
source('test_basic.R')          # Step 3
```

### Test basic (quick, after setup):
```r
source('test_basic.R')
```

### Test everything (full, after setup):
```r
source('GlycofragR_test_UPDATED.R')
```

---

## Troubleshooting

If you're having issues with your R installation:

1. Make sure R version is 4.0+
2. Check Python paths are correct
3. Verify glycofrag Python package is installed in your python environment
4. Try running tests one step at a time

**Quick diagnostics:**
```r
# Check devtools is available
library(devtools)

# Check reticulate is installed
library(reticulate)

# Check Python is found
py_config()

# Check glycofrag is available in Python
reticulate::py_evaluate("import glycofrag; print(glycofrag.__version__)")
```

---

**Start with the steps above in order.** 🚀

1️⃣ `source('install_packages.R')`  
2️⃣ `source('install_glycofragr.R')`  
3️⃣ `source('test_basic.R')`

