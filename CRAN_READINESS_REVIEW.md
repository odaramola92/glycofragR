# CRAN Readiness Review for glycofragR

Date: February 9, 2026
Package Version: 0.0.0.9000

## ✅ STRENGTHS

### Core Package Structure
- ✅ Proper R package structure (DESCRIPTION, NAMESPACE, R/, man/, tests/)
- ✅ MIT License with proper LICENSE file
- ✅ Good documentation with roxygen2 comments
- ✅ Comprehensive test suite (5 test files)
- ✅ Vignettes included (glycofragR_intro.Rmd)
- ✅ Proper .Rbuildignore for CRAN compliance
- ✅ Global variables declared to avoid R CMD CHECK notes

### Dependencies
- ✅ Only depends on reticulate (single, well-maintained dependency)
- ✅ Python package (glycofrag) is external, not bundled
- ✅ Reasonable suggests: testthat, knitr, rmarkdown

### Code Quality
- ✅ Functions properly exported in NAMESPACE
- ✅ Functions have @export roxygen tags
- ✅ Proper use of importFrom(reticulate, ...) in NAMESPACE
- ✅ Package-level documentation exists

---

## ⚠️ ISSUES TO ADDRESS BEFORE CRAN SUBMISSION

### 1. **LICENSE File is Incomplete** 🔴 HIGH PRIORITY
**Current:** Only has Year and Copyright Holder
**Required:** Full MIT license text
**Solution:** Replace with:
```
YEAR: 2026
COPYRIGHT HOLDER: Oluwatosin Daramola
```
Plus full MIT license text (get from https://opensource.org/licenses/MIT)

### 2. **Missing / Incomplete Function Documentation** 🔴 HIGH PRIORITY
**Issue:** Some exported functions may lack @param @return @examples
**Solution:** Check that ALL exported functions have:
- @title (brief description)
- @param (all parameters documented)
- @return (what it returns)
- @examples (working examples)
- @export tag

**Exported but need review:**
- glycan_analysis_theoretical_df()
- glycan_analysis_clean_df()
- glycan_analysis_summary()
- glycan_analysis_structures()
- (And others - check all 22 exported functions)

### 3. **README Section on Python Setup** 🟡 MEDIUM PRIORITY
**Current:** Mentions conda but doesn't address CRAN users
**Issue:** CRAN users expect auto-installation or clear setup
**Solution:** Add section explaining:
```
## Important: Python Installation Required

This package requires glycofrag Python package. CRAN does not allow bundled 
Python packages, so users must install it themselves:

### Option 1: Using conda (recommended)
```r
reticulate::install_miniconda()
reticulate::conda_create("r-glycofrag", packages = "python=3.10")
reticulate::conda_install("r-glycofrag", "glycofrag")
glycofragR::set_glycofrag_conda_env("r-glycofrag")
```

### Option 2: Using pip with existing Python
```r
reticulate::virtualenv_create("r-glycofrag")
reticulate::virtualenv_install("r-glycofrag", "glycofrag")
glycofragR::set_glycofrag_conda_env("r-glycofrag")
```
```

### 4. **DESCRIPTION Fields** 🟡 MEDIUM PRIORITY
**Current:**
```
Version: 0.0.0.9000
```
**Issue:** CRAN discourages .9000 in initial submissions. For CRAN, use 0.1.0
**Also:** Description mentions reticulate but license doesn't mention dependency clearly

### 5. **Test Coverage** 🟡 MEDIUM PRIORITY
**Issue:** Need to verify tests:
- Run without requiring Python installation? (They should?)
- Use mocking/fixtures if tests require glycofrag
- All tests pass without errors

**Action:** Run `devtools::test()` and ensure no failures

### 6. **R CMD CHECK** 🔴 HIGH PRIORITY
**Must pass without errors or warnings:**
```
cd glycofragR
R CMD check .
```
Look for:
- Undeclared global variables ✅ (already handled)
- Missing documentation
- Missing imports
- Vignette build errors

### 7. **URL and Repository Info** 🟡 MEDIUM PRIORITY
**Missing from DESCRIPTION:**
```
URL: https://github.com/yourusername/glycofragR
BugReports: https://github.com/yourusername/glycofragR/issues
```
**Solution:** Add these to DESCRIPTION file

### 8. **Author/Maintainer Email** 🟡 MEDIUM PRIORITY
**Current:** oluwatosindaramola@gmail.com
**CRAN prefers:** A professional email or institutional email
**Consider:** Using university email if available

---

## 📋 CRAN TESTING OPTIONS (Before Final Submission)

### Option 1: **R-Hub** (Recommended - Official CRAN testing)
```r
rhub::check_for_cran()  # Comprehensive CRAN check
```
This tests on:
- Windows (Release, Devel, oldrelease)
- Linux (multiple OS versions)
- macOS

### Option 2: **Win-Builder** (Windows-specific)
```r
devtools::check_win_devel()
devtools::check_win_release()
```
Tests your package on Windows with R development/release versions

### Option 3: **Local CRAN-like Check**
```r
devtools::check(args = c("--as-cran"))
```
Runs the strictest checks locally (not 100% CRAN equivalent but close)

---

## ⚡ RECOMMENDED CRAN SUBMISSION WORKFLOW

1. **Fix Issues (TODAY)**
   - [ ] Replace LICENSE file with complete MIT text
   - [ ] Verify ALL function documentation is complete
   - [ ] Update DESCRIPTION with URL, BugReports, and Version
   - [ ] Add Python setup section to README
   - [ ] Run `R CMD check` locally - fix any warnings

2. **Test on R-Hub (1-2 days)**
   - [ ] Run `rhub::check_for_cran()`
   - [ ] Fix any platform-specific issues
   - [ ] Verify tests pass on Windows, Linux, macOS

3. **Test on Win-Builder (optional but recommended)**
   - [ ] Run `devtools::check_win_devel()`
   - [ ] Verify Windows compatibility

4. **Submit to CRAN**
   - [ ] Create GitHub release with version tag
   - [ ] Submit to CRAN with submission letter
   - [ ] Be prepared for initial review feedback (usually fixed quickly)

---

## 📊 CRAN TIMELINE EXPECTATIONS

- **Initial submission → First feedback:** 1-3 weeks
- **Revision turnaround:** 2-5 days per revision
- **Total time to published:** 2-4 weeks typical
- **Common issues:** Missing documentation, vignette issues, Windows compatibility

---

## 🚀 GITHUB READINESS CHECKLIST

Before pushing to GitHub:

- [ ] Add .gitignore (Python/R standard)
- [ ] Remove DEVELOPMENT PLAN and changelog if internal
- [ ] Add GitHub-specific documentation (CONTRIBUTING.md)
- [ ] Create GitHub Actions CI (optional but good practice)
- [ ] Initial GitHub release with version 0.0.0.9000

---

## FINAL NOTES

**This package is ~85% CRAN-ready. Main blockers:**
1. Complete LICENSE file (2 minutes)
2. Document all exported functions (1-2 hours)
3. Run R CMD CHECK locally (10 minutes)
4. Update DESCRIPTION with metadata (5 minutes)
5. Test with R-Hub (1-2 days)

**Estimated time to CRAN submission: 3-5 days**

Once GitHub is public and you've passed R-Hub checks, CRAN submission is straightforward.
