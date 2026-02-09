# ======= EASIEST WAY TO GET STARTED =======

# Option 1: Just load glycofragR (fastest)
source('load_glycofragr.R')

# Option 2: Full setup with test
source('setup_simple.R')

# Option 3: Just run the test
source('test_simple.R')

# ========================================
# Then use it:

library(reticulate)
library(writexl)

# Create analysis
analysis <- glycan_analysis_new('4502', 'permethylated_reduced')

# Get data (native R data.frames, no conversion needed!)
theoretical <- glycan_analysis_theoretical_df(analysis)
clean <- glycan_analysis_clean_df(analysis)
summary <- glycan_analysis_summary(analysis)

# Use native R operations
head(clean)
nrow(clean)
colnames(clean)

# Export
write_xlsx(list(clean = clean, summary = summary), path = 'results.xlsx')

# Filter with dplyr
library(dplyr)
filtered <- clean %>%
  filter(mass > 500, mass < 2000) %>%
  arrange(mass)

# ========================================
# Files available:

# Quick load:
# source('load_glycofragr.R')

# Simple setup:
# source('setup_simple.R')

# Simple test:
# source('test_simple.R')

# Quickstart examples:
# source('quickstart.R')

# Full test (requires all packages):
# source('GlycofragR_test_UPDATED.R')
