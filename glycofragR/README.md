# glycofragR

R wrapper for the glycofrag Python package. Handles Python environment setup automatically.

## Quick Start

**First time only:**
```r
library(glycofragR)
glycofrag_install()
```

**After that:**
```r
library(glycofragR)

# Glycan example
glycan <- glycan("4501", "N")
structures <- glycan_predict_structures(glycan)
frags <- glycan_generate_fragments(glycan, structures[[1]])

# Glycopeptide example
gp <- glycopeptide(
  peptide_sequence = "LCPDCPLLAPLNDSR",
  glycan_code = "4501",
  glycosylation_site = 12,
  glycan_type = "N"
)
frags_gp <- glycopeptide_fragments(gp)
```

## Installation

### Prerequisites
- R 4.0+
- Windows, Mac, or Linux

### Setup
```r
# Install the R package (once)
remotes::install_github("your_org/glycofragR")

# Initialize Python environment (once per system)
library(glycofragR)
glycofrag_install()  # This creates conda env and installs glycofrag
```

## Requirements

- R 4.0+
- Python (automatically installed via Miniconda during `glycofrag_install()`)
- glycofrag (automatically installed during setup)

## How It Works

- `.onLoad()`: Automatically configures the glycofrag conda environment
- `glycofrag_install()`: One-time setup (creates conda env, installs dependencies)
- Error messages: Guide you if something isn't set up

## Advanced Usage

### Use a different conda environment
```r
set_glycofrag_conda_env("my_custom_env")
```

### Reinstall glycofrag (troubleshooting)
```r
glycofrag_install(force_reinstall = TRUE)
```

### Manual conda setup (if needed)
```bash
conda create -n glycofrag python=3.11
conda activate glycofrag
pip install glycofrag
```

Then in R:
```r
library(glycofragR)
set_glycofrag_conda_env("glycofrag")
```

## Notes

- All peptide and glycosylation positions are 1-based
- This package respects the GLYCOFRAG_DISABLE_AUTOENV environment variable if you need to override conda env detection
- The Python environment is isolated to avoid conflicts with other projects

## Troubleshooting

If you get an error about glycofrag not being installed:
```r
glycofrag_install()
```

If your conda installation is broken:
```r
glycofrag_install(force_reinstall = TRUE)
```
