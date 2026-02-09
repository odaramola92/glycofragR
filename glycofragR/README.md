# glycofragR

R wrapper for the glycofrag Python package. This package uses reticulate to
access glycofrag from a conda environment.

## Requirements

- R 4.0+
- Python (via conda)
- glycofrag installed in a conda env

## Setup (CRAN-safe)

Create a conda env and install glycofrag:

```bash
conda create -n glycofrag python=3.10
conda activate glycofrag
pip install glycofrag
```

In R:

```r
library(glycofragR)
set_glycofrag_conda_env("glycofrag")
```

## Example

```r
library(glycofragR)
set_glycofrag_conda_env("glycofrag")

# Glycan example
glycan <- glycan_new("4501", "N")
structures <- glycan_predict_structures(glycan)
frags <- glycan_generate_fragments(glycan, structures[[1]])

# Glycopeptide example
gp <- glycopeptide_new(
  peptide_sequence = "LCPDCPLLAPLNDSR",
  glycan_code = "4501",
  glycosylation_site = 12,
  glycan_type = "N"
)
frags_gp <- glycopeptide_fragments(gp)
```

## Notes

- All peptide and glycosylation positions are 1-based.
- This package does not auto-install Python or conda environments.
