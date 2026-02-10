# glycofragR - Glycan and Glycopeptide Fragmentation Analysis in R

[![CRAN status](https://www.r-pkg.org/badges/version/glycofragR)](https://CRAN.R-project.org/package=glycofragR)
[![R-CMD-check](https://github.com/odaramola92/glycofragR/workflows/R-CMD-check/badge.svg)](https://github.com/odaramola92/glycofragR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R interface to the [glycofrag](https://github.com/odaramola92/glycofrag) Python package that enables glycan structure prediction, theoretical fragment generation, and visualization directly in R using [reticulate](https://rstudio.github.io/reticulate/).

## Features

✨ **Core Capabilities:**
- **Glycan Structure Prediction** - Predict N-glycan and O-glycan structures from mass codes with biosynthetic rules
- **Fragment Ion Generation** - Generate theoretical fragments (B, Y, BY, YY, CZ series) for glycan and glycopeptide MS/MS analysis
- **Mass Calculations** - Accurate monoisotopic mass computation with modification support
- **Peptide Modifications** - 29 modification types (carbamidomethylation, oxidation, phosphorylation, etc.)
- **Visualization** - SNFG-compliant glycan structure diagrams with publication-quality output
- **Modification Support** - Fixed (CAM, TMT, iTRAQ) and variable modifications (Oxidation, Phosphorylation, etc.)

📊 **Fragment Types Supported:**
- **Glycan Fragments**: B, Y, BY, YY, BYY, YYY (BY-series), C, Z, CZ, ZZ (CZ-series), A (oxonium/diagnostic ions)
- **Peptide Fragments**: b, y, c, z ions with configurable charge states
- **Glycopeptide Fragments**: Y0 (peptide only), Y1 (peptide + glycan), Y-glycan (intact), hybrid ions

## Monosaccharide Nomenclature

**Important**: glycofragR distinguishes between monosaccharide types:
- **N-Glycans**: Use **GlcNAc** (N-acetyl-D-glucosamine) as the base monosaccharide
  - Motif: Asparagine (N) in N-X-Ser/Thr context (where X ≠ Pro)
  - Core structure: Two GlcNAc residues linked to asparagine
- **O-Glycans**: Use **GalNAc** (N-acetyl-D-galactosamine) as the base monosaccharide
  - Motif: Directly attached to Serine (S) or Threonine (T)
  - Core structure: Single GalNAc attached to S/T hydroxyl group

## Installation

### Prerequisites

1. **R** (≥ 4.0)
2. **Python** (3.8+) via conda or pip
3. **glycofrag** Python package

### Step 1: Install glycofragR from CRAN (when available)

```r
install.packages("glycofragR")
```

Or install from GitHub:

```r
remotes::install_github("odaramola92/glycofragR")
```

### Step 2: Set up Python environment

Use the one-time installer to create a dedicated environment and install `glycofrag`:

```r
library(glycofragR)
glycofrag_install()
```

**Optional: Manual environment setup**

If you prefer to manage your own environment, you can still configure it explicitly.

**Option A: Using conda**

```r
# In R console:
reticulate::install_miniconda()
reticulate::conda_create("r-glycofrag", packages = "python=3.11")
reticulate::conda_install("r-glycofrag", "glycofrag")

# Configure reticulate to use this conda environment
glycofragR::set_glycofrag_conda_env("r-glycofrag")
```

**Option B: Using pip with existing Python**

```r
# In R console:
reticulate::virtualenv_create("r-glycofrag")
reticulate::virtualenv_install("r-glycofrag", "glycofrag")
glycofragR::set_glycofrag_conda_env("r-glycofrag")
```

## Quick Start

### Basic glycan analysis

```r
library(glycofragR)

# Create glycan from composition code
# Format: HexNAc(n)Hex(n)Fuc(n)NeuAc(n)
glycan <- glycan("4501", "N")

# Predict structures
structures <- glycan_predict_structures(glycan)

# Generate fragments for first structure  
fragments <- glycan_generate_fragments(
  glycan, 
  structures[[1]],
  fragment_types = c("BY", "CZ")
)

# Calculate mass
mass <- glycan_mass("4501")
```

### Advanced: Full glycan analysis

```r
# All-in-one analysis with visualization and Excel export
analysis <- glycan_analysis_new(
  glycan_code = "6632",           # Glycan composition code
  glycan_type = "N",              # N or O glycan
  fragment_types = c("BY", "CZ"), # Fragment series to generate
  visualize_structure = "best"    # Show structure visualization
)

# Access results as native R data frames
theoretical_fragments <- analysis$theoretical_df       # All fragments
clean_fragments <- analysis$clean_df                   # Deduplicated
summary <- analysis$summary_df                         # Summary stats
structures <- analysis$structures_df                   # Structure predictions

# Export to Excel
glycan_analysis_export(analysis, "glycan_results.xlsx")
```

### Glycopeptide Analysis

```r
# Full glycopeptide fragmentation
gp <- glycopeptide(
  peptide_sequence = "NFSTKTR",  # Sequence with N-glycan site at position 1 (N-F-S)
  glycan_code = "4501",
  glycosylation_site = 1,        # Position 1 (1-based indexing)
  glycan_type = "N",
  modifications = list(CAM = c(6))  # Carbamidomethyl on Cys at position 6
)

# Generate all fragment types
fragments <- glycopeptide_fragments(gp)
print(sprintf("Total fragments: %d", nrow(fragments)))

# Get summary
summary <- glycopeptide_summary(gp)
cat("Glycopeptide Summary:\n")
cat(sprintf("Peptide: %s\n", summary$peptide_sequence))
cat(sprintf("Glycan: %s\n", summary$glycan_code))
cat(sprintf("Site: Position %d\n", summary$glycosylation_site))
```

### O-Glycosylation Analysis

```r
# O-glycopeptide with Ser/Thr site
o_glyco <- glycopeptide(
  peptide_sequence = "PEPTSDES",
  glycan_code = "2201",          # Simple O-glycan (HexNAc-1, Hex-1)
  glycosylation_site = 4,        # Ser at position 4 (1-based)
  glycan_type = "O"
)

fragments <- glycopeptide_fragments(o_glyco)

# Extract diagnostic ions
oxonium <- fragments[fragments$type == "A", ]
y0_ions <- fragments[fragments$type == "Y0", ]

cat(sprintf("Diagnostic oxonium ions: %d\n", nrow(oxonium)))
cat(sprintf("Peptide backbone ions (Y0): %d\n", nrow(y0_ions)))
```

## API Reference

### Core Functions

#### `glycan(glycan_code, glycan_type = "N")`

Creates a glycan object for structure prediction and fragmentation.

**Parameters:**
- `glycan_code` (string): Composition code (e.g., "4501" = HexNAc-4, Hex-5, Fuc-0, NeuAc-1)
- `glycan_type` (string): "N" for N-glycan or "O" for O-glycan

**Returns:** Glycan object

**Example:**
```r
glycan <- glycan("4501", "N")
structures <- glycan_predict_structures(glycan)
```

#### `glycan_predict_structures(glycan_object)`

Predicts possible 2D structures from glycan composition.

**Parameters:**
- `glycan_object`: Result from `glycan()`

**Returns:** List of predicted structures

---

#### `glycan_generate_fragments(glycan_object, structure, fragment_types = c("BY", "CZ", "A"))`

Generates theoretical fragment ions for a glycan structure.

**Parameters:**
- `glycan_object`: Glycan object created via `glycan()`
- `structure`: Structure from structure prediction
- `fragment_types`: Vector of fragment series ("BY", "CZ", "A")

**Returns:** Data frame with fragments (m/z, charge, fragment name, composition)

---

#### `peptide(sequence, modifications = NULL)`

Creates a peptide object for fragmentation analysis.

**Parameters:**
- `sequence` (string): Amino acid sequence (single letter code)
- `modifications` (list): Optional modifications like `list(CAM = c(1, 5))`

**Returns:** Peptide object

**Example:**
```r
peptide <- peptide("LCPDCPLLAPLNDSR", modifications = list(CAM = c(1, 5)))
fragments <- peptide_generate_all_fragments(peptide)
```

---

#### `glycopeptide(peptide_sequence, glycan_code, glycosylation_site, glycan_type = "N", modifications = NULL)`

Creates a glycopeptide object combining peptide and glycan.

**Parameters:**
- `peptide_sequence` (string): Amino acid sequence
- `glycan_code` (string): Glycan composition code
- `glycosylation_site` (integer): Position of glycosylation (1-based indexing, where 1 = first amino acid)
- `glycan_type` (string): "N" or "O"
- `modifications` (list): Peptide modifications (1-based positions)

**Returns:** Glycopeptide object

---

#### `glycan_analysis_new(glycan_code, glycan_type = "N", ...)`

All-in-one analysis: generates structures, fragments, and optionally visualizes/exports.

**Parameters:**
- `glycan_code` (string): Glycan composition code
- `glycan_type` (string): "N" or "O"
- `fragment_types` (vector): Fragment series to generate
- `visualize_structure` (logical or string): Show visualization

**Returns:** Analysis list with `$theoretical_df`, `$clean_df`, `$summary_df`, `$structures_df`

## Structure Prediction Algorithms

glycofragR uses the same biosynthetically-informed algorithms as the Python package:

### High Mannose (HexNAc = 2)
- Used for early N-glycans (oligomannose, Man5-Man9)  
- Only mannose (Man) residues on antenna
- Maximum 3 direct branching points

### Hybrid (HexNAc = 3)
- Two distinct arms: mannose arm + antenna arm
- Antenna arm has branching HexNAc with Gal/NeuAc
- NeuAc only on antenna Gal, never on Man arm

### Complex (HexNAc > 4)
- Multi-antennary (bi, tri, etc.)
- Trimannosyl core (Man³)
- All antenna Hex labeled as Gal

### O-Glycans
- No mannose backbone (all Hex = Gal)
- Multiple core types based on HexNAc count
- Simpler structure than N-glycans

## Visualization

SNFG-compliant glycan structure diagrams are generated automatically:

```r
# Create and visualize glycan
glycan <- glycan("4501", "N")
structures <- glycan_predict_structures(glycan)

# Visualize first structure
glycan_visualize(glycan, structures[[1]], 
                 output_path = "glycan_structure.png")
```

**SNFG Monosaccharide Symbols & Colors:**
- Blue square: **GlcNAc** (N-glycan core)
- Yellow square: **GalNAc** (O-glycan base)
- Green circle: **Man** (Mannose - core)
- Yellow circle: **Gal** (Galactose - antenna)
- Magenta diamond: **NeuAc** (Sialic acid)
- Gray diamond: **NeuGc** (Alternative sialic acid)
- Red triangle: **Fuc** (Fucose)

## Peptide Modifications

glycofragR supports **29 peptide modifications** via `mod_string` parameter:

### Fixed Modifications (9)
- **CAM**: Carbamidomethylation (+57.021 Da) on Cys
- **PAM**: Propionamide (+71.037 Da) on Cys
- **Palm**: Palmitoylation (+238.230 Da) on Cys
- **Carbamyl**: Carbamylation (+43.006 Da) on K, N-term
- **TMT6/TMT10/TMT16**: Tandem Mass Tags on K, N-term  
- **iTRAQ4/iTRAQ8**: iTRAQ reagents on K, N-term

### Variable Modifications (20)
- **Ox**: Oxidation (+15.995 Da) on Met
- **Phos**: Phosphorylation (+79.966 Da) on Ser, Thr, Tyr
- **Ac**: Acetylation (+42.011 Da) on Lys, N-term
- **Deam**: Deamidation (+0.984 Da) on Asn, Gln
- **Methyl/DiMethyl/TriMethyl**: Methylation on Lys, Arg
- **HexNAc**: O-GlcNAcylation (+203.079 Da) on Ser, Thr
- Plus: Formyl, Biotin, Malonyl, Succinyl, Farnesyl, SUMO1-GG, nitrosylation, sulfation, etc.

**Usage:**
```r
# Apply modifications via mod_string parameter
peptide <- peptide("MPEPTIDE", mod_string = "M:Ox")

# Or with glycopeptides
gp <- glycopeptide(
  "MKPEPTIDE", "4501", 5, "N",
  mod_string = "M:Ox; K:Ac"
)
```

## Fragment Notation

**Peptide fragments** use lowercase:
- b, y (backbone)
- c, z (alternative cleavage)

**Glycan fragments** use uppercase:
- B, Y (BY-series)
- C, Z (CZ-series)
- A (oxonium ions)

**Glycopeptide hybrid** compounds both:
- Y0 (peptide-only ions)
- Y1 (peptide + glycan)
- Y1-B2, Y1-Y3 (hybrid ions)

## Batch Processing

```r
# Analyze multiple glycopeptides efficiently
peptides <- c("NFSTKTR", "NGSTQQSR", "TPEPTIDES")
glycans <- c("4501", "3410", "2201")
sites <- c(1, 1, 1)

results <- data.frame()
for (i in seq_along(peptides)) {
  gp <- glycopeptide(
    peptides[i],
    glycans[i],
    sites[i],
    "N"
  )
  frags <- glycopeptide_fragments(gp)
  results <- rbind(results, data.frame(
    peptide = peptides[i],
    glycan = glycans[i],
    fragment_count = nrow(frags)
  ))
}

print(results)
```

## Troubleshooting  modification_type = "permethylated_reduced",
  max_structures = 100,           # Predict up to 100 structures
  isomer_sensitive = FALSE,       # Treat isomers as same
  fragment_types = c("BY", "CZ"), # Generate both series
  charges = c(1, 2, 3),          # Charge states for m/z
  visualize_structure = "all"     # Visualize all structures
)

# Access results as native R data.frames
theoretical <- glycan_analysis_theoretical_df(analysis)
clean <- glycan_analysis_clean_df(analysis)
summary <- glycan_analysis_summary(analysis)
structures <- glycan_analysis_structures(analysis)

# Export to Excel
glycan_analysis_export(analysis, "glycan_results.xlsx")
```

### Glycopeptide analysis

```r
# Analyze glycopeptide with attached glycan
gp <- glycopeptide(
  peptide_sequence = "LCPDCPLLAPLNDSR",
  glycan_code = "4501",
  glycosylation_site = 12,        # 1-based position
  glycan_type = "N",
  modifications = list(C = "CAM") # Carbamidomethyl on cysteines
)

# Generate all fragment ions
fragments <- glycopeptide_fragments(gp)

# Get summary info
summary <- glycopeptide_summary(gp)
```

## Available Functions

### Glycan Analysis
- `glycan()` - Create glycan object
- `glycan_predict_structures()` - Predict structures
- `glycan_generate_fragments()` - Generate fragments
- `glycan_mass()` - Calculate mass
- `glycan_visualize()` - Visualize structure

### High-Level Glycan Analysis
- `glycan_analysis_new()` - Complete analysis in one call
- `glycan_analysis_theoretical_df()` - Get all fragments as data.frame
- `glycan_analysis_clean_df()` - Get deduplicated fragments
- `glycan_analysis_summary()` - Get summary data
- `glycan_analysis_structures()` - Get predicted structures
- `glycan_analysis_export()` - Export to Excel
- `glycan_analysis_visualize()` - Visualize structures

### Peptide Analysis
- `peptide()` - Create peptide object
- `peptide_fragments()` - Generate peptide fragments

### Glycopeptide Analysis
- `glycopeptide()` - Create glycopeptide object
- `glycopeptide_fragments()` - Generate all fragments
- `glycopeptide_summary()` - Get summary info

### Utilities
- `glycofrag_install()` - One-time Python setup
- `set_glycofrag_conda_env()` - Configure a custom Python environment
- `list_supported_modifications()` - List available modifications

## Output Formats

All functions return native R objects:
- **data.frames** for tabular results (fragments, structures, summaries)
- **Lists** for nested/complex data
- **Numeric vectors** for masses and m/z values

No Python objects are returned - everything is automatically converted to R-compatible formats.

## Supported Modifications

### Reducing End (glycans)
- `free` - Free reducing end (default)
- `reduced` - Reduced/alditol end
- `permethylated_free` - Permethylated free
- `permethylated_reduced` - Permethylated reduced
- `2ab` - 2-AB labeled
- `2ab_permethylated` - 2-AB + permethylated

### Peptide Modifications
- `C:CAM` - Carbamidomethyl on cysteine
- `M:Ox` - Oxidation on methionine
- `S:Phos` - Phosphorylation on serine
- And 50+ more (see `list_supported_modifications()`)

## Fragment Series

- **BY**: B-ions (C-terminal) and Y-ions (N-terminal)
- **CZ**: C-ions (N-terminal) and Z-ions (C-terminal)
- **Custom diagnostic ions**: Oxonium and marker ions

## Troubleshooting

### "Python module 'glycofrag' not found"

Your Python environment doesn't have glycofrag installed:

```r
library(glycofragR)
glycofrag_install()
```

### "reticulate not available"

Install reticulate if not already installed:

```r
install.packages("reticulate")
library(reticulate)
```

### Python environment issues

```r
# Check current Python path
reticulate::py_available()
reticulate::py_config()

# List available conda environments
reticulate::conda_list()

# Set to specific environment
set_glycofrag_conda_env("r-glycofrag")  # or Sys.setenv(RETICULATE_PYTHON_ENV = "r-glycofrag")
```

## Citation

### Methodology Paper

If you use glycofragR in your research, please cite the GlypPRM methodology paper:

```bibtex
@article{daramola2026glypprm,
  title={{GlypPRM: An Automated Analyzer and Quantification Tool for Glycopeptides Parallel Reaction Monitoring}},
  author={Daramola, Oluwatosin and Onigbinde, Sherifdeen and Adeniyi, Moyinoluwa and Gutierrez-Reyes, Cristian D. and Fowowe, Mojibola and Sandilya, Vishal and Mechref, Yehia},
  journal={Analytical Chemistry},
  volume={98},
  number={2},
  pages={1528--1540},
  year={2026},
  publisher={ACS Publications},
  doi={10.1021/acs.analchem.5c06005}
}
```

### Software Citation (Optional but Recommended)

To cite the software packages, we recommend registering them on [Zenodo](https://zenodo.org) for persistent DOIs. This enables proper attribution of the software implementations separately from the methodology paper.

#### How to Generate a Citation

**In R:**
```r
citation("glycofragR")
```

**From a CITATION.cff file:**
Both glycofragR and the underlying glycofrag Python package include `CITATION.cff` files (Citation File Format) that can be used by reference management tools.

#### Zenodo Registration (Coming Soon)

Once registered on Zenodo, the packages will have DOIs like:

```bibtex
@software{daramola2026glycofragR,
  title={glycofragR: R Interface to Glycan Fragmentation Analysis},
  author={Daramola, Oluwatosin},
  year={2026},
  url={https://github.com/odaramola92/glycofragR},
  doi={10.5281/zenodo.XXXXXXX(R)}
}

@software{daramola2026glycofrag,
  title={glycofrag: Glycan and Glycopeptide Fragmentation Analysis},
  author={Daramola, Oluwatosin},
  year={2026},
  url={https://github.com/odaramola92/glycofrag},
  doi={10.5281/zenodo.XXXXXXX(Python)}
}
```

**Note:** glycofragR is an R interface to the glycofrag Python package, which is a separate implementation of the GlypPRM algorithms. All three (methodology paper, Python package, and R package) should ideally be cited to give proper credit.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Commit (`git commit -m 'Add amazing feature'`)
5. Push to branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## Support

For issues, questions, or suggestions:
- Open an issue on [GitHub Issues](https://github.com/odaramola92/glycofragR/issues)
- Check the [Python glycofrag documentation](https://github.com/odaramola92/glycofrag/README.md)
- Review test cases and examples in the package

## Acknowledgments

- R wrapper implementation for [glycofrag](https://github.com/odaramola92/glycofrag) Python package
- Uses [reticulate](https://rstudio.github.io/reticulate/) for Python/R interoperability
- SNFG nomenclature for glycan visualization
- Monosaccharide masses from standard biochemistry references

## Changelog

### Version 0.1.0 (2026-02-08)
- ✅ Initial release
- ✅ Core glycan functions (new, predict, generate_fragments)
- ✅ Peptide fragmentation interface
- ✅ Glycopeptide analysis (complete pipeline)
- ✅ SNFG visualization support
- ✅ 29 modification types supported
- ✅ Comprehensive R documentation
- ✅ Integration tests with Python glycofrag (354 tests total)

---

**Last Updated**: February 9, 2026  
**Version**: 0.1.0  
**R Requirement**: ≥ 4.0  
**Status**: Production Ready ✅
