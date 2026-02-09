# Glycofrag - Glycan and Glycopeptide Fragmentation Analysis

A comprehensive Python package for analyzing glycan structures, peptide fragmentation, and glycopeptide mass spectrometry data.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests: 193/193 ✅](https://img.shields.io/badge/tests-193%2F193%20%E2%9C%85-brightgreen)](https://github.com/glypprm/glycofrag)

## Features

✨ **Core Capabilities:**
- **Glycan Mass Calculation** - Accurate mass computation for various glycan compositions
- **Glycan Structure Prediction** - Predict N-glycan and O-glycan structures from mass codes
- **Glycan Fragmentation** - Generate theoretical fragment ions (B, Y, BY, YY, CZ series)
- **Peptide Fragmentation** - Generate peptide fragment ions (b, y, c, z series)
- **Glycopeptide Analysis** - Complete analysis of glycosylated peptides
- **Modification Support** - Carbamidomethylation, oxidation, permethylation, and more
- **Reducing End Modifications** - Support for 6 different reducing end types

📊 **Fragment Types:**
- **Glycan Fragments**: B, Y, BY, YY, BYY, YYY, BYYY (BY-series), C, Z, CZ, ZZ, CZZ, BZZ, ZZZ, CZZZ, BZZZ (CZ-series)
- **Oxonium Ions**: A-type diagnostic ions for O-glycan and N-glycan detection
- **Peptide Fragments**: b, y, c, z ions with configurable charge states
- **Glycopeptide Fragments**: Y0 (peptide only), Y1 (peptide + glycan), Y-glycan (peptide + glycan fragments), intact ions

## Installation

```bash
# From PyPI (recommended)
pip install glycofrag

# From source (development)
git clone https://github.com/odaramola92/glycofrag.git
cd glycofrag
pip install -e .

# Requirements: Python 3.8+, numpy, pandas, networkx, matplotlib
```

## Quick Start

### 1. Basic Glycan Analysis

```python
from glycofrag import Glycan, GlycanMassCalculator

# Calculate glycan mass
calc = GlycanMassCalculator()
mass = calc.calculate_glycan_mass("4501")
print(f"Glycan mass: {mass:.4f} Da")

# Predict glycan structures
glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()
print(f"Found {len(structures)} possible N-glycan structures")

# Generate fragments
fragments, info = glycan.generate_fragments(structures[0])
print(f"Generated {len(fragments)} fragment ions")
for fragment in fragments[:5]:
    print(f"  {fragment['name']}: {fragment['mass']:.4f} m/z at charge {fragment['charge']}")
```

### 2. Peptide Fragmentation

```python
from glycofrag import Peptide

# Create peptide object (positions are 1-based)
peptide = Peptide("LCPDCPLLAPLNDSR", modifications={'CAM': [1, 5]})  # L at position 1, C at position 5

# Generate all fragments
fragments = peptide.generate_all_fragments()
print(f"Generated {len(fragments)} peptide fragments")

# Get specific fragment type
b_ions = [f for f in fragments if f['type'] == 'b']
y_ions = [f for f in fragments if f['type'] == 'y']
print(f"b-ions: {len(b_ions)}, y-ions: {len(y_ions)}")
```

### 3. Complete Glycopeptide Analysis

```python
from glycofrag import Glycopeptide

# Create glycopeptide object
glycopeptide = Glycopeptide(
    peptide_sequence="LCPDCPLLAPLNDSR",
    glycan_code="4501",
    glycosylation_site=12,  # N at position 12 (1-based indexing)
    glycan_type="N",
    modifications={'CAM': [2, 5]}  # Carbamidomethylation at 2 and 5
)

# Generate all fragment types
fragments = glycopeptide.generate_fragments()
print(f"Total fragments: {len(fragments)}")

# Analyze fragment types
fragment_types = {}
for frag in fragments:
    ftype = frag['type']
    fragment_types[ftype] = fragment_types.get(ftype, 0) + 1

print("\nFragment distribution:")
for ftype, count in sorted(fragment_types.items()):
    print(f"  {ftype}: {count}")

# Get specific ions
y0_ions = [f for f in fragments if f['type'] == 'Y0']
y1_ions = [f for f in fragments if f['type'] == 'Y1']
print(f"\nY0 ions (peptide only): {len(y0_ions)}")
print(f"Y1 ions (peptide + glycan): {len(y1_ions)}")
```

## API Reference

### Core Classes

#### `Glycan(glycan_code, glycan_type="N")`

Represents a glycan structure and provides analysis methods.

**Parameters:**
- `glycan_code` (str): Composition code (e.g., "4501" for HexNAc-4, Hex-5, Fuc-0, NeuAc-1)
- `glycan_type` (str): Type of glycan - "N" or "O" (default: "N")

**Methods:**
- `predict_structures()` → list: Predict possible structures from composition
- `generate_fragments(structure, modification_type=0, permethylated=False, fragment_types=['BY','CZ','A'])` → tuple: Generate fragment ions

**Example:**
```python
glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()
if structures:
    fragments, metadata = glycan.generate_fragments(structures[0])
```

---

#### `Peptide(sequence, modifications=None)`

Represents a peptide sequence and provides fragmentation analysis.

**Parameters:**
- `sequence` (str): Amino acid sequence (single letter code)
- `modifications` (dict, optional): Modification dictionary with amino acid positions (1-based indexing)
  - Format: `{'CAM': [1, 3, 5], 'Oxidation': [2]}`  # Position 1 = first amino acid
  - Supported: CAM (carbamidomethylation), Oxidation

**Methods:**
- `generate_all_fragments(charge_states=[1, 2])` → list: Generate all ion types
- `get_fragment_by_name(name)` → dict: Retrieve specific fragment by name

**Example:**
```python
peptide = Peptide("MEPEPTIDE")
fragments = peptide.generate_all_fragments(charge_states=[1, 2, 3])
```

---

#### `Glycopeptide(peptide_sequence, glycan_code, glycosylation_site, glycan_type="N", modifications=None)`

Represents a glycosylated peptide with comprehensive fragmentation analysis.

**Parameters:**
- `peptide_sequence` (str): Amino acid sequence
- `glycan_code` (str): Glycan composition code
- `glycosylation_site` (int): Position of glycosylation in peptide (1-based indexing, where position 1 = first amino acid)
- `glycan_type` (str): "N" for N-glycan, "O" for O-glycan (default: "N")
- `modifications` (dict, optional): Peptide modifications (positions are 1-based)

**Fragment Types Generated:**
- **Peptide fragments**: b, y, c, z ions (backbone only)
- **Glycan fragments**: B, Y, BY, YY, C, Z ions (all with -PEP suffix)
- **Intact**: Complete glycopeptide (peptide + glycan)

**Methods:**
- `generate_fragments()` → list: Generate all fragment ions

**Example:**
```python
gp = Glycopeptide("EPEPTIDES", "4501", 2, glycan_type="N")  # Glycosylation at position 2 (P = second amino acid)
fragments = gp.generate_fragments()
```

---

#### `GlycanMassCalculator()`

Low-level mass calculation utilities for glycan components.

**Methods:**
- `calculate_glycan_mass(code)` → float: Calculate mass from composition code
- `calculate_monosaccharide_mass(name, count)` → float: Calculate mass of monosaccharide

**Example:**
```python
calc = GlycanMassCalculator()
mass = calc.calculate_glycan_mass("4501")  # HexNAc-4, Hex-5, Fuc-0, NeuAc-1
```

---

## Structure Prediction Algorithms

Glycofrag uses biosynthetically-informed algorithms to predict glycan structures from composition alone. Each glycan type uses different rules to ensure realistic structure generation.

### High Mannose (HexNAc = 2)

**Rules:**
- **Core constraint**: Maximum 3 total direct children across nodes 4 and 5 combined
- **Fill strategy**: Depth-first filling with minimal children balancing
- **All Hex labeled as Man** (mannose - no branching HexNAc allowed)
- Used for early-stage N-glycans (oligomannose, Man5 through Man9)

**Example:**
```python
# High Mannose: HexNAc(2)Man(5) - captured as code "2500"
glycan = Glycan("2500", glycan_type="N")
structures = glycan.predict_structures()
# All structures will have only Man residues, no branching
```

### Hybrid (HexNAc = 3)

**Rules:**
- **Two distinct arms**:
  - Mannose arm: unbranched Man chain from core (nodes 4 or 5)
  - Antenna arm: branching HexNAc with Gal → NeuAc architecture
- **NeuAc placement**: Only on antenna Gal, never on mannose arm Man
- **Gal reservation**: When NeuAc_count > antenna_Gal_count, only branch HexNAc receives additional Gal residues
  - Ensures biosynthetic realism (each NeuAc typically paired with Gal)
- **Core Hex labeled as Man**, antenna Hex labeled as Gal

**Example:**
```python
# Hybrid with 2 NeuAc: HexNAc(3)Hex(6)NeuAc(2) - code "3602"
glycan = Glycan("3602", glycan_type="N")
structures = glycan.predict_structures()
# Structures show clear separation between mannose and antenna arms
# Each NeuAc attached to antenna Gal only
```

### Complex (HexNAc > 3)

**Rules:**
- **Multi-antennary**: Multiple branching HexNAc arms radiate from core
- **Core nodes 3, 4, 5**: Labeled as Man (trimannosyl core)
- **All branching Hex**: Labeled as Gal (no mannose on antenna arms)
- **Supports**: Bi-antennary (HexNAc=4), tri-antennary (HexNAc=5+), etc.

**Example:**
```python
# Complex bi-antennary: HexNAc(4)Hex(4)NeuAc(2) - code "4402"
glycan = Glycan("4402", glycan_type="N")
structures = glycan.predict_structures()
# Typical N-glycan serum proteins (IgG, transferrin, etc.)
```

### O-Glycans (vs N-Glycans)

**Rules:**
- **No mannose backbone** - all Hex labeled as Gal (galactose)
- **Multiple core types** (Core 0, 1, 2, 3) determined by available HexNAc
- **No tri-mannosyl core** - O-glycans lack the Man₃ architecture

**Example:**
```python
# O-glycan Core 1: HexNAc(1)Hex(1) - code "1100"
glycan = Glycan("1100", glycan_type="O")
structures = glycan.predict_structures()
# Simpler T antigen structure (GalNAc-Gal)
```

### Visualization of Predictions

Structures can be visualized automatically with SNFG-compliant symbols:

```python
from glycofrag import Glycan
from glycofrag.io.visualizer import GlycanVisualizer

glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()

# Visualize first structure
GlycanVisualizer.visualize(
    structures[0],
    title="N-Glycan Structure 1",
    output_path="structure.png",
    show=True
)

# View as ASCII tree for terminal output
print(glycan.format_structure_tree(structures[0]))
```

**Output Colors (SNFG Standard):**
- 🟦 Blue box: HexNAc
- 🟢 Green circle: Man (mannose, core)
- 🟡 Yellow circle: Gal (galactose, antenna)
- 🔴 Red triangle: Fuc (fucose)
- 🟣 Magenta diamond: NeuAc (sialic acid)

---

### Constants and Masses

Available in `glycofrag.core.constants`:

```python
from glycofrag import MONOSACCHARIDE_MASSES, AMINO_ACID_MASSES, MODIFICATION_MASSES

# Access standard masses
hex_mass = MONOSACCHARIDE_MASSES['Hex']
ala_mass = AMINO_ACID_MASSES['A']
cam_mass = MODIFICATION_MASSES['CAM']
```

## Peptide Modification Support

Glycofrag supports **29 peptide modifications** with targeted or global application using the `mod_string` parameter.

### Quick Start

```python
from glycofrag import Peptide, Glycopeptide, list_supported_modifications

# List all supported modifications
mods = list_supported_modifications(verbose=True)

# Apply modifications to peptides
peptide = Peptide("PEPTIDE", mod_string="M:Ox")  # Oxidize all Met
peptide = Peptide("PEPTIDE", mod_string="M4,24:Ox")  # Oxidize Met at positions 4, 24
peptide = Peptide("PEPTIDE", mod_string="M:Ox; S3:Phos")  # Multiple mods

# Apply to glycopeptides
gp = Glycopeptide("EEQYNSTYR", "4501", 5, "N", mod_string="M:Ox")
```

### Supported Modifications (29 total)

#### Fixed Modifications (9)
Typically applied to all occurrences of target residues:

| Modification | Mass (Da) | Targets | Description |
|-------------|-----------|---------|-------------|
| **CAM** | +57.021 | C | Carbamidomethylation (iodoacetamide alkylation) |
| **PAM** | +71.037 | C | Propionamide (acrylamide alkylation) |
| **Palm** | +238.230 | C | Palmitoylation (lipid modification) |
| **Carbamyl** | +43.006 | K, N-term | Carbamylation |
| **TMT6** | +229.163 | K, N-term | Tandem Mass Tag 6-plex |
| **TMT10** | +229.163 | K, N-term | Tandem Mass Tag 10-plex |
| **TMT16** | +304.207 | K, N-term | TMTpro 16-plex |
| **iTRAQ4** | +144.102 | K, N-term | iTRAQ 4-plex reagent |
| **iTRAQ8** | +304.205 | K, N-term | iTRAQ 8-plex reagent |

#### Variable Modifications (20)
Applied selectively based on MS/MS evidence:

| Modification | Mass (Da) | Targets | Description |
|-------------|-----------|---------|-------------|
| **Ox** | +15.995 | M | Oxidation (methionine sulfoxide) |
| **Phos** | +79.966 | S, T, Y | Phosphorylation |
| **Ac** | +42.011 | K, N-term | Acetylation |
| **Deam** | +0.984 | N, Q | Deamidation (Asn→Asp, Gln→Glu) |
| **Methyl** | +14.016 | K, R | Monomethylation |
| **DiMethyl** | +28.031 | K, R | Dimethylation |
| **TriMethyl** | +42.047 | K | Trimethylation |
| **Formyl** | +27.995 | K, N-term | Formylation |
| **Biotin** | +226.078 | K | Biotinylation (affinity tag) |
| **Malonyl** | +86.000 | K | Malonylation |
| **Succinyl** | +100.016 | K | Succinylation |
| **Nitration** | +44.985 | Y | Nitration (nitrotyrosine) |
| **Sulf** | +79.957 | Y | Sulfation (tyrosine sulfate) |
| **GG** | +114.043 | K | Ubiquitin remnant (GlyGly) |
| **HexNAc** | +203.079 | S, T | O-GlcNAcylation |
| **Pyro-glu** | -17.027 | N-term-Q | Pyroglutamic acid from Gln |
| **Pyro-cmC** | -17.027 | N-term-C-CAM | Pyroglutamic acid from CAM-Cys |
| **Myristoyl** | +210.198 | N-term-G | N-terminal myristoylation |
| **Farnesyl** | +204.188 | C | Farnesylation (prenylation) |
| **SUMO1-GG** | +213.144 | K | SUMOylation remnant |

### Modification String Formats

The `mod_string` parameter supports multiple formats:

```python
# Format 1: Apply to all occurrences of amino acid
mod_string = "M:Ox"  # Oxidize ALL methionines

# Format 2: Target specific positions (1-based indexing)
mod_string = "M4,24:Ox"  # Oxidize Met at positions 4 and 24 (4th and 24th amino acids)

# Format 3: Multiple modifications (semicolon-separated)
mod_string = "M:Ox; S3:Phos; K:Ac"  # Multiple mods

# Format 4: Custom mass modification
mod_string = "Custom:A, 45.8978"  # Add +45.8978 Da to all Alanines

# Format 5: Named modifications with positions
mod_string = "Ox:M6; Phos:S3,T15"  # Explicit position syntax
```

### Usage Examples

```python
from glycofrag import Peptide, Glycopeptide

# Example 1: Simple oxidation
peptide = Peptide("MPEPTIDE", mod_string="M:Ox")
print(f"Mass with oxidation: {peptide.mass:.4f} Da")

# Example 2: Phosphorylation at specific sites (1-based positions)
peptide = Peptide("KESMKDES", mod_string="S3:Phos; S8:Phos")  # S at position 3 and S at position 8
b_ions = peptide.generate_b_ions()

# Example 3: Multiple modification types
peptide = Peptide("MKPEPTIDE", mod_string="M:Ox; K:Ac")

# Example 4: Glycopeptide with modifications
gp = Glycopeptide(
    peptide_sequence="MKPEPTIDE",
    glycan_code="4501",
    glycosylation_site=5,  # P at position 5 (1-based)
    glycan_type="N",
    mod_string="M:Ox; K:Ac"
)
fragments = gp.generate_fragments()

# Example 5: TMT-labeled peptide
peptide = Peptide("PEPTIDEK", mod_string="K:TMT10")  # TMT on lysine and N-term
```

### Programmatic Access

Query modifications programmatically:

```python
from glycofrag import MODIFICATION_MASSES, MODIFICATION_TARGETS, list_supported_modifications

# Get all modifications as dict
mods = list_supported_modifications(verbose=False)

# Access specific modification
ox_mass = MODIFICATION_MASSES['Ox']
ox_targets = MODIFICATION_TARGETS['Ox']
print(f"Oxidation: {ox_mass} Da on {ox_targets}")

# Check if modification is supported
if 'Phos' in MODIFICATION_MASSES:
    print(f"Phosphorylation mass: {MODIFICATION_MASSES['Phos']} Da")
```

### Reducing End Modifications (for glycans)
- **Type 0**: Free reducing end (H₂O)
- **Type 1**: Reduced end (alditol form)
- **Type 2**: Permethylated free end
- **Type 3**: Permethylated reduced end
- **Type 4**: 2-AB labeled
- **Type 5**: 2-AB labeled and permethylated
- **Type 6**: Glycopeptide mode (peptide replaces reducing end)

## Fragment Ion Information

Each fragment ion contains:
- `mass` (float): m/z value
- `charge` (int): Charge state
- `name` (str): Fragment designation (e.g., "B1", "Y3", "Y1-b2")
- `type` (str): Fragment category (e.g., "B", "Y", "BY")
- `composition` (dict): Monosaccharide composition (for glycan fragments)
- `metadata` (dict): Additional information (modifications, source, etc.)

## Visualization Features

Glycofrag provides publication-quality glycan structure visualization with intelligent layout and SNFG-compliant symbols.

### Automatic Layout

- **Hierarchical tree layout**: Root HexNAc centered, branches spread downward
- **Dynamic branch spacing**: Automatically adjusted based on total branches
  - High Mannose (2 branches): Standard spacing
  - Hybrid (4 branches): Increased spacing for clarity (3.9 units between core arms)
  - Complex (3+ antenna): Adaptive spacing
- **Intelligent node positioning**: Prevents overlaps while preserving biosynthetic relationships

### Fucose Positioning Rules

- **Core HexNAc (node 1)**: Straight horizontal alignment
- **Far-left/far-right branches**: Straight, pointing outward (left or right respectively)
- **Middle branches**: Bent upward for unobstructed visibility
- **Branching HexNAc nodes**: Straight, pointing away from center

### Visualization Example

```python
from glycofrag import Glycan
from glycofrag.io.visualizer import GlycanVisualizer

# Predict structures
glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()

# Method 1: Save to file
GlycanVisualizer.visualize(
    structures[0],
    title="N-Glycan Complex Bi-antennary",
    output_path="glycan_structure.png",
    show=False  # Don't open in browser
)

# Method 2: Get ASCII tree for terminal
tree = glycan.format_structure_tree(structures[0])
print(tree)
```

### Customizing Layout Parameters

Edit `GlycanVisualizer` class constants to customize spacing:

```python
# Larger spacing for better readability
GlycanVisualizer.CORE_BRANCH_SPACING_INCREMENT = 1.2  # default: 1.0

# More vertical separation
GlycanVisualizer.VERTICAL_GAP = 1.5  # default: 1.2

# Adjust fucose bend angle
GlycanVisualizer.FUC_VERTICAL_OFFSET = -1.2  # default: -0.9 (more negative = sharper bend)
```

## Usage Examples

### Example 1: Compare Experimental vs Theoretical Fragments

```python
from glycofrag import Glycopeptide

# Define your glycopeptide
gp = Glycopeptide(
    "EPTETTQSVK",      # Peptide sequence
    "4501",             # N-glycan composition
    2,                  # Glycan attached at position 2 (P = second amino acid, 1-based)
    glycan_type="N"
)

# Generate theoretical fragments
theoretical = gp.generate_fragments()

# Compare with experimental m/z values
experimental_mz = [300.5, 400.3, 520.8, 1200.4]

print("Experimental vs Theoretical Matches:")
for exp_mz in experimental_mz:
    matches = [f for f in theoretical if abs(f['mass'] - exp_mz) < 0.05]
    if matches:
        for match in matches[:3]:
            print(f"  {exp_mz:.1f} m/z → {match['name']} ({match['mass']:.4f} m/z)")
    else:
        print(f"  {exp_mz:.1f} m/z → No match found")
```

### Example 2: Analyze O-Glycosylation

```python
from glycofrag import Glycopeptide

# O-glycopeptide with core 1 structure
o_glyco = Glycopeptide(
    "TPEPTIDE",        # Note: Threonine at position 1 (1-based indexing)
    "2201",            # Simple O-glycan (HexNAc-2, Hex-2, Fuc-0)
    1,                 # Glycosylated at position 1 (T = first amino acid)
    glycan_type="O"
)

fragments = o_glyco.generate_fragments()

# Extract diagnostic ions
oxonium = [f for f in fragments if f['type'] == 'A']
y0 = [f for f in fragments if f['type'] == 'Y0']

print(f"Diagnostic oxonium ions: {len(oxonium)}")
print(f"Peptide backbone ions: {len(y0)}")
```

### Example 3: Batch Processing Multiple Glycopeptides

```python
from glycofrag import Glycopeptide

# Define multiple glycopeptides (positions are 1-based)
glycoptides = [
    ("EPTETTQSVK", "4501", 2),  # Position 2 = P
    ("EPTETTQSVK", "3410", 2),  # Position 2 = P
    ("TPEPTIDES", "2201", 1),   # Position 1 = T (changed from 0 to 1-based)
]

results = []
for seq, code, pos in glycopeptides:
    gp = Glycopeptide(seq, code, pos)
    fragments = gp.generate_fragments()
    results.append({
        'sequence': seq,
        'glycan': code,
        'fragment_count': len(fragments),
        'mass': sum(f['mass'] for f in fragments if f['type'] == 'Intact')
    })

# Display results
for result in results:
    print(f"{result['sequence']} + {result['glycan']}: {result['fragment_count']} fragments")
```

## Testing

Run the test suite to verify functionality:

```bash
# All tests
pytest glycofrag/tests/

# Specific module
pytest glycofrag/tests/test_glycopeptide.py -v

# With coverage
pytest glycofrag/tests/ --cov=glycofrag --cov-report=html
```

**Test Coverage:**
- Phase 1 (Core): 48 tests
- Phase 2 (Glycan Structure): 22 tests
- Phase 2.5 (Glycan Fragmentation): 58 tests
- Phase 3 (Peptide): 34 tests
- Phase 4 (Glycopeptide): 31 tests
- **Total: 193 tests** ✅

## Citation

If you use Glycofrag in your research, please cite:

```bibtex
@article{daramola2026glypprm,
  title={{GlypPRM: Structural Elucidation of N-Glycans by Integrating Topology Analysis with Molecular Fragment Ions}},
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

**Plain Text Citation:**
> Daramola, O.; Onigbinde, S.; Adeniyi, M.; Gutierrez-Reyes, C. D.; Fowowe, M.; Sandilya, V.; Mechref, Y. "GlypPRM: Structural Elucidation of N-Glycans by Integrating Topology Analysis with Molecular Fragment Ions." *Analytical Chemistry* **2026**, *98* (2), 1528-1540. DOI: [10.1021/acs.analchem.5c06005](https://doi.org/10.1021/acs.analchem.5c06005)

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Add tests for new functionality
4. Ensure all tests pass (`pytest`)
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Based on the GlypPRM algorithm
- Built with contributions from the glycoproteomics community
- Mass data sourced from NIST and published literature

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check existing documentation
- Review test cases for usage examples

## Changelog

### Version 0.1.0 (2026-02-01)
- ✅ Initial release
- ✅ Core module with mass calculations
- ✅ Glycan structure prediction (N and O)
- ✅ Glycan fragmentation (BY, CZ, A-type series)
- ✅ Peptide fragmentation (b, y, c, z ions)
- ✅ Glycopeptide analysis (Y0, Y1, Y-glycan, intact)
- ✅ 193 unit tests
- ✅ Comprehensive documentation

---

**Last Updated**: February 8, 2026  
**Package Version**: 0.1.0  
**Status**: Production Ready ✅
