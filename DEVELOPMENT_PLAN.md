# GlycoFrag Development Plan

## Project Overview
**Package Name:** `glycofrag`  
**Purpose:** A Python package for glycan and glycopeptide fragmentation analysis  
**Source:** Extracted and refactored from GlypPRM_v01.py  
**Start Date:** February 1, 2026  
**Status:** 🟡 In Progress (Phase 1 Complete ✅ | Phase 2 Complete ✅ | Phase 3 Pending)

---

## Executive Summary

**Phase 1 Completion Status:** ✅ **COMPLETE**
- Package structure: ✅ Created
- Constants module: ✅ 900+ lines, all masses validated
- Mass calculator: ✅ 700+ lines, full GlycanMassCalculator refactored
- Modifications: ✅ 400+ lines, all modification types supported
- Utilities: ✅ Caching and validation modules
- Tests: ✅ 48 test cases, all passing (0.56s execution)
- Package config: ✅ pyproject.toml complete with metadata, dependencies, tool configs
- Documentation: ✅ README.md with comprehensive usage examples

**Phase 2 Completion Status:** ✅ **COMPLETE**
- Glycan module: ✅ Created glycofrag/glycan/ with __init__.py
- Glycan class: ✅ 700+ lines with structure prediction methods
- N-glycan core building: ✅ Complete implementation
- O-glycan core building: ✅ Multiple core types supported
- Structure prediction: ✅ Recursive building algorithms implemented
- Structure fingerprinting: ✅ Deduplication and canonicalization
- Tests: ✅ 22 test cases, all passing (5.19s execution)
- Documentation: ✅ README updated with Glycan examples

**Installation Status:** ✅ **SUCCESS**
- Package installed in development mode: `pip install -e .`
- All dependencies resolved: numpy, pandas, networkx, scipy
- Development tools: pytest, black, mypy, isort installed
- Test suite: **70/70 tests passing ✅** (48 Phase 1 + 22 Phase 2)

**Next Phase:** Phase 3 - Peptide Module & Fragmentation

---

## Package Structure

```
glycofrag/
├── __init__.py                    # Package exports and version
├── core/
│   ├── __init__.py
│   ├── constants.py               # Mass constants, amino acids
│   ├── mass_calculator.py         # GlycanMassCalculator class
│   └── modifications.py           # Modification parsing and handling
├── glycan/
│   ├── __init__.py
│   ├── glycan.py                  # Glycan class (structure + fragmentation)
│   ├── structure.py               # Structure prediction algorithms
│   └── fragmentation.py           # Glycan fragmentation logic
├── peptide/
│   ├── __init__.py
│   ├── peptide.py                 # Peptide class
│   └── fragmentation.py           # Peptide b/y/c/z ion generation
├── glycopeptide/
│   ├── __init__.py
│   ├── glycopeptide.py            # Glycopeptide class
│   └── fragmentation.py           # Combined fragmentation
├── utils/
│   ├── __init__.py
│   ├── caching.py                 # Cache key generation utilities
│   └── validation.py              # Input validation helpers
├── io/
│   ├── __init__.py
│   └── tables.py                  # Fragment table generation (DataFrame)
├── tests/
│   ├── __init__.py
│   ├── test_constants.py
│   ├── test_mass_calculator.py
│   ├── test_glycan.py
│   ├── test_peptide.py
│   └── test_glycopeptide.py
├── pyproject.toml                 # Package configuration
├── README.md                      # Package documentation
└── LICENSE                        # MIT License
```

---

## Phase 1: Core Module - Constants & Mass Calculator
**Target:** Extract fundamental mass constants and the mass calculation engine

### Task 1.1: Create Package Structure
- [x] Create `glycofrag/` directory
- [x] Create `glycofrag/__init__.py`
- [x] Create `glycofrag/core/` directory
- [x] Create `glycofrag/core/__init__.py`

### Task 1.2: Extract Constants (`core/constants.py`)
- [x] Create `core/constants.py`
- [x] Extract `BASE_MASSES` (monosaccharide masses)
- [x] Extract `OH_GROUPS` (hydroxyl groups per monosaccharide)
- [x] Extract `MOD_MASSES` (modification masses)
- [x] Extract `REDUCING_END_MASSES` (reducing end options)
- [x] Extract `ADDITIONAL_MODIFICATIONS`
- [x] Extract `PROTON_MASS`
- [x] Extract `AMINO_ACID_MASSES` dictionary
- [x] Extract `MODIFICATION_MASSES` dictionary (CAM, Ox, Phos, etc.)
- [x] Extract `MODIFICATION_TARGETS` dictionary
- [x] Add docstrings for all constants

**Source Lines in GlypPRM_v01.py:**
- Lines 114-166: Monosaccharide and reducing end masses
- Lines 335-371: Amino acid masses
- Lines 415-497: Modification masses and targets

### Task 1.3: Extract Mass Calculator (`core/mass_calculator.py`)
- [x] Create `core/mass_calculator.py`
- [x] Extract `GlycanMassCalculator` class
- [x] Clean `__init__` method
- [x] Extract `calculate_MONO_MASSES()` method
- [x] Extract `parse_glycan_code()` method
- [x] Extract `calculate_mz()` method
- [x] Extract `get_amino_acid_mass()` method
- [x] Extract `calculate_peptide_mass()` method
- [x] Extract `calculate_fragment_mass()` method
- [x] Extract `_calculate_base_fragment_mass()` method
- [x] Extract `_create_peptide_cache_key()` method
- [x] Remove debug print statements
- [x] Remove GUI-specific code
- [x] Add comprehensive docstrings
- [x] Add type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 109-1040: Full GlycanMassCalculator class

### Task 1.4: Extract Modifications (`core/modifications.py`)
- [x] Create `core/modifications.py`
- [x] Extract `parse_modification_string()` function
- [x] Extract `parse_custom_mass_mod()` function
- [x] Create `apply_modifications()` helper function
- [x] Add docstrings and type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 1149-1218: Modification parsing functions

### Task 1.5: Core Module Tests
- [x] Create `tests/test_constants.py`
- [x] Create `tests/test_mass_calculator.py`
- [x] Test glycan mass calculation
- [x] Test peptide mass calculation
- [x] Test modification parsing
- [x] Test m/z calculation

### Task 1.6: Package Configuration & Documentation
- [x] Create `pyproject.toml` with metadata, dependencies, tool configs
- [x] Create `README.md` with usage examples and API reference
- [x] Create `DEVELOPMENT_PLAN.md` master roadmap
- [x] Install package in development mode: `pip install -e .`
- [x] Install development dependencies: pytest, black, mypy, isort
- [x] Run test suite: **48/48 tests passing ✅** (0.56s execution)

**Phase 1 Summary:**
- Lines of code written: ~4,500 (constants, calculator, modifications, utils, tests)
- Files created: 8 Python modules + 3 config/doc files
- Test coverage: 48 test cases across 2 files
- Installation status: ✅ Success - all dependencies resolved
- Ready for Phase 2: ✅ Yes

---

## Phase 2: Glycan Module - Structure & Fragmentation
**Target:** Extract glycan structure prediction and fragmentation

### Task 2.1: Create Glycan Module Structure
- [ ] Create `glycofrag/glycan/` directory
- [ ] Create `glycofrag/glycan/__init__.py`

### Task 2.2: Extract Glycan Class (`glycan/glycan.py`)
### Task 2.1: Create Glycan Module Structure
- [x] Create `glycofrag/glycan/` directory
- [x] Create `glycofrag/glycan/__init__.py`

### Task 2.2: Extract Glycan Class (`glycan/glycan.py`)
- [x] Create `glycan/glycan.py`
- [x] Extract `Glycan` class initialization
- [x] Extract `_build_n_glycan_core()` method
- [x] Extract `_build_o_glycan_core()` method
- [x] Extract `predict_structures()` method
- [x] Extract `_build_structures()` method (N-glycan)
- [x] Extract `_build_complete_o_glycan_structures()` method
- [x] Extract `classify_structure()` method
- [x] Extract `_structure_fingerprint()` method
- [x] Extract `_get_subtree_hash()` method
- [x] Add docstrings and type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 1222-2100: Glycan class

### Task 2.3: Glycan Module Tests
- [x] Create `tests/test_glycan.py`
- [x] Test glycan initialization
- [x] Test N-glycan core building
- [x] Test O-glycan core building
- [x] Test structure prediction
- [x] Test structure classification
- [x] Test structure fingerprinting
- [x] Test residue counting

### Task 2.4: Update Package Exports
- [x] Update `glycofrag/__init__.py` to export Glycan class
- [x] Update `README.md` with Glycan usage examples
- [x] Test package imports

**Phase 2 Summary:**
- Lines of code written: ~700 (glycan class + tests)
- Files created: 2 (glycan.py, test_glycan.py)
- Test coverage: 22 test cases, all passing ✅
- Total test suite: 70 tests passing (48 Phase 1 + 22 Phase 2)
- Execution time: 6.10 seconds

---

## Phase 3: Peptide Module - Fragmentation (FUTURE)
**Target:** Extract peptide fragmentation logic

### Task 3.1: Create Peptide Module Structure
- [ ] Create `glycofrag/peptide/` directory
- [ ] Create `glycofrag/peptide/__init__.py`

### Task 3.2: Extract Peptide Class
- [ ] Create `peptide/peptide.py`
- [ ] Extract Peptide class
- [ ] Extract b/y/c/z ion generation
- [ ] Add docstrings

### Task 3.3: Extract Structure Predictors (Deferred)
- [ ] Create `glycan/structure.py` (for future complex fragmentation)
- [ ] Extract `GlycanStructurePredictor` base class
- [ ] Extract `OGlycanStructurePredictor` class
- [ ] Extract helper methods
- [ ] Add docstrings

**Source Lines in GlypPRM_v01.py:**
- Lines 3153-3500: Structure predictor classes

### Task 3.4: Extract Fragmentation Logic (Deferred)
- [ ] Create `glycan/fragmentation.py` (for future)
- [ ] Extract `generate_fragments()` method
- [ ] Extract B-ion generation logic
- [ ] Extract Y-ion generation logic
- [ ] Extract BY-ion generation logic
- [ ] Extract YY-ion generation logic
- [ ] Extract `_generate_core_by_ions()` method
- [ ] Extract `_generate_core_fucose_yy_ions()` method
- [ ] Extract `_generate_branch_specific_y_ions()` method
- [ ] Extract `_generate_oglycan_yy_ions()` method
- [ ] Extract helper methods (`_format_fragment`, `_format_string`, `_count_residues`)
- [ ] Add docstrings and type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 1950-3100: Fragment generation methods

### Task 3.5: Peptide Module Tests
- [ ] Create `tests/test_glycan.py`
- [ ] Test N-glycan structure prediction
- [ ] Test O-glycan structure prediction
- [ ] Test fragmentation for various glycan codes
- [ ] Test fragment mass accuracy

---

## Phase 3: Peptide Module - Fragmentation
**Target:** Extract peptide fragmentation logic

### Task 3.1: Create Peptide Module Structure
- [ ] Create `glycofrag/peptide/` directory
- [ ] Create `glycofrag/peptide/__init__.py`

### Task 3.2: Extract Peptide Class (`peptide/peptide.py`)
- [ ] Create `peptide/peptide.py`
- [ ] Create `Peptide` class
- [ ] Implement `calculate_mass()` method
- [ ] Implement `get_sequence()` method
- [ ] Implement `apply_modifications()` method
- [ ] Add docstrings and type hints

### Task 3.3: Extract Peptide Fragmentation (`peptide/fragmentation.py`)
- [ ] Create `peptide/fragmentation.py`
- [ ] Extract `generate_peptide_fragments()` function
- [ ] Extract b-ion generation logic
- [ ] Extract y-ion generation logic
- [ ] Extract Y0/Y1/Y1F generation (glycan stubs)
- [ ] Extract c/z ion generation (optional)
- [ ] Add docstrings and type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 4800-5200: Peptide fragmentation

### Task 3.4: Peptide Module Tests
- [ ] Create `tests/test_peptide.py`
- [ ] Test peptide mass calculation
- [ ] Test b-ion generation
- [ ] Test y-ion generation
- [ ] Test modifications on peptide fragments

---

## Phase 4: Glycopeptide Module - Combined Analysis
**Target:** Combine glycan and peptide for glycopeptide analysis

### Task 4.1: Create Glycopeptide Module Structure
- [ ] Create `glycofrag/glycopeptide/` directory
- [ ] Create `glycofrag/glycopeptide/__init__.py`

### Task 4.2: Create Glycopeptide Class (`glycopeptide/glycopeptide.py`)
- [ ] Create `glycopeptide/glycopeptide.py`
- [ ] Create `Glycopeptide` class
- [ ] Implement `__init__()` with peptide + glycan
- [ ] Implement `calculate_mass()` method
- [ ] Implement `get_peptide()` method
- [ ] Implement `get_glycan()` method
- [ ] Add docstrings and type hints

### Task 4.3: Combined Fragmentation (`glycopeptide/fragmentation.py`)
- [ ] Create `glycopeptide/fragmentation.py`
- [ ] Implement combined fragment generation
- [ ] Handle peptide+glycan Y-ions
- [ ] Handle BYY ions (optional)
- [ ] Handle charge state generation
- [ ] Add docstrings and type hints

### Task 4.4: Glycopeptide Module Tests
- [ ] Create `tests/test_glycopeptide.py`
- [ ] Test glycopeptide mass calculation
- [ ] Test combined fragmentation
- [ ] Test various glycan types (N/O)

---

## Phase 5: Utilities Module
**Target:** Extract helper functions and utilities

### Task 5.1: Create Utils Module Structure
- [ ] Create `glycofrag/utils/` directory
- [ ] Create `glycofrag/utils/__init__.py`

### Task 5.2: Extract Caching Utilities (`utils/caching.py`)
- [ ] Create `utils/caching.py`
- [ ] Extract `create_comprehensive_peptide_cache_key()` function
- [ ] Extract `create_comprehensive_fragment_cache_key()` function
- [ ] Extract `create_peptide_cache_key()` function
- [ ] Add docstrings

**Source Lines in GlypPRM_v01.py:**
- Lines 1042-1147: Cache key functions

### Task 5.3: Create Validation Utilities (`utils/validation.py`)
- [ ] Create `utils/validation.py`
- [ ] Create `validate_glycan_code()` function
- [ ] Create `validate_peptide_sequence()` function
- [ ] Create `validate_modifications()` function
- [ ] Add docstrings

### Task 5.4: Utils Module Tests
- [ ] Test cache key generation
- [ ] Test validation functions

---

## Phase 6: I/O Module - Tables & Export
**Target:** Fragment table generation and export

### Task 6.1: Create I/O Module Structure
- [ ] Create `glycofrag/io/` directory
- [ ] Create `glycofrag/io/__init__.py`

### Task 6.2: Extract Table Generation (`io/tables.py`)
- [ ] Create `io/tables.py`
- [ ] Extract `generate_fragment_table()` function
- [ ] Extract `generate_extended_fragment_table()` function
- [ ] Extract `add_custom_fragments()` function
- [ ] Extract `filter_fragments_by_generation_settings()` function
- [ ] Extract `deduplicate_fragments_by_mass()` function
- [ ] Add DataFrame column documentation
- [ ] Add docstrings and type hints

**Source Lines in GlypPRM_v01.py:**
- Lines 4700-4800: Table generation functions
- Lines 3800-4700: Fragment addition and filtering

### Task 6.3: I/O Module Tests
- [ ] Test fragment table generation
- [ ] Test DataFrame output format
- [ ] Test deduplication

---

## Phase 7: Package Configuration & Documentation
**Target:** Make package installable and documented

### Task 7.1: Create Package Configuration
- [ ] Create `pyproject.toml`
- [ ] Define package metadata
- [ ] Define dependencies (numpy, pandas, networkx, scipy)
- [ ] Define optional dependencies
- [ ] Define entry points (if CLI needed)

### Task 7.2: Create Documentation
- [ ] Create `README.md` with:
  - [ ] Installation instructions
  - [ ] Quick start guide
  - [ ] API overview
  - [ ] Example usage
- [ ] Create `CHANGELOG.md`
- [ ] Create `CONTRIBUTING.md` (optional)
- [ ] Add docstrings to all public APIs

### Task 7.3: Final Testing
- [ ] Run all unit tests
- [ ] Test installation via pip
- [ ] Test import and basic usage
- [ ] Validate mass calculations against known values

### Task 7.4: Package Release (Optional)
- [ ] Register on PyPI (test.pypi.org first)
- [ ] Create GitHub repository
- [ ] Set up CI/CD (GitHub Actions)
- [ ] Publish to PyPI

---

## Module Specifications

### `core/constants.py` Specification

```python
"""
Fundamental mass constants for glycan and peptide analysis.

This module contains all mass constants used throughout the glycofrag package,
including monosaccharide masses, amino acid masses, modification masses,
and physical constants.
"""

# Monosaccharide monoisotopic masses (Da)
MONOSACCHARIDE_MASSES: dict[str, float]

# Hydroxyl groups per monosaccharide (for permethylation)
MONOSACCHARIDE_OH_GROUPS: dict[str, int]

# Modification masses (Da)
MODIFICATION_MASSES: dict[str, float]

# Reducing end masses by type
REDUCING_END_MASSES: dict[int, float]

# Standard amino acid masses (Da)
AMINO_ACID_MASSES: dict[str, float]

# Modification target amino acids
MODIFICATION_TARGETS: dict[str, list[str]]

# Physical constants
PROTON_MASS: float  # 1.00727646688 Da
WATER_MASS: float   # 18.010564684 Da
```

### `core/mass_calculator.py` Specification

```python
"""
Mass calculator for glycans, peptides, and glycopeptides.

The GlycanMassCalculator class provides methods to calculate:
- Glycan masses from composition codes
- Peptide masses with modifications
- Glycopeptide masses
- Fragment ion masses
- m/z values for various charge states
"""

class GlycanMassCalculator:
    """
    Calculator for glycan and glycopeptide masses.
    
    Attributes:
        modification_type (int): Type of reducing end modification (0-6)
        use_cam (bool): Whether to apply carbamidomethylation to cysteines
        fixed_mods (list): List of fixed modifications
        variable_mods (list): List of variable modifications
        mod_string (str): Modification string from external source
        peptide (str): Peptide sequence
    
    Example:
        >>> calc = GlycanMassCalculator()
        >>> mass = calc.calculate_glycan_mass("4501")
        >>> print(f"Glycan mass: {mass:.4f} Da")
    """
    
    def __init__(self, modification_type=6, use_cam=True, 
                 fixed_mods=None, variable_mods=None, 
                 mod_string=None, peptide=None):
        """Initialize the mass calculator."""
        
    def parse_glycan_code(self, code: str) -> tuple[int, int, int, int, int]:
        """
        Parse a glycan composition code.
        
        Args:
            code: Glycan code in format "4501" or "HexNAc(4)Hex(5)Fuc(0)NeuAc(1)"
            
        Returns:
            Tuple of (HexNAc, Hex, Fuc, NeuAc, NeuGc) counts
        """
        
    def calculate_glycan_mass(self, glycan_code: str) -> float:
        """
        Calculate the monoisotopic mass of a glycan.
        
        Args:
            glycan_code: Glycan composition code
            
        Returns:
            Monoisotopic mass in Daltons
        """
        
    def calculate_peptide_mass(self, peptide: str, 
                               use_cam: bool = None,
                               fixed_mods: list = None,
                               variable_mods: list = None,
                               mod_string: str = None) -> float:
        """
        Calculate the monoisotopic mass of a peptide with modifications.
        
        Args:
            peptide: Amino acid sequence (single letter codes)
            use_cam: Apply carbamidomethylation to cysteines
            fixed_mods: List of fixed modifications (e.g., ["CAM:C"])
            variable_mods: List of variable modifications (e.g., ["Ox:M"])
            mod_string: Modification string (e.g., "Ox:M6")
            
        Returns:
            Monoisotopic mass in Daltons
        """
        
    def calculate_mz(self, mass: float, charge: int) -> float:
        """
        Calculate m/z for a given mass and charge state.
        
        Args:
            mass: Neutral mass in Daltons
            charge: Charge state (positive integer)
            
        Returns:
            m/z value
        """
        
    def calculate_fragment_mass(self, composition: dict, 
                                fragment_type: str,
                                peptide: str = None) -> float:
        """
        Calculate the mass of a glycan/glycopeptide fragment.
        
        Args:
            composition: Dict with monosaccharide counts
            fragment_type: Type of fragment ('b_ions', 'y_ions', etc.)
            peptide: Peptide sequence (for glycopeptide fragments)
            
        Returns:
            Fragment mass in Daltons
        """
```

### `glycan/glycan.py` Specification

```python
"""
Glycan structure representation and analysis.

The Glycan class represents a glycan by its monosaccharide composition
and provides methods for structure prediction and fragmentation.
"""

class Glycan:
    """
    Represents a glycan structure with composition and fragmentation capabilities.
    
    Attributes:
        glycan_type (str): 'N' for N-glycan, 'O' for O-glycan
        hexnac_total (int): Total HexNAc count
        hex_total (int): Total Hex count
        fuc_total (int): Total Fuc count
        neuac_total (int): Total NeuAc count
        neugc_total (int): Total NeuGc count
        possible_structures (list): List of predicted graph structures
    
    Example:
        >>> glycan = Glycan("4501", glycan_type="N")
        >>> structures = glycan.predict_structures()
        >>> print(f"Found {len(structures)} possible structures")
    """
    
    def __init__(self, code: str, glycan_type: str = "N", max_structures: int = 100):
        """
        Initialize a Glycan from a composition code.
        
        Args:
            code: Glycan composition code
            glycan_type: 'N' for N-glycan, 'O' for O-glycan
            max_structures: Maximum structures to predict
        """
        
    def predict_structures(self) -> list:
        """
        Predict possible glycan structures based on composition.
        
        Returns:
            List of networkx DiGraph objects representing structures
        """
        
    def generate_fragments(self, graph, modification_type: int = 6,
                          peptide: str = None, **kwargs) -> tuple[dict, dict]:
        """
        Generate fragment ions for a glycan structure.
        
        Args:
            graph: NetworkX DiGraph of glycan structure
            modification_type: Reducing end modification type
            peptide: Peptide sequence (for glycopeptide mode)
            
        Returns:
            Tuple of (fragments dict, cleavage_info dict)
            fragments contains: b_ions, y_ions, by_ions, yy_ions
        """
        
    def classify_structure(self, graph) -> str:
        """
        Classify glycan structure type.
        
        Args:
            graph: NetworkX DiGraph of glycan structure
            
        Returns:
            Classification string (e.g., "Complex", "High Mannose", "Hybrid")
        """
```

### `peptide/fragmentation.py` Specification

```python
"""
Peptide fragmentation for mass spectrometry analysis.

This module provides functions to generate peptide fragment ions
(b, y, c, z series) with support for modifications.
"""

def generate_peptide_fragments(peptide_sequence: str,
                               glycan_code: str = None,
                               calculator = None,
                               modification_type: int = 6,
                               use_cam: bool = False,
                               fixed_mods: list = None,
                               variable_mods: list = None,
                               mod_string: str = None) -> dict:
    """
    Generate peptide b and y ion fragments.
    
    Args:
        peptide_sequence: Amino acid sequence
        glycan_code: Glycan code (for Y0/Y1/Y1F generation)
        calculator: GlycanMassCalculator instance
        modification_type: Modification type
        use_cam: Apply carbamidomethylation
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Modification string
        
    Returns:
        Dict with 'b_ions' and 'y_ions' lists
    """
```

### `io/tables.py` Specification

```python
"""
Fragment table generation and export.

This module provides functions to generate pandas DataFrames
containing fragment information for mass spectrometry analysis.
"""

def generate_fragment_table(fragments: dict,
                           glycan_code: str,
                           modification_type: int = 6,
                           peptide: str = None,
                           use_cam: bool = False,
                           fixed_mods: list = None,
                           variable_mods: list = None,
                           mod_string: str = None,
                           glycan_type: str = None) -> pd.DataFrame:
    """
    Generate a DataFrame of fragment ions with masses and metadata.
    
    Args:
        fragments: Dict of fragment lists by type
        glycan_code: Glycan composition code
        modification_type: Reducing end modification type
        peptide: Peptide sequence
        use_cam: Apply carbamidomethylation
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Modification string
        glycan_type: 'N' or 'O' glycan type
        
    Returns:
        DataFrame with columns:
        - FragmentName: Fragment identifier
        - FragmentType: Ion type (B, Y, BY, YY, etc.)
        - Composition: Monosaccharide composition
        - Mass: Neutral mass
        - Charge: Charge state
        - m/z: Mass-to-charge ratio
    """
```

---

## Progress Tracking

| Phase | Task | Status | Date Completed |
|-------|------|--------|----------------|
| 1 | Package Structure | ✅ Completed | 2026-02-01 |
| 1 | Constants Module | ✅ Completed | 2026-02-01 |
| 1 | Mass Calculator | ✅ Completed | 2026-02-01 |
| 1 | Modifications | ✅ Completed | 2026-02-01 |
| 1 | Core Tests | ✅ Completed | 2026-02-01 |
| 2 | Glycan Class | ⬜ Not Started | |
| 2 | Structure Predictors | ⬜ Not Started | |
| 2 | Fragmentation | ⬜ Not Started | |
| 2 | Glycan Tests | ⬜ Not Started | |
| 3 | Peptide Class | ⬜ Not Started | |
| 3 | Peptide Fragmentation | ⬜ Not Started | |
| 3 | Peptide Tests | ⬜ Not Started | |
| 4 | Glycopeptide Class | ⬜ Not Started | |
| 4 | Combined Fragmentation | ⬜ Not Started | |
| 4 | Glycopeptide Tests | ⬜ Not Started | |
| 5 | Caching Utils | ✅ Completed | 2026-02-01 |
| 5 | Validation Utils | ✅ Completed | 2026-02-01 |
| 6 | Table Generation | ⬜ Not Started | |
| 7 | Package Config | ✅ Completed | 2026-02-01 |
| 7 | Documentation | 🟡 In Progress | |
| 7 | Release | ⬜ Not Started | |

**Legend:**
- ⬜ Not Started
- 🟡 In Progress
- ✅ Completed
- ❌ Blocked

---

## Dependencies

### Required
- `numpy>=1.20.0` - Numerical operations
- `pandas>=1.3.0` - DataFrame operations
- `networkx>=2.6.0` - Graph structures for glycan representation
- `scipy>=1.7.0` - Scientific computing utilities

### Development
- `pytest>=7.0.0` - Testing framework
- `pytest-cov>=3.0.0` - Coverage reporting
- `black>=22.0.0` - Code formatting
- `mypy>=0.950` - Type checking
- `sphinx>=4.0.0` - Documentation generation

---

## Notes & Decisions

### Design Decisions
1. **NetworkX for structures**: Continue using networkx DiGraph for glycan structures as it provides good graph manipulation capabilities
2. **Caching strategy**: Keep LRU caching for performance-critical calculations
3. **Modification handling**: Support both string-based (e.g., "Ox:M6") and dict-based modifications
4. **Type hints**: Add type hints to all public APIs for better IDE support

### Known Issues to Address
1. Remove all `print()` debug statements during extraction
2. Replace GUI-specific error handling with proper exceptions
3. Ensure thread-safety for potential multi-threaded use

### Future Enhancements
1. Add CLI interface for command-line usage
2. Add visualization functions for glycan structures
3. Add mzML/mzXML file reading capabilities (optional module)
4. Add spectral matching capabilities

---

*Last Updated: February 1, 2026*
