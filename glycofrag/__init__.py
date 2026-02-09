"""
GlycoFrag - A Python package for glycan and glycopeptide fragmentation analysis.

This package provides tools for:
- Glycan mass calculation
- Peptide mass calculation with modifications
- Glycan structure prediction (N-glycans and O-glycans)
- Fragment ion generation (B, Y, BY, YY ions)
- Glycopeptide fragmentation analysis
- MS/MS structure prediction from experimental data

Quick Start (MS/MS Analysis - Complete Pipeline):
    >>> from glycofrag import GlycoPeptideAnalysis
    >>> 
    >>> # All-in-one initialization (fragments generated automatically)
    >>> analysis = GlycoPeptideAnalysis(
    ...     peptide_sequence='EEQYNSTYR',
    ...     glycan_code='4501',
    ...     glycosylation_site=5,
    ...     glycan_type='N',
    ...     peptide_fragment_types=['by', 'cz'],
    ...     glycan_fragment_types=['BY'],
    ...     visualize_structure='best'
    ... )
    >>>
    >>> # Access generated data
    >>> print(analysis.theoretical_df)  # 45 fragments
    >>> print(analysis.clean_theoretical_df)  # 43 deduplicated
    >>> print(analysis.summary_df)  # masses & mods
    >>>
    >>> # Match experimental peaks
    >>> matched = analysis.match_fragments(spectrum, ppm_tolerance=20)
    >>>
    >>> # Rank structures
    >>> prediction = analysis.predict_structures(matched)
    >>> analysis.visualize_best_structure(output_image='structure.png')
    >>> analysis.export_to_excel('analysis.xlsx', matched)

Basic Usage (Fragment Generation):
    >>> from glycofrag import Glycan, Glycopeptide, GlycanVisualizer
    >>>
    >>> # Released glycan fragments
    >>> glycan = Glycan("4501", glycan_type="N")
    >>> structures = glycan.predict_structures()
    >>> fragments, info = glycan.generate_fragments(structures[0])
    >>>
    >>> # Glycopeptide fragments
    >>> gp = Glycopeptide("EEQYNSTYR", "4501", 5, "N")
    >>> frags = gp.generate_fragments()
    >>>
    >>> # Visualize structures
    >>> GlycanVisualizer.draw(gp.glycan_structures, 'all', glycan_code='4501')

Author: Oluwatosin Daramola
License: MIT
"""

__version__ = "0.1.0"
__author__ = "Oluwatosin Daramola"

# Core imports
from glycofrag.core.constants import (
    MONOSACCHARIDE_MASSES,
    AMINO_ACID_MASSES,
    MODIFICATION_MASSES,
    MODIFICATION_TARGETS,
    PROTON_MASS,
    WATER_MASS,
)

from glycofrag.core.modifications import ReducingEndType, get_modification_type

from glycofrag.core.mass_calculator import GlycanMassCalculator

# New clean import - Glycan facade with unified API
from glycofrag.glycan.facade import Glycan

# Peptide module
from glycofrag.peptide import Peptide

# Glycopeptide module
from glycofrag.glycopeptide import Glycopeptide

# High-level analysis APIs
from glycofrag.analysis import GlycoPeptideAnalysis
from glycofrag.glycan_analysis import GlycanAnalysis

# Visualization
from glycofrag.io.visualizer import GlycanVisualizer

# I/O utilities
from glycofrag.io.tables import (
    format_fragment_string,
    extract_fragment_composition,
    generate_fragment_table,
    generate_all_fragments_table,
    deduplicate_fragments_by_mass,
    export_fragment_table_to_excel
)

__all__ = [
    # Constants
    "MONOSACCHARIDE_MASSES",
    "AMINO_ACID_MASSES",
    "MODIFICATION_MASSES",
    "MODIFICATION_TARGETS",
    "PROTON_MASS",
    "WATER_MASS",
    # Core classes
    "GlycanMassCalculator",
    # Glycan classes
    "Glycan",
    "ReducingEndType",
    "get_modification_type",
    # Peptide classes
    "Peptide",
    # Glycopeptide classes
    "Glycopeptide",
    # High-level analysis APIs (recommended for MS/MS workflows)
    "GlycoPeptideAnalysis",  # For glycopeptide MS/MS analysis
    "GlycanAnalysis",         # For released glycan MS/MS analysis
    # Visualization
    "GlycanVisualizer",
    # I/O utilities
    "format_fragment_string",
    "extract_fragment_composition",
    "generate_fragment_table",
    "generate_all_fragments_table",
    "deduplicate_fragments_by_mass",
    "export_fragment_table_to_excel",
    # Utility functions
    "list_supported_modifications",
]


def list_supported_modifications(verbose: bool = True) -> dict:
    """
    List all supported peptide modifications with their masses and target residues.
    
    Args:
        verbose: If True, prints formatted output. If False, returns dict only.
    
    Returns:
        Dictionary with modification info: {mod_name: {'mass': float, 'targets': list}}
    
    Example:
        >>> from glycofrag import list_supported_modifications
        >>> mods = list_supported_modifications(verbose=True)
        >>> # Use specific modification
        >>> mass = mods['Ox']['mass']
        >>> targets = mods['Ox']['targets']
    """
    result = {}
    for mod_name in MODIFICATION_MASSES:
        result[mod_name] = {
            'mass': MODIFICATION_MASSES[mod_name],
            'targets': MODIFICATION_TARGETS.get(mod_name, [])
        }
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"GLYCOFRAG SUPPORTED MODIFICATIONS ({len(MODIFICATION_MASSES)} total)")
        print(f"{'='*70}")
        
        # Group modifications
        fixed = ['CAM', 'PAM', 'Palm', 'Carbamyl', 'TMT6', 'TMT10', 'TMT16', 'iTRAQ4', 'iTRAQ8']
        variable = [m for m in MODIFICATION_MASSES if m not in fixed]
        
        print(f"\n[FIXED MODIFICATIONS] ({len(fixed)}):")
        print("-" * 70)
        for mod in fixed:
            if mod in result:
                mass = result[mod]['mass']
                targets_str = ', '.join(result[mod]['targets'])
                print(f"  {mod:15s} {mass:10.6f} Da  ->  {targets_str}")
        
        print(f"\n[VARIABLE MODIFICATIONS] ({len(variable)}):")
        print("-" * 70)
        for mod in sorted(variable):
            mass = result[mod]['mass']
            targets_str = ', '.join(result[mod]['targets'])
            print(f"  {mod:15s} {mass:10.6f} Da  ->  {targets_str}")
        
        print(f"\n{'='*70}")
        print("USAGE:")
        print(f"{'='*70}")
        print('  peptide = Peptide("PEPTIDE", mod_string="M:Ox")')
        print('  peptide = Peptide("PEPTIDE", mod_string="M4,24:Ox")')
        print('  peptide = Peptide("PEPTIDE", mod_string="M:Ox; S3:Phos")')
        print('  gp = Glycopeptide("PEPTIDE", "4501", 5, "N", mod_string="M:Ox")')
        print(f"{'='*70}\n")
    
    return result
