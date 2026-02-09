"""
GlycoFrag Core Module.

This module contains fundamental components:
- constants: Mass constants and reference data
- mass_calculator: Mass calculation engine
- modifications: Modification parsing and handling
"""

from glycofrag.core.constants import (
    MONOSACCHARIDE_MASSES,
    MONOSACCHARIDE_OH_GROUPS,
    MODIFICATION_MASSES,
    MODIFICATION_TARGETS,
    REDUCING_END_MASSES,
    ADDITIONAL_MODIFICATIONS,
    AMINO_ACID_MASSES,
    PROTON_MASS,
    WATER_MASS,
)

from glycofrag.core.mass_calculator import GlycanMassCalculator

from glycofrag.core.modifications import (
    parse_modification_string,
    parse_custom_mass_mod,
)

__all__ = [
    # Constants
    "MONOSACCHARIDE_MASSES",
    "MONOSACCHARIDE_OH_GROUPS",
    "MODIFICATION_MASSES",
    "MODIFICATION_TARGETS",
    "REDUCING_END_MASSES",
    "ADDITIONAL_MODIFICATIONS",
    "AMINO_ACID_MASSES",
    "PROTON_MASS",
    "WATER_MASS",
    # Classes
    "GlycanMassCalculator",
    # Functions
    "parse_modification_string",
    "parse_custom_mass_mod",
]
