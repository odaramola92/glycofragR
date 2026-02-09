"""
Fundamental mass constants for glycan and peptide analysis.

This module contains all mass constants used throughout the glycofrag package,
including monosaccharide masses, amino acid masses, modification masses,
and physical constants.

All masses are monoisotopic values in Daltons (Da).

References:
    - Unimod: https://www.unimod.org/
    - GlycoMod: https://web.expasy.org/glycomod/
"""

from typing import Dict, List


# =============================================================================
# Physical Constants
# =============================================================================

#: Proton mass in Daltons (high precision for m/z calculations)
PROTON_MASS: float = 1.00727646688

#: Water mass in Daltons (H2O, monoisotopic)
WATER_MASS: float = 18.010564684

#: Hydrogen mass in Daltons
HYDROGEN_MASS: float = 1.00782503207

# =============================================================================
# Monosaccharide Masses
# =============================================================================

#: Monoisotopic masses of common monosaccharides in Daltons
#: These are the residue masses (mass of monosaccharide - H2O)
#: Permethylation will be added based on modification type
MONOSACCHARIDE_MASSES: Dict[str, float] = {
    'HexNAc': 203.079373,    # N-acetylhexosamine (GlcNAc, GalNAc)
    'Hex': 162.052824,       # Hexose (Glucose, Galactose, Mannose)
    'Fuc': 146.057909,       # Fucose (deoxyhexose)
    'NeuAc': 291.095417,     # N-acetylneuraminic acid (sialic acid)
    'NeuGc': 307.090331,     # N-glycolylneuraminic acid
}

#: Number of hydroxyl (-OH) groups per monosaccharide
#: Used for calculating permethylation mass shifts
#: These are BASE values - B-ions add one extra methylation for terminal protection
#:   HexNAc: 3 OH → 245.1265 Da (base and B-ion)
#:   Hex: 3 OH → 204.0999 Da (base); 4 OH → 218.1156 Da (B-ion with extra methylation)
#:   Fuc: 2 OH → 188.1050 Da (base) 
#:   NeuAc: 5 OH → 361.1739 Da (base); 6 OH → 375.1896 Da (B-ion with extra methylation)
#:   NeuGc: 6 OH (base)
MONOSACCHARIDE_OH_GROUPS: Dict[str, int] = {
    'HexNAc': 3,
    'Hex': 3,
    'Fuc': 2,
    'NeuAc': 5,
    'NeuGc': 6,
}

# =============================================================================
# Modification Masses
# =============================================================================

#: Mass shift for permethylation per hydroxyl group (CH3 group mass)
PERMETHYLATION_MASS: float = 14.0157

#: Mass of methyl group (CH3)
METHYL_MASS: float = 14.0157

#: Reducing end masses for different modification types
#: Key: modification_type (int), Value: mass (float)
REDUCING_END_MASSES: Dict[int, float] = {
    0: 18.010564684,   # Free reducing end (H2O)
    1: 20.0262,        # Reduced end (alditol)
    2: 18.010564684,   # Permethylated free end (before permethylation adjustment)
    3: 20.0262,        # Permethylated reduced end
    4: 18.010564684,   # 2-AB labeled (before 2-AB addition)
    5: 18.010564684,   # 2-AB labeled and permethylated
    6: 0.0,            # Glycopeptide mode - peptide replaces reducing end
    7: 0.0,            # Custom mode - user supplies reducing end mass via custom_reducing_end_mass
}

#: Additional mass modifications for specific reducing end types
ADDITIONAL_MODIFICATIONS: Dict[int, float] = {
    0: 0.0,            # Free end - no additional modification
    1: 0.0,            # Reduced end - no additional modification
    2: 28.0313,        # Permethylation mass adjustment
    3: 42.0471,        # Permethylation mass adjustment (reduced end)
    4: 120.0688,       # 2-AB label mass (2-aminobenzamide - H2O)
    5: 190.1451,       # 2-AB + permethylation combined
    6: 18.010564684,   # Glycopeptide - water mass for glycosidic bond
    7: 0.0,            # Custom mode - user supplies additional mass via custom_reducing_end_mass
}

# =============================================================================
# Amino Acid Masses
# =============================================================================

#: Monoisotopic residue masses for standard amino acids in Daltons
#: These are residue masses (amino acid - H2O)
AMINO_ACID_MASSES: Dict[str, float] = {
    'A': 71.037114,    # Alanine
    'C': 103.009185,   # Cysteine
    'D': 115.026943,   # Aspartic acid
    'E': 129.042593,   # Glutamic acid
    'F': 147.068414,   # Phenylalanine
    'G': 57.021464,    # Glycine
    'H': 137.058912,   # Histidine
    'I': 113.084064,   # Isoleucine
    'K': 128.094963,   # Lysine
    'L': 113.084064,   # Leucine
    'M': 131.040485,   # Methionine
    'N': 114.042927,   # Asparagine
    'P': 97.052764,    # Proline
    'Q': 128.058578,   # Glutamine
    'R': 156.101111,   # Arginine
    'S': 87.032028,    # Serine
    'T': 101.047679,   # Threonine
    'V': 99.068414,    # Valine
    'W': 186.079313,   # Tryptophan
    'Y': 163.063329,   # Tyrosine
}

# =============================================================================
# Peptide Modification Masses
# =============================================================================

#: Mass changes for common peptide modifications in Daltons
MODIFICATION_MASSES: Dict[str, float] = {
    # Fixed modifications (typically applied to all occurrences)
    'CAM': 57.021464,       # Carbamidomethylation (Cys) - iodoacetamide
    'PAM': 71.03711,        # Propionamide (Cys) - acrylamide
    'Palm': 238.2297,       # Palmitoylation (Cys)
    'Carbamyl': 43.0058,    # Carbamylation (N-term, Lys)
    
    # Isobaric labeling reagents
    'TMT6': 229.16293,      # TMT 6-plex (Lys, N-term)
    'TMT10': 229.16293,     # TMT 10-plex (Lys, N-term)
    'TMT16': 304.20710,     # TMTpro 16-plex (Lys, N-term)
    'iTRAQ4': 144.10206,    # iTRAQ 4-plex (Lys, N-term)
    'iTRAQ8': 304.20536,    # iTRAQ 8-plex (Lys, N-term)
    
    # Variable modifications
    'Ox': 15.994915,        # Oxidation (Met)
    'Deam': 0.984016,       # Deamidation (Asn, Gln)
    'Phos': 79.966331,      # Phosphorylation (Ser, Thr, Tyr)
    'Ac': 42.010565,        # Acetylation (Lys, N-term)
    'Methyl': 14.015650,    # Methylation (Lys, Arg)
    'DiMethyl': 28.031300,  # Dimethylation (Lys, Arg)
    'TriMethyl': 42.046950, # Trimethylation (Lys)
    'Pyro-glu': -17.026549, # Pyroglutamic acid (N-term Gln)
    'Pyro-cmC': -17.026549, # Pyroglutamic acid (N-term CAM-Cys)
    'GG': 114.042927,       # GlyGly (Lys) - ubiquitination remnant
    'HexNAc': 203.079373,   # O-GlcNAc (Ser, Thr)
    'Formyl': 27.994915,    # Formylation (Lys, N-term)
    'Nitration': 44.985078, # Nitration (Tyr)
    'Sulf': 79.956815,      # Sulfation (Tyr)
    'Biotin': 226.0776,     # Biotinylation (Lys)
    'Malonyl': 86.000394,   # Malonylation (Lys)
    'Succinyl': 100.016044, # Succinylation (Lys)
    'Myristoyl': 210.198366,# Myristoylation (N-term Gly)
    'Farnesyl': 204.187801, # Farnesylation (Cys)
    'SUMO1-GG': 213.14430,  # SUMOylation remnant (Lys)
}

#: Target amino acids for each modification type
MODIFICATION_TARGETS: Dict[str, List[str]] = {
    # Fixed modifications
    'CAM': ['C'],                         # Carbamidomethylation on Cys
    'PAM': ['C'],                         # Propionamide on Cys
    'Palm': ['C'],                        # Palmitoylation on Cys
    'Carbamyl': ['K', 'N-term'],          # Carbamylation on Lys, N-terminus
    'TMT6': ['K', 'N-term'],              # TMT 6-plex
    'TMT10': ['K', 'N-term'],             # TMT 10-plex
    'TMT16': ['K', 'N-term'],             # TMTpro 16-plex
    'iTRAQ4': ['K', 'N-term'],            # iTRAQ 4-plex
    'iTRAQ8': ['K', 'N-term'],            # iTRAQ 8-plex
    
    # Variable modifications
    'Ox': ['M'],                          # Oxidation on Met
    'Deam': ['N', 'Q'],                   # Deamidation on Asn, Gln
    'Phos': ['S', 'T', 'Y'],              # Phosphorylation on Ser, Thr, Tyr
    'Ac': ['K', 'N-term'],                # Acetylation on Lys, N-terminus
    'Methyl': ['K', 'R'],                 # Methylation on Lys, Arg
    'DiMethyl': ['K', 'R'],               # Dimethylation on Lys, Arg
    'TriMethyl': ['K'],                   # Trimethylation on Lys
    'Pyro-glu': ['N-term-Q'],             # Pyroglutamic from N-term Gln
    'Pyro-cmC': ['N-term-C-CAM'],         # Pyroglutamic from N-term CAM-Cys
    'GG': ['K'],                          # GlyGly (ubiquitin) on Lys
    'HexNAc': ['S', 'T'],                 # O-GlcNAc on Ser, Thr
    'Formyl': ['K', 'N-term'],            # Formylation on Lys, N-terminus
    'Nitration': ['Y'],                   # Nitration on Tyr
    'Sulf': ['Y'],                        # Sulfation on Tyr
    'Biotin': ['K'],                      # Biotinylation on Lys
    'Malonyl': ['K'],                     # Malonylation on Lys
    'Succinyl': ['K'],                    # Succinylation on Lys
    'Myristoyl': ['N-term-G'],            # Myristoylation on N-term Gly
    'Farnesyl': ['C'],                    # Farnesylation on Cys
    'SUMO1-GG': ['K'],                    # SUMOylation remnant on Lys
}

# =============================================================================
# Glycan Fragment Constants
# =============================================================================

#: Common glycan fragment neutral losses
GLYCAN_NEUTRAL_LOSSES: Dict[str, float] = {
    'H2O': 18.010565,       # Water loss
    'HexNAc': 203.079373,   # Loss of HexNAc
    'Hex': 162.052824,      # Loss of Hex
    'Fuc': 146.057909,      # Loss of Fuc
    'NeuAc': 291.095417,    # Loss of NeuAc
}

#: Oxonium ion masses (common glycan diagnostic ions)
OXONIUM_IONS: Dict[str, float] = {
    'HexNAc': 204.086646,           # [HexNAc + H]+
    'HexNAc-H2O': 186.076081,       # [HexNAc - H2O + H]+
    'HexNAc-2H2O': 168.065516,      # [HexNAc - 2H2O + H]+
    'Hex': 163.060097,              # [Hex + H]+
    'Fuc': 147.065182,              # [Fuc + H]+
    'NeuAc': 292.102690,            # [NeuAc + H]+
    'NeuAc-H2O': 274.092126,        # [NeuAc - H2O + H]+
    'HexHexNAc': 366.139472,        # [Hex-HexNAc + H]+
    'HexHexNAc-H2O': 348.128907,    # [Hex-HexNAc - H2O + H]+
}
