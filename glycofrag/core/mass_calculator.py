"""
Mass calculator for glycans, peptides, and glycopeptides.

The GlycanMassCalculator class provides methods to calculate:
- Glycan masses from composition codes
- Peptide masses with modifications
- Glycopeptide masses
- Fragment ion masses
- m/z values for various charge states
"""

import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Union

from glycofrag.core.constants import (
    MONOSACCHARIDE_MASSES,
    MONOSACCHARIDE_OH_GROUPS,
    AMINO_ACID_MASSES,
    MODIFICATION_MASSES,
    MODIFICATION_TARGETS,
    REDUCING_END_MASSES,
    ADDITIONAL_MODIFICATIONS,
    PROTON_MASS,
    WATER_MASS,
    PERMETHYLATION_MASS,
)
from glycofrag.core.modifications import (
    parse_modification_string,
    parse_custom_mass_mod,
    apply_modifications_to_positions,
)


class GlycanMassCalculator:
    """
    Calculator for glycan and glycopeptide masses.
    
    This class provides comprehensive mass calculation capabilities for:
    - Glycans (from composition codes like "4501" or "HexNAc(4)Hex(5)...")
    - Peptides (with support for various modifications)
    - Glycopeptides (combined glycan + peptide)
    - Fragment ions (B, Y, BY, YY series)

    Public API:
        - parse_glycan_code(...)
        - calculate_glycan_mass(...)
        - calculate_composition_mass(...)
        - calculate_peptide_mass(...)
        - calculate_glycopeptide_mass(...)
        - calculate_fragment_mass(...)
        - calculate_mz(...)
    
    Attributes:
        modification_type (int): Type of reducing end modification (0-6)
            0: Free reducing end (H2O)
            1: Reduced end (alditol)
            2: Permethylated free end
            3: Permethylated reduced end
            4: 2-AB labeled
            5: 2-AB labeled and permethylated
            6: Glycopeptide mode (peptide at reducing end)
        use_cam (bool): Whether to apply carbamidomethylation to cysteines
        fixed_mods (list): List of fixed modifications (e.g., ["CAM:C"])
        variable_mods (list): List of variable modifications (e.g., ["Ox:M"])
        mod_string (str): Modification string from external source
        peptide (str): Associated peptide sequence
    
    Example:
        >>> calc = GlycanMassCalculator()
        >>> 
        >>> # Calculate glycan mass
        >>> mass = calc.calculate_glycan_mass("4501")
        >>> print(f"Glycan mass: {mass:.4f} Da")
        >>>
        >>> # Calculate peptide mass
        >>> pep_mass = calc.calculate_peptide_mass("LCPDCPLLAPLNDSR", use_cam=True)
        >>> print(f"Peptide mass: {pep_mass:.4f} Da")
        >>>
        >>> # Calculate m/z
        >>> mz = calc.calculate_mz(mass, charge=2)
        >>> print(f"m/z at z=2: {mz:.4f}")
    """
    
    def __init__(
        self,
        modification_type: int = 6,
        use_cam: bool = True,
        fixed_mods: Optional[List[str]] = None,
        variable_mods: Optional[List[str]] = None,
        mod_string: Optional[str] = None,
        peptide: Optional[str] = None,
        custom_reducing_end_mass: Optional[float] = None,
        custom_mono_masses: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize the mass calculator.
        
        Args:
            modification_type: Reducing end modification type (0-7)
            use_cam: Apply carbamidomethylation to cysteines
            fixed_mods: List of fixed modifications
            variable_mods: List of variable modifications
            mod_string: Modification string from external source
            peptide: Associated peptide sequence
            custom_reducing_end_mass: Custom mass (Da) to add at the reducing end.
                Only used when modification_type=7 ('custom').
            custom_mono_masses: Dict of additional mass to add per monosaccharide type.
                Only used when modification_type=7 ('custom').
                Example: {'HexNAc': 14.02, 'Hex': 14.02} adds 14.02 Da to each
                HexNAc and Hex residue.
        """
        # Store configuration
        self.modification_type = modification_type
        self.use_cam = bool(use_cam)
        self.fixed_mods = [str(mod) for mod in (fixed_mods or [])]
        self.variable_mods = [str(mod) for mod in (variable_mods or [])]
        self.mod_string = mod_string
        self.peptide = peptide
        self.custom_reducing_end_mass = custom_reducing_end_mass
        self.custom_mono_masses = custom_mono_masses or {}
        
        # Initialize mass dictionaries (may be modified by permethylation)
        self.MONO_MASSES = MONOSACCHARIDE_MASSES.copy()
        self._calculate_mono_masses()
        
        # Caches for performance
        self._peptide_mass_cache: Dict[str, float] = {}
        self._applied_mods: Dict[str, List[str]] = {}
    
    def _calculate_mono_masses(self) -> None:
        """Calculate monosaccharide masses based on modification type.
        
        For permethylated glycans (modification_type 2, 3, 5):
        - Add permethylation mass for each OH group
        - Reducing end mass is added separately for Y-ions only
        
        For custom modifications (modification_type 7):
        - Add user-specified mass per monosaccharide type
        """
        self.MONO_MASSES = MONOSACCHARIDE_MASSES.copy()
        
        # Apply permethylation if needed
        if self.modification_type in [2, 3, 5]:
            permethyl_mass = 14.0157
            for mono in self.MONO_MASSES:
                oh_groups = MONOSACCHARIDE_OH_GROUPS.get(mono, 0)
                self.MONO_MASSES[mono] += oh_groups * permethyl_mass
        
        # Apply custom monosaccharide mass additions
        if self.modification_type == 7 and self.custom_mono_masses:
            for mono, extra_mass in self.custom_mono_masses.items():
                if mono in self.MONO_MASSES:
                    self.MONO_MASSES[mono] += extra_mass
    
    def parse_glycan_code(self, code: str) -> Tuple[int, int, int, int, int]:
        """
        Parse a glycan composition code.
        
        Supports multiple formats:
        - Numeric: "4501" -> HexNAc(4)Hex(5)Fuc(0)NeuAc(1)
        - Named: "HexNAc(4)Hex(5)Fuc(0)NeuAc(1)"
        
        Args:
            code: Glycan composition code
        
        Returns:
            Tuple of (HexNAc, Hex, Fuc, NeuAc, NeuGc) counts
        
        Examples:
            >>> calc = GlycanMassCalculator()
            >>> calc.parse_glycan_code("4501")
            (4, 5, 0, 1, 0)
            >>> calc.parse_glycan_code("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)")
            (4, 5, 1, 2, 0)
        """
        code = str(code).replace(" ", "")
        
        # Check if named format (starts with letter)
        if code and code[0].isalpha():
            return self._parse_named_glycan_code(code)
        else:
            return self._parse_numeric_glycan_code(code)
    
    def _parse_named_glycan_code(self, code: str) -> Tuple[int, int, int, int, int]:
        """Parse named glycan code format like 'HexNAc(4)Hex(5)...' or 'HexNAc4Hex5...'"""
        mono_counts = {
            'HexNAc': 0,
            'Hex': 0,
            'Fuc': 0,
            'NeuAc': 0,
            'NeuGc': 0,
        }
        
        # Name mappings for various conventions
        name_mappings = {
            'HEXNAC': 'HexNAc',
            'HEXN': 'HexNAc',
            'GLCNAC': 'HexNAc',
            'GALNAC': 'HexNAc',
            'HEX': 'Hex',
            'HEXOSE': 'Hex',
            'GLC': 'Hex',
            'GAL': 'Hex',
            'MAN': 'Hex',
            'FUC': 'Fuc',
            'FUCOSE': 'Fuc',
            'NEUAC': 'NeuAc',
            'SIA': 'NeuAc',
            'SIALIC': 'NeuAc',
            'NEUGC': 'NeuGc',
        }
        
        # Try format with parentheses first: "HexNAc(3)" or "Hex(4)"
        pattern_paren = r'([A-Za-z]+)\((\d+)\)'
        matches_paren = re.findall(pattern_paren, code)
        
        if matches_paren:
            # Format with parentheses found
            for mono_name, count in matches_paren:
                mono_upper = mono_name.upper()
                
                # Map to standard name
                if mono_upper in name_mappings:
                    standard_name = name_mappings[mono_upper]
                    mono_counts[standard_name] += int(count)
                else:
                    # Try prefix matching
                    for known_name, standard in name_mappings.items():
                        if mono_upper.startswith(known_name):
                            mono_counts[standard] += int(count)
                            break
        else:
            # Try format without parentheses: "HexNAc4Hex5Fuc1NeuGc1"
            # Match monosaccharide names followed by digits
            # Order matters: NeuGc before NeuAc, HexNAc before Hex
            pattern_no_paren = r'(NeuGc|NeuAc|HexNAc|Hex|Fuc)(\d+)'
            matches_no_paren = re.findall(pattern_no_paren, code, re.IGNORECASE)
            
            for mono_name, count in matches_no_paren:
                mono_upper = mono_name.upper()
                
                # Map to standard name
                if mono_upper in name_mappings:
                    standard_name = name_mappings[mono_upper]
                    mono_counts[standard_name] += int(count)
        
        return (
            mono_counts['HexNAc'],
            mono_counts['Hex'],
            mono_counts['Fuc'],
            mono_counts['NeuAc'],
            mono_counts['NeuGc'],
        )
    
    def _parse_numeric_glycan_code(self, code: str) -> Tuple[int, int, int, int, int]:
        """Parse numeric glycan code format like '4501'"""
        if code.isdigit() and len(code) >= 4:
            hexnac = int(code[0])
            hex_count = int(code[1])
            fuc = int(code[2])
            neuac = int(code[3])
            neugc = int(code[4]) if len(code) > 4 else 0
            return (hexnac, hex_count, fuc, neuac, neugc)
        else:
            # Fallback pattern matching
            pattern = r'\((\d+)\)|(\d)'
            matches = re.findall(pattern, code)
            components = [int(m[0] or m[1]) for m in matches]
            
            if len(components) >= 5:
                return tuple(components[:5])
            elif len(components) == 4:
                return (components[0], components[1], components[2], components[3], 0)
            else:
                raise ValueError(f"Invalid glycan code format: {code}")
    
    def calculate_glycan_mass(self, glycan_code: str) -> float:
        """
        Calculate the monoisotopic mass of a glycan.
        
        Args:
            glycan_code: Glycan composition code (e.g., "4501")
        
        Returns:
            Monoisotopic mass in Daltons
        
        Example:
            >>> calc = GlycanMassCalculator()
            >>> mass = calc.calculate_glycan_mass("4501")
            >>> print(f"{mass:.4f}")  # ~1622.58 Da
        """
        hexnac, hex_count, fuc, neuac, neugc = self.parse_glycan_code(glycan_code)
        
        mass = 0.0
        mass += hexnac * self.MONO_MASSES['HexNAc']
        mass += hex_count * self.MONO_MASSES['Hex']
        mass += fuc * self.MONO_MASSES['Fuc']
        mass += neuac * self.MONO_MASSES['NeuAc']
        mass += neugc * self.MONO_MASSES['NeuGc']
        
        # Add reducing end mass
        mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
        
        # Add additional modifications if applicable
        mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
        
        # Add custom reducing end mass for type 7
        if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
            mass += self.custom_reducing_end_mass
        
        return mass
    
    def calculate_composition_mass(self, composition: Dict[str, int]) -> float:
        """
        Calculate the monoisotopic mass of a glycan from its composition dictionary.
        
        Args:
            composition: Dictionary with monosaccharide single-letter codes and counts
                        (e.g., {'H': 5, 'N': 4, 'S': 2} for H5N4S2)
        
        Returns:
            Monoisotopic mass in Daltons
        
        Example:
            >>> calc = GlycanMassCalculator()
            >>> mass = calc.calculate_composition_mass({'H': 5, 'N': 4, 'S': 2})
            >>> print(f"{mass:.4f}")
        """
        # Mapping from single-letter codes to mass constants
        code_to_mass = {
            'H': self.MONO_MASSES['Hex'],
            'N': self.MONO_MASSES['HexNAc'],
            'F': self.MONO_MASSES['Fuc'],
            'S': self.MONO_MASSES['NeuAc'],  # S for Sialic acid
            'G': self.MONO_MASSES.get('NeuGc', 0.0),  # G for NeuGc
            'P': self.MONO_MASSES.get('Phos', 0.0),  # P for Phosphate
        }
        
        mass = 0.0
        for code, count in composition.items():
            if code in code_to_mass:
                mass += count * code_to_mass[code]
        
        # Add reducing end mass (water for Y-type fragments)
        mass += WATER_MASS
        
        return mass
    
    def calculate_mz(self, mass: float, charge: int) -> float:
        """
        Calculate m/z for a given mass and charge state.
        
        Args:
            mass: Neutral mass in Daltons
            charge: Charge state (positive integer)
        
        Returns:
            m/z value
        
        Example:
            >>> calc = GlycanMassCalculator()
            >>> mz = calc.calculate_mz(1000.0, 2)
            >>> print(f"{mz:.4f}")  # ~501.01 m/z
        """
        if charge <= 0:
            raise ValueError("Charge must be a positive integer")
        return (mass + (charge * PROTON_MASS)) / charge
    
    def calculate_peptide_mass(
        self,
        peptide: str,
        use_cam: Optional[bool] = None,
        fixed_mods: Optional[List[str]] = None,
        variable_mods: Optional[List[str]] = None,
        mod_string: Optional[str] = None,
    ) -> float:
        """
        Calculate the monoisotopic mass of a peptide with modifications.
        
        Args:
            peptide: Amino acid sequence (single letter codes)
            use_cam: Apply carbamidomethylation to cysteines (uses instance default if None)
            fixed_mods: List of fixed modifications (uses instance default if None)
            variable_mods: List of variable modifications (uses instance default if None)
            mod_string: Modification string (uses instance default if None)
        
        Returns:
            Monoisotopic mass in Daltons
        
        Example:
            >>> calc = GlycanMassCalculator()
            >>> mass = calc.calculate_peptide_mass("LCPDCPLLAPLNDSR", use_cam=True)
            >>> print(f"{mass:.4f}")
        """
        # Ensure peptide is a string
        peptide = str(peptide or "")
        if not peptide:
            return 0.0
        
        # Use instance defaults if not provided
        use_cam = use_cam if use_cam is not None else self.use_cam
        fixed_mods = fixed_mods if fixed_mods is not None else self.fixed_mods
        variable_mods = variable_mods if variable_mods is not None else self.variable_mods
        mod_string = mod_string if mod_string is not None else self.mod_string
        
        # Create cache key
        cache_key = self._create_cache_key(peptide, use_cam, fixed_mods, variable_mods, mod_string)
        
        # Check cache
        if cache_key in self._peptide_mass_cache:
            return self._peptide_mass_cache[cache_key]
        
        # Calculate mass
        mass = self._calculate_peptide_mass_internal(
            peptide, use_cam, fixed_mods, variable_mods, mod_string
        )
        
        # Cache result
        self._peptide_mass_cache[cache_key] = mass
        
        return mass
    
    def _calculate_peptide_mass_internal(
        self,
        peptide: str,
        use_cam: bool,
        fixed_mods: List[str],
        variable_mods: List[str],
        mod_string: Optional[str],
    ) -> float:
        """Internal peptide mass calculation."""
        # Start with water mass (for complete peptide)
        mass = WATER_MASS
        
        # Determine which positions have modifications
        modifications = apply_modifications_to_positions(
            peptide,
            mod_string=mod_string,
            fixed_mods=fixed_mods,
            variable_mods=variable_mods,
            use_cam=use_cam,
        )
        
        # Store applied modifications for later reference
        self._applied_mods = self._group_modifications(modifications, peptide)
        
        # Calculate mass from amino acids
        for i, aa in enumerate(peptide.upper()):
            # Add amino acid mass
            aa_mass = AMINO_ACID_MASSES.get(aa, 0.0)
            mass += aa_mass
            
            # Add modification mass if present
            if i in modifications:
                mod = modifications[i]
                if isinstance(mod, tuple) and mod[0] == 'CUSTOM':
                    mass += mod[1]
                elif mod in MODIFICATION_MASSES:
                    mass += MODIFICATION_MASSES[mod]
        
        # Handle terminal modifications
        for term in ['N-term', 'C-term']:
            if term in modifications:
                mod = modifications[term]
                if isinstance(mod, tuple) and mod[0] == 'CUSTOM':
                    mass += mod[1]
                elif mod in MODIFICATION_MASSES:
                    mass += MODIFICATION_MASSES[mod]
        
        return mass
    
    def _group_modifications(self, modifications: dict, peptide: str) -> Dict[str, List[str]]:
        """Group modifications by type for reporting."""
        grouped = {}
        for pos, mod in modifications.items():
            if isinstance(mod, tuple):
                mod_name = f"Custom({mod[1]:+.2f})"
            else:
                mod_name = mod
            
            if mod_name not in grouped:
                grouped[mod_name] = []
            
            if isinstance(pos, int):
                aa = peptide[pos] if pos < len(peptide) else '?'
                grouped[mod_name].append(f"{aa}{pos+1}")
            else:
                grouped[mod_name].append(pos)
        
        return grouped
    
    def _create_cache_key(
        self,
        peptide: str,
        use_cam: bool,
        fixed_mods: List[str],
        variable_mods: List[str],
        mod_string: Optional[str],
    ) -> str:
        """Create a cache key for peptide mass calculations."""
        components = [peptide, f"cam:{use_cam}"]
        
        if fixed_mods:
            components.append("fixed:" + ",".join(sorted(fixed_mods)))
        if variable_mods:
            components.append("var:" + ",".join(sorted(variable_mods)))
        if mod_string:
            components.append("mod:" + str(mod_string))
        
        return "||".join(components)
    
    def get_amino_acid_mass(
        self,
        aa: str,
        position: int,
        peptide_sequence: Optional[str] = None,
    ) -> float:
        """
        Get the mass of an amino acid at a specific position, including modifications.
        
        Args:
            aa: Amino acid (single letter code)
            position: Position in the peptide (0-indexed)
            peptide_sequence: Full peptide sequence for context
        
        Returns:
            Mass of the amino acid with any applicable modifications
        """
        base_mass = AMINO_ACID_MASSES.get(aa.upper(), 0.0)
        
        # Check if this position has a modification
        if hasattr(self, '_applied_mods'):
            aa_pos_key = f"{aa.upper()}{position+1}"
            for mod_name, positions in self._applied_mods.items():
                if aa_pos_key in positions:
                    # Extract modification mass
                    if mod_name.startswith("Custom("):
                        # Parse custom mass from string like "Custom(+15.99)"
                        try:
                            custom_val = float(mod_name[7:-1])
                            base_mass += custom_val
                        except ValueError:
                            pass
                    elif mod_name in MODIFICATION_MASSES:
                        base_mass += MODIFICATION_MASSES[mod_name]
                    break
        
        return base_mass
    
    @lru_cache(maxsize=1024)
    def _calculate_base_fragment_mass(self, composition_tuple: tuple) -> float:
        """Calculate base fragment mass with caching."""
        composition = dict(composition_tuple)
        base_mass = 0.0
        
        for mono in ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']:
            if mono in composition:
                base_mass += self.MONO_MASSES[mono] * composition[mono]
        
        return base_mass
    
    def calculate_fragment_mass(
        self,
        composition: Dict[str, int],
        fragment_type: str,
        peptide: Optional[str] = None,
    ) -> float:
        """
        Calculate the mass of a glycan/glycopeptide fragment.
        
        Args:
            composition: Dictionary with monosaccharide counts
                e.g., {'HexNAc': 2, 'Hex': 3}
            fragment_type: Type of fragment:
                'b_ions' - B-type ions (non-reducing end)
                'y_ions' - Y-type ions (reducing end, includes peptide)
                'by_ions' - Internal BY fragments
                'yy_ions' - Double Y cleavage fragments
            peptide: Peptide sequence (for glycopeptide fragments)
        
        Returns:
            Fragment mass in Daltons
        """
        # Use stored peptide if not provided
        peptide = peptide or self.peptide
        
        # Convert composition to hashable tuple for caching
        comp_tuple = tuple(sorted(
            (k, v) for k, v in composition.items() 
            if not k.startswith('_')
        ))
        
        # Get base glycan mass
        base_mass = self._calculate_base_fragment_mass(comp_tuple)
        
        # Check for special flags
        is_glycan_only = composition.get('_glycan_only', False)
        
        # Handle custom mass adjustment
        if '_mass_adjustment' in composition:
            base_mass += composition['_mass_adjustment']
        
        # Adjust based on fragment type
        if is_glycan_only:
            # Glycan-only fragments get reducing end mass
            base_mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
            base_mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
            if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
                base_mass += self.custom_reducing_end_mass
        elif fragment_type in ['b_ions']:
            # B-ions (single cleavage at non-reducing end) - add terminal methylation for permethylated
            if self.modification_type in [2, 3, 5]:
                base_mass += PERMETHYLATION_MASS
        elif fragment_type in ['by_ions']:
            # BY-ions (double cleavage, internal fragment) - no reducing end
            # No additional mass adjustment here; base_mass already accounts for
            # monosaccharide permethylation state for BY fragments.
            pass
        elif fragment_type == 'y_ions':
            # Y-ions: single Y cleavage
            # For permethylated: subtract one methylation
            if self.modification_type in [2, 3, 5]:
                base_mass -= PERMETHYLATION_MASS
            if self.modification_type == 6 and peptide:
                # Glycopeptide mode - peptide replaces reducing end; remove water at linkage
                peptide_mass = self.calculate_peptide_mass(peptide)
                base_mass += peptide_mass
                base_mass -= WATER_MASS
            else:
                # Add reducing end mass and additional modifications
                base_mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
                base_mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
                if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
                    base_mass += self.custom_reducing_end_mass
        elif fragment_type == 'yy_ions':
            # YY-ions: double Y cleavage
            # For permethylated: subtract two methylations
            if self.modification_type in [2, 3, 5]:
                base_mass -= 2 * PERMETHYLATION_MASS
            if self.modification_type == 6 and peptide:
                # Glycopeptide mode - peptide replaces reducing end; remove water at linkage
                peptide_mass = self.calculate_peptide_mass(peptide)
                base_mass += peptide_mass
                base_mass -= WATER_MASS
            else:
                # Add reducing end mass and additional modifications
                base_mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
                base_mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
                if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
                    base_mass += self.custom_reducing_end_mass
        elif fragment_type == 'yyy_ions':
            # YYY-ions: triple Y cleavage
            # For permethylated: subtract three methylations
            if self.modification_type in [2, 3, 5]:
                base_mass -= 3 * PERMETHYLATION_MASS
            if self.modification_type == 6 and peptide:
                # Glycopeptide mode - peptide replaces reducing end; remove water at linkage
                peptide_mass = self.calculate_peptide_mass(peptide)
                base_mass += peptide_mass
                base_mass -= WATER_MASS
            else:
                # Add reducing end mass and additional modifications
                base_mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
                base_mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
                if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
                    base_mass += self.custom_reducing_end_mass
        elif fragment_type == 'byy_ions':
            # BYY-ions: Internal fragments (both ends cleaved, no reducing end)
            # For permethylated glycans: subtract ONE methylation at the cleavage site
            if self.modification_type in [2, 3, 5]:
                base_mass -= PERMETHYLATION_MASS
        elif fragment_type == 'byyy_ions':
            # BYYY-ions: Internal fragments (no reducing end)
            # For permethylated glycans: subtract ONE methylation at the cleavage site
            if self.modification_type in [2, 3, 5]:
                base_mass -= PERMETHYLATION_MASS
        elif fragment_type in ['cz_ions', 'czz_ions', 'czzz_ions', 'bzzz_ions']:
            # CZ-series internal fragments (CZ, CZZ, CZZZ, BZZZ)
            # Same as BY-series: subtract ONE methylation at the cleavage site
            if self.modification_type in [2, 3, 5]:
                base_mass -= PERMETHYLATION_MASS
        elif fragment_type in ['c_ions', 'z_ions', 'zz_ions', 'zzz_ions']:
            # C and Z-series fragments parallel B and Y-series
            if fragment_type == 'c_ions':
                # C-ions (single cleavage at non-reducing end, like B-ions)
                if self.modification_type in [2, 3, 5]:
                    base_mass += PERMETHYLATION_MASS
            elif fragment_type in ['z_ions', 'zz_ions', 'zzz_ions']:
                # Z-ions (reducing end like Y-ions)
                # Subtract methylations for multiple cleavages
                if fragment_type == 'z_ions' and self.modification_type in [2, 3, 5]:
                    base_mass -= PERMETHYLATION_MASS
                elif fragment_type == 'zz_ions' and self.modification_type in [2, 3, 5]:
                    base_mass -= 2 * PERMETHYLATION_MASS
                elif fragment_type == 'zzz_ions' and self.modification_type in [2, 3, 5]:
                    base_mass -= 3 * PERMETHYLATION_MASS
                
                # Add reducing end mass for Z-ions (like Y-ions)
                if self.modification_type == 6 and peptide:
                    # Glycopeptide mode
                    peptide_mass = self.calculate_peptide_mass(peptide)
                    base_mass += peptide_mass
                    base_mass -= WATER_MASS
                else:
                    base_mass += REDUCING_END_MASSES.get(self.modification_type, 0.0)
                    base_mass += ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
                    if self.modification_type == 7 and self.custom_reducing_end_mass is not None:
                        base_mass += self.custom_reducing_end_mass
        
        return base_mass
    
    def calculate_glycopeptide_mass(
        self,
        peptide: str,
        glycan_code: str,
    ) -> float:
        """
        Calculate the mass of a glycopeptide (peptide + glycan).
        
        Args:
            peptide: Peptide sequence
            glycan_code: Glycan composition code
        
        Returns:
            Glycopeptide mass in Daltons
        
        Note:
            The glycopeptide mass is: peptide_mass + glycan_mass - water
            (water is lost during glycosidic bond formation)
        """
        peptide_mass = self.calculate_peptide_mass(peptide)
        
        # Calculate glycan mass without reducing end (since peptide provides it)
        hexnac, hex_count, fuc, neuac, neugc = self.parse_glycan_code(glycan_code)
        
        glycan_mass = 0.0
        glycan_mass += hexnac * self.MONO_MASSES['HexNAc']
        glycan_mass += hex_count * self.MONO_MASSES['Hex']
        glycan_mass += fuc * self.MONO_MASSES['Fuc']
        glycan_mass += neuac * self.MONO_MASSES['NeuAc']
        glycan_mass += neugc * self.MONO_MASSES['NeuGc']
        
        # Glycopeptide: peptide + glycan (no additional water since glycan
        # masses are already residue masses)
        return peptide_mass + glycan_mass
    
    def reset_state(self) -> None:
        """Reset calculator state (clears caches)."""
        self._peptide_mass_cache.clear()
        self._applied_mods.clear()
    
    def get_applied_modifications(self) -> Dict[str, List[str]]:
        """
        Get the modifications that were applied in the last calculation.
        
        Returns:
            Dictionary mapping modification names to list of positions
        """
        return self._applied_mods.copy()


# =============================================================================
# Convenience Functions
# =============================================================================

def calculate_glycan_mass(glycan_code: str, modification_type: int = 0) -> float:
    """
    Convenience function to calculate glycan mass.
    
    Args:
        glycan_code: Glycan composition code
        modification_type: Reducing end modification type (0-6)
    
    Returns:
        Glycan mass in Daltons
    """
    calc = GlycanMassCalculator(modification_type=modification_type)
    return calc.calculate_glycan_mass(glycan_code)


def calculate_peptide_mass(
    peptide: str,
    use_cam: bool = True,
    modifications: Optional[List[str]] = None,
) -> float:
    """
    Convenience function to calculate peptide mass.
    
    Args:
        peptide: Peptide sequence
        use_cam: Apply carbamidomethylation to cysteines
        modifications: List of modifications (e.g., ["Ox:M"])
    
    Returns:
        Peptide mass in Daltons
    """
    calc = GlycanMassCalculator(
        use_cam=use_cam,
        variable_mods=modifications,
    )
    return calc.calculate_peptide_mass(peptide)


def calculate_mz(mass: float, charge: int) -> float:
    """
    Convenience function to calculate m/z.
    
    Args:
        mass: Neutral mass in Daltons
        charge: Charge state
    
    Returns:
        m/z value
    """
    return (mass + (charge * PROTON_MASS)) / charge
