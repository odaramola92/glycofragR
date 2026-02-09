"""
Peptide class for peptide fragmentation analysis.

This module provides the main Peptide class for representing peptide sequences
and generating fragment ions (b, y, c, z series).
"""

from typing import Dict, List, Optional, Tuple, Any

from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.constants import PROTON_MASS, WATER_MASS


class Peptide:
    """
    Represents a peptide sequence with methods for fragmentation.
    
    This class handles:
    - Parsing peptide sequences
    - Applying modifications (CAM, variable mods, etc.)
    - Generating b-ions and y-ions
    - Generating c-ions and z-ions
    - Calculating fragment masses

    Public API:
        - generate_b_ions()
        - generate_y_ions(...)
        - generate_c_ions()
        - generate_z_ions()
        - generate_all_fragments(...)
        - generate_all_fragments_with_losses()
    
    Attributes:
        sequence (str): The peptide amino acid sequence
        modifications (dict): Applied modifications by position
        use_cam (bool): Whether CAM is applied to cysteines
        calculator (GlycanMassCalculator): Mass calculator instance
    
    Example:
        >>> peptide = Peptide('PEPTIDE', use_cam=True)
        >>> b_ions = peptide.generate_b_ions()
        >>> y_ions = peptide.generate_y_ions()
        >>> print(f"Generated {len(b_ions)} b-ions and {len(y_ions)} y-ions")
    """
    
    def __init__(
        self,
        sequence: str,
        use_cam: bool = True,
        fixed_mods: Optional[List[str]] = None,
        mod_string: Optional[str] = None,
        modification_type: int = 6
    ):
        """
        Initialize a Peptide object.
        
        Args:
            sequence: Amino acid sequence (single-letter codes)
            use_cam: Whether to apply carbamidomethylation to cysteines
            fixed_mods: List of fixed modifications (for internal use)
            mod_string: Targeted modification string. Supports 29 modifications.
                        Use list_supported_modifications() to see all available modifications.
                        
                        Examples:
                            - "M:Ox" (oxidize all M)
                            - "M4,24:Ox" (oxidize M at positions 4 and 24)
                            - "N-term:Ac" (N-terminal acetylation)
                            - "M:Ox; K:Ac" (multiple modifications)
                            - "Custom:A, 45.8978" (add custom mass to all A)
                            
            modification_type: Reducing end modification type (default 6 for glycopeptide)
        
        Note:
            To discover all supported modifications:
            >>> from glycofrag import list_supported_modifications
            >>> list_supported_modifications(verbose=True)
        """
        self.sequence = sequence.upper()
        self.use_cam = use_cam
        self.fixed_mods = fixed_mods or []
        self.mod_string = mod_string
        self.modification_type = modification_type
        
        # Create calculator instance
        self.calculator = GlycanMassCalculator(
            modification_type=modification_type,
            use_cam=use_cam,
            fixed_mods=fixed_mods,
            mod_string=mod_string,
            peptide=sequence
        )
        
        # Calculate peptide mass to apply modifications
        self.mass = self.calculator.calculate_peptide_mass(
            sequence,
            use_cam=use_cam,
            fixed_mods=fixed_mods,
            mod_string=mod_string
        )
    
    @property
    def residue_mass(self) -> float:
        """
        Get the peptide mass without water (residue mass only).
        
        Useful for glycopeptides where water is released at the glycosidic bond.
        
        Returns:
            Mass of peptide residues without water molecule
        """
        from glycofrag.core.constants import WATER_MASS
        return self.mass - WATER_MASS
        
        # Store applied modifications (if available from calculator)
        self.modifications = {}
    
    def generate_b_ions(self) -> List[Dict[str, any]]:
        """
        Generate b-ion fragments (N-terminal fragments).
        
        b-ions are formed by cleavage at the peptide bond with charge
        retention on the N-terminal fragment.
        
        Returns:
            List of dictionaries containing fragment information:
            - fragment_name: e.g., 'b1', 'b2', etc.
            - fragment_mass: Mass in Da (with proton)
            - fragment_sequence: Amino acid sequence
            - fragment_position: Position in sequence
            - fragment_type: 'b'
        
        Example:
            >>> peptide = Peptide('PEPTIDE')
            >>> b_ions = peptide.generate_b_ions()
            >>> print(b_ions[0])  # First b-ion (b1)
            {'fragment_name': 'b1', 'fragment_mass': 98.06..., ...}
        """
        fragments = []
        cumulative_mass = 0.0
        
        for i, aa in enumerate(self.sequence):
            # Get amino acid mass with modifications
            aa_mass = self.calculator.get_amino_acid_mass(aa, i, self.sequence)
            cumulative_mass += aa_mass
            
            # Create b-ion
            b_ion = {
                'fragment_name': f"b{i+1}",
                'fragment_mass': cumulative_mass + PROTON_MASS,
                'fragment_sequence': self.sequence[:i+1],
                'fragment_position': i+1,
                'fragment_type': 'b',
                'position': i+1,
                'mass': cumulative_mass + PROTON_MASS,
                'sequence': self.sequence[:i+1]
            }
            fragments.append(b_ion)
        
        return fragments
    
    def generate_y_ions(self, include_full_peptide: bool = True) -> List[Dict[str, any]]:
        """
        Generate y-ion fragments (C-terminal fragments).
        
        y-ions are formed by cleavage at the peptide bond with charge
        retention on the C-terminal fragment (includes +H2O).
        
        Args:
            include_full_peptide: Whether to include Y0 (full peptide + H2O)
        
        Returns:
            List of dictionaries containing fragment information:
            - fragment_name: e.g., 'y1', 'y2', 'Y0' (full peptide)
            - fragment_mass: Mass in Da (with proton)
            - fragment_sequence: Amino acid sequence
            - fragment_position: Position from C-terminus
            - fragment_type: 'y'
        
        Example:
            >>> peptide = Peptide('PEPTIDE')
            >>> y_ions = peptide.generate_y_ions()
            >>> print(y_ions[0])  # First y-ion from C-terminus
        """
        fragments = []
        cumulative_mass = 18.010564684  # H2O mass
        
        for i in range(len(self.sequence)-1, -1, -1):
            aa = self.sequence[i]
            aa_mass = self.calculator.get_amino_acid_mass(aa, i, self.sequence)
            cumulative_mass += aa_mass
            
            y_position = len(self.sequence) - i
            
            # Name Y0 for full peptide, otherwise y<n>
            if y_position == len(self.sequence):
                fragment_name = "Y0" if include_full_peptide else f"y{y_position}"
            else:
                fragment_name = f"y{y_position}"
            
            y_ion = {
                'fragment_name': fragment_name,
                'fragment_mass': cumulative_mass + PROTON_MASS,
                'fragment_sequence': self.sequence[i:],
                'fragment_position': y_position,
                'fragment_type': 'y',
                'position': y_position,
                'mass': cumulative_mass + PROTON_MASS,
                'sequence': self.sequence[i:]
            }
            
            # Skip Y0 if not including full peptide
            if y_position == len(self.sequence) and not include_full_peptide:
                continue
            
            fragments.append(y_ion)
        
        return fragments
    
    def generate_c_ions(self) -> List[Dict[str, any]]:
        """
        Generate c-ion fragments (N-terminal fragments with +NH3).
        
        c-ions are similar to b-ions but include an additional NH3 group.
        
        Returns:
            List of dictionaries containing fragment information
        
        Example:
            >>> peptide = Peptide('PEPTIDE')
            >>> c_ions = peptide.generate_c_ions()
        """
        fragments = []
        cumulative_mass = 17.026549  # NH3 mass
        
        for i, aa in enumerate(self.sequence):
            aa_mass = self.calculator.get_amino_acid_mass(aa, i, self.sequence)
            cumulative_mass += aa_mass
            
            c_ion = {
                'fragment_name': f"c{i+1}",
                'fragment_mass': cumulative_mass + PROTON_MASS,
                'fragment_sequence': self.sequence[:i+1],
                'fragment_position': i+1,
                'fragment_type': 'c',
                'position': i+1,
                'mass': cumulative_mass + PROTON_MASS,
                'sequence': self.sequence[:i+1]
            }
            fragments.append(c_ion)
        
        return fragments
    
    def generate_z_ions(self) -> List[Dict[str, any]]:
        """
        Generate z-ion fragments (C-terminal fragments with -NH).
        
        z-ions are similar to y-ions but with loss of NH group.
        
        Returns:
            List of dictionaries containing fragment information
        
        Example:
            >>> peptide = Peptide('PEPTIDE')
            >>> z_ions = peptide.generate_z_ions()
        """
        fragments = []
        cumulative_mass = 18.010564684 - 15.010899  # H2O - NH
        
        for i in range(len(self.sequence)-1, -1, -1):
            aa = self.sequence[i]
            aa_mass = self.calculator.get_amino_acid_mass(aa, i, self.sequence)
            cumulative_mass += aa_mass
            
            z_position = len(self.sequence) - i
            
            z_ion = {
                'fragment_name': f"z{z_position}",
                'fragment_mass': cumulative_mass + PROTON_MASS,
                'fragment_sequence': self.sequence[i:],
                'fragment_position': z_position,
                'fragment_type': 'z',
                'position': z_position,
                'mass': cumulative_mass + PROTON_MASS,
                'sequence': self.sequence[i:]
            }
            fragments.append(z_ion)
        
        return fragments
    
    def generate_all_fragments(self, fragment_types: Optional[List[str]] = None) -> Dict[str, List[Dict[str, Any]]]:
        """
        Generate peptide fragment ions based on specified series.
        
        Args:
            fragment_types: List of fragment series to generate. Options:
                - 'by' - Generate b and y ions (default)
                - 'cz' - Generate c and z ions
                Can specify both: ['by', 'cz'] to get all types
        
        Returns:
            Dictionary of fragment ion lists organized by type:
            - 'b_ions': List of b-ion fragments
            - 'y_ions': List of y-ion fragments
            - 'c_ions': List of c-ion fragments (if 'cz' requested)
            - 'z_ions': List of z-ion fragments (if 'cz' requested)
            
            Each fragment dictionary contains:
            - name (str): Fragment name (e.g., 'b1', 'y3', 'c2', 'z5')
            - type (str): Fragment type ('b', 'y', 'c', or 'z')
            - mass (float): m/z value of fragment
            - charge (int): Charge state (typically 1)
            - sequence (str): Amino acid sequence of fragment
            - position (int): Position in original peptide
            - metadata (dict): Additional information
        
        Example:
            >>> peptide = Peptide('LCPDCPLLAPLNDSR', modifications={'CAM': [1, 5]})
            >>> # Generate only b and y ions
            >>> fragments = peptide.generate_all_fragments(fragment_types=['by'])
            >>> # Generate all types
            >>> all_frags = peptide.generate_all_fragments(fragment_types=['by', 'cz'])
        """
        # Default to 'by' series
        if fragment_types is None:
            fragment_types = ['by']
        
        # Normalize to lowercase
        fragment_types = [ft.lower() for ft in fragment_types]
        
        result = {}
        
        # Generate requested series
        if 'by' in fragment_types:
            result['b_ions'] = self.generate_b_ions()
            result['y_ions'] = self.generate_y_ions()
        
        if 'cz' in fragment_types:
            result['c_ions'] = self.generate_c_ions()
            result['z_ions'] = self.generate_z_ions()
        
        return result
    
    def get_fragment_by_name(self, fragment_name: str) -> Optional[Dict[str, any]]:
        """
        Get a specific fragment by name.
        
        Args:
            fragment_name: Fragment name (e.g., 'b3', 'y2', 'Y0')
        
        Returns:
            Fragment dictionary or None if not found
        
        Example:
            >>> peptide = Peptide('PEPTIDE')
            >>> fragment = peptide.get_fragment_by_name('b3')
            >>> print(fragment['fragment_mass'])
        """
        all_fragments = self.generate_all_fragments(fragment_types=['by', 'cz'])
        
        for fragment_type in ['b_ions', 'y_ions', 'c_ions', 'z_ions']:
            for fragment in all_fragments.get(fragment_type, []):
                if fragment.get('fragment_name') == fragment_name:
                    return fragment
        
        return None
    
    def _generate_neutral_losses(self) -> List[Dict[str, Any]]:
        """
        Generate peptide fragment neutral loss variants.
        
        Generates the following neutral losses:
        - H2O loss (18.0106 Da): from S, T, E, D residues
        - NH3 loss (17.0265 Da): from K, R, N, Q residues
        - CO2 loss (43.9898 Da): from E, D residues
        - CO loss (27.9949 Da): from b-ions only (a-ion series)
        - CH4OS loss (64.0034 Da): from oxidized M residues only
        
        Returns:
            List of neutral loss fragment dictionaries
        """
        neutral_losses = []
        
        # Define neutral loss masses (in Da)
        NEUTRAL_LOSS_MASSES = {
            'H2O': WATER_MASS,  # 18.0105647
            'NH3': 17.0265491,
            'CO2': 43.9898922,
            'CO': 27.9949146,
            'CH4OS': 64.00324,  # For oxidized methionine
        }
        
        # Residues that can lose H2O
        h2o_residues = {'S', 'T', 'E', 'D'}
        
        # Residues that can lose NH3
        nh3_residues = {'K', 'R', 'N', 'Q'}
        
        # Residues that can lose CO2
        co2_residues = {'E', 'D'}
        
        # Residues that can lose CH4OS (methylsulfenic acid from oxidized Met)
        ch4os_residues = {'M'}
        
        all_fragments = self.generate_all_fragments()

        # Track oxidized Met positions from applied modifications (1-indexed)
        oxidized_positions = set()
        try:
            applied_mods = self.calculator.get_applied_modifications()
        except Exception:
            applied_mods = getattr(self.calculator, '_applied_mods', {}) or {}

        for mod_name, positions in (applied_mods or {}).items():
            if mod_name == 'Ox':
                for pos in positions:
                    if isinstance(pos, str):
                        digits = "".join(ch for ch in pos if ch.isdigit())
                        if digits:
                            oxidized_positions.add(int(digits))

        def _fragment_has_oxidized_m(frag_type: str, frag_pos: Optional[int]) -> bool:
            if not oxidized_positions or not frag_pos:
                return False
            seq_len = len(self.sequence)
            if frag_type == 'b_ions':
                # b-ions cover positions 1..frag_pos
                return any(pos <= frag_pos for pos in oxidized_positions)
            if frag_type == 'y_ions':
                # y-ions cover positions (seq_len - frag_pos + 1)..seq_len
                start = seq_len - frag_pos + 1
                return any(pos >= start for pos in oxidized_positions)
            return False
        
        def _make_loss(fragment, frag_type, loss_name, loss_mass):
            """Create a neutral loss variant of a fragment."""
            ion_prefix = frag_type.replace('_ions', '')
            frag_seq = fragment.get('fragment_sequence', fragment.get('sequence', ''))
            loss_frag = fragment.copy()
            loss_frag['fragment_name'] = f"{fragment.get('fragment_name')}(-{loss_name})"
            loss_frag['fragment_mass'] = fragment.get('fragment_mass', 0) - loss_mass
            loss_frag['mass'] = loss_frag['fragment_mass']
            loss_frag['fragment_type'] = f"{ion_prefix}-{loss_name}"
            loss_frag['type'] = loss_frag['fragment_type']
            # Store loss in fragment_sequence for Composition column (matches GlypPRM)
            loss_frag['fragment_sequence'] = f"{frag_seq}(-{loss_name})"
            loss_frag['sequence'] = loss_frag['fragment_sequence']
            loss_frag['_is_custom'] = True
            loss_frag['_custom_source'] = f"{loss_name.lower()}_loss"
            loss_frag['metadata'] = {'neutral_loss': loss_name, 'source': frag_type}
            return loss_frag
        
        # Process all b and y ions
        for frag_type in ['b_ions', 'y_ions']:
            for fragment in all_fragments[frag_type]:
                seq = fragment.get('fragment_sequence', fragment.get('sequence', ''))
                if not seq:
                    continue
                
                # 1. H2O loss — from fragments containing S, T, E, or D (b and y)
                if any(aa in seq for aa in h2o_residues):
                    neutral_losses.append(_make_loss(fragment, frag_type, 'H2O', NEUTRAL_LOSS_MASSES['H2O']))
                
                # 2. NH3 loss — from b-ions containing K, R, N, or Q (b-only)
                if frag_type == 'b_ions' and any(aa in seq for aa in nh3_residues):
                    neutral_losses.append(_make_loss(fragment, frag_type, 'NH3', NEUTRAL_LOSS_MASSES['NH3']))
                
                # 3. CO2 loss — from fragments containing E or D (b and y)
                if any(aa in seq for aa in co2_residues):
                    neutral_losses.append(_make_loss(fragment, frag_type, 'CO2', NEUTRAL_LOSS_MASSES['CO2']))
                
                # 4. CH4OS loss — only when oxidized Met is present in this fragment
                frag_pos = fragment.get('position', fragment.get('fragment_position'))
                if _fragment_has_oxidized_m(frag_type, frag_pos):
                    neutral_losses.append(_make_loss(fragment, frag_type, 'CH4OS', NEUTRAL_LOSS_MASSES['CH4OS']))

                # 5. CO loss — b-ions only (a-ion series)
                if frag_type == 'b_ions':
                    neutral_losses.append(_make_loss(fragment, frag_type, 'CO', NEUTRAL_LOSS_MASSES['CO']))
        
        return neutral_losses
    
    def generate_all_fragments_with_losses(self) -> List[Dict[str, Any]]:
        """
        Generate all fragment types including neutral loss variants.
        
        This method generates all b, y, c, z ions plus neutral loss variants.
        
        Returns:
            List of all fragment dictionaries including neutral loss variants
        
        Example:
            >>> peptide = Peptide('LCPDCPLLAPLNDSR')
            >>> frags = peptide.generate_all_fragments_with_losses()
            >>> h2o_losses = [f for f in frags if 'H2O' in f['fragment_name']]
            >>> print(f"Total fragments: {len(frags)}")
            >>> print(f"H2O loss variants: {len(h2o_losses)}")
        """
        all_frags = self.generate_all_fragments()
        
        # Flatten the dictionary into a list
        flat_frags = []
        for frag_list in all_frags.values():
            flat_frags.extend(frag_list)
        
        # Add neutral loss variants
        neutral_loss_frags = self._generate_neutral_losses()
        flat_frags.extend(neutral_loss_frags)
        
        return flat_frags
    
    def __repr__(self) -> str:
        """String representation of the peptide."""
        return f"Peptide('{self.sequence}', mass={self.mass:.4f}, use_cam={self.use_cam})"
    
    def __len__(self) -> int:
        """Return the length of the peptide sequence."""
        return len(self.sequence)
