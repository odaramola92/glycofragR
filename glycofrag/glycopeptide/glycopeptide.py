"""
Glycopeptide class for combined peptide and glycan fragmentation.

This module provides comprehensive glycopeptide analysis including:
- Peptide backbone fragments (b, y, c, z ions)
- Glycan fragments with peptide (B, Y, BY, YY, C, Z ions with -PEP suffix)
- Glycan structure visualization (predicted 3D structures)
- Intact glycopeptide ions

Public API:
    - Glycopeptide(...)
    - generate_fragments()
    - visualize_structure()

Typical usage:
    >>> gp = Glycopeptide("LCPDCPLLAPLNDSR", "4501", 12, "N")
    >>> fragments = gp.generate_fragments()
    >>> gp.visualize_structure(structure_number=1)  # Visualize glycan structure
"""

from typing import Dict, List, Tuple, Optional, Any, Union
import networkx as nx

from glycofrag import Peptide, Glycan
from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.constants import WATER_MASS
from glycofrag.core.modifications import get_modification_type
from glycofrag.io.visualizer import GlycanVisualizer


class Glycopeptide:
    """
    Represents a glycopeptide with combined peptide and glycan fragmentation.
    
    This class integrates peptide and glycan analysis to generate:
    - Peptide fragments: b, y, c, z ions (backbone only)
    - Glycan fragments: B, Y, BY, YY, C, Z ions (all with -PEP suffix indicating peptide attachment)
    - Intact: Complete glycopeptide
    
    Design Standard:
        The reducing end modification is ALWAYS 'glycopeptide' mode, meaning the reducing end
        HexNAc is modified by attachment to the peptide. This is the standard for glycopeptides.
        This is NOT configurable - future versions may support other glycan modifications (v2+).
    
    Public API:
        - generate_fragments(...)     -> dict of fragment lists
        - visualize_structure(...)    -> visualize predicted glycan structure
        - __repr__()                  -> readable summary of the object

    To get fragments, call `generate_fragments()` on a `Glycopeptide` instance.
    
    Attributes:
        peptide_sequence (str): Amino acid sequence
        glycan_code (str): Glycan composition code
        glycosylation_site (int): Position of glycan attachment (1-indexed)
        peptide (Peptide): Peptide object
        glycan (Glycan): Glycan object
        mass_calculator (GlycanMassCalculator): Mass calculator instance
    
    Example:
        >>> gp = Glycopeptide(
        ...     peptide_sequence='EEQYNSTYR',
        ...     glycan_code='4501',
        ...     glycosylation_site=5,  # N in NxS/T motif
        ...     glycan_type='N'
        ... )
        >>> fragments = gp.generate_fragments()
        >>> print(f"Peptide fragments: {len(fragments['peptide_b_ions']) + len(fragments['peptide_y_ions'])}")
        >>> print(f"Glycan fragments: {len(fragments['b_ions']) + len(fragments['y_ions']) + len(fragments['by_ions'])}")
    """
    
    def __init__(
        self,
        peptide_sequence: str,
        glycan_code: str,
        glycosylation_site: int,
        glycan_type: str,
        use_cam: bool = True,
        max_structures: int = 100,
        mod_string: Optional[str] = None,
        isomer_sensitive: bool = False,
        custom_mono_masses: Optional[Dict[str, float]] = None,
        preferred_core: Optional[int] = None
    ):
        """
        Initialize a Glycopeptide.
        
        Args:
            peptide_sequence: Amino acid sequence (e.g., 'EEQYNSTYR')
            glycan_code: Glycan composition code (e.g., '4501')
            glycosylation_site: 1-indexed position of glycan attachment
            glycan_type: 'N' for N-glycan or 'O' for O-glycan
            use_cam: Whether to apply carbamidomethylation to cysteines
            max_structures: Maximum glycan structures to predict
            mod_string: Targeted peptide modification string. Supports 29 modifications.
                        Use list_supported_modifications() to see all available modifications.
                        
                        Examples:
                            - "M:Ox" (oxidize all M)
                            - "M4,24:Ox" (oxidize M at positions 4 and 24)
                            - "N-term:Ac" (N-terminal acetylation)
                            - "M:Ox; K:Ac; N-term:Ac" (multiple modifications)
                            - "S3:Phos; T8:Phos" (phosphorylation at specific sites)
            isomer_sensitive: If True, treats mirror images as distinct structures.
                            If False (default), deduplicates topologically identical
                            mirror images for more concise results.
            custom_mono_masses: Dict of mass (Da) to ADD to specific monosaccharides.
                              Keys: 'HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc'
                              Values: Mass to add (can be positive or negative)
                              Examples:
                                - {'NeuAc': 42.0106} for O-acetylated NeuAc
                                - {'NeuAc': 50.0657} for custom modification
                                - {'Hex': -18.0106} for mass loss (rare)
            preferred_core: For O-glycans with structurally ambiguous cores.
                           When fragments cannot distinguish between cores (e.g., Core 1 vs 8),
                           specify which to use: 1, 3, 5, 6, 7, or 8.
                           Defaults to lowest core number if None.
        
        Note:
            Glycopeptide fragments always use the "glycopeptide" reducing end mode,
            which means the reducing end HexNAc is modified by attachment to the peptide.
            This is NOT configurable - it is the standard for glycopeptide analysis.
            
            Future versions may support other glycan modifications (v2+).
            
            To discover all supported modifications and their targets:
            >>> from glycofrag import list_supported_modifications
            >>> list_supported_modifications(verbose=True)
            
        Example:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N')
            >>> # With oxidation:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N', mod_string='M:Ox')
            >>> # With N-terminal acetylation:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N', mod_string='N-term:Ac')
            >>> # Generate BY series (b/y peptide ions + B/Y glycan ions)
            >>> fragments = gp.generate_fragments(structure_index=1, peptide_fragment_types=['by'], glycan_fragment_types=['BY'])
            >>> # Generate all series (BY+CZ glycan + by+cz peptide)
            >>> fragments = gp.generate_fragments(structure_index=1, peptide_fragment_types=['by', 'cz'], glycan_fragment_types=['BY', 'CZ'])
            >>> # Visualize the glycan structure
            >>> gp.visualize_structure(structure_index=1)
        """
        self.peptide_sequence = peptide_sequence.upper()
        self.glycan_code = glycan_code
        self.glycosylation_site = glycosylation_site
        self._glycosylation_index0 = glycosylation_site - 1
        self.glycan_type = glycan_type
        self.use_cam = use_cam
        
        # Glycopeptides always use glycopeptide reducing end mode
        self.modification_type = get_modification_type('glycopeptide')
        
        # Validate glycosylation site (1-indexed)
        if glycosylation_site < 1 or glycosylation_site > len(peptide_sequence):
            raise ValueError(f"Glycosylation site {glycosylation_site} out of range for peptide length {len(peptide_sequence)}")
        
        # Initialize peptide and glycan
        self.peptide = Peptide(
            peptide_sequence,
            use_cam=use_cam,
            mod_string=mod_string
        )
        self.glycan = Glycan(
            glycan_code,
            glycan_type=glycan_type,
            max_structures=max_structures,
            modification_type=self.modification_type,
            isomer_sensitive=isomer_sensitive,
            custom_mono_masses=custom_mono_masses,
            _allow_glycopeptide=True,
            preferred_core=preferred_core
        )
        self.mass_calculator = GlycanMassCalculator(modification_type=self.modification_type)
        
        # Predict glycan structures
        self.glycan_structures = self.glycan.predict_structures()
    
    @property
    def mass(self) -> float:
        """
        Calculate the monoisotopic mass of the glycopeptide.
        
        When peptide and glycan combine at the glycosidic bond, a water molecule is released.
        Therefore: Glycopeptide Mass = Peptide Mass + Glycan Mass - Water
        
        Returns:
            Neutral monoisotopic mass in Daltons
        """
        from glycofrag.core.constants import WATER_MASS
        peptide_mass = self.peptide.mass
        glycan_mass = self.mass_calculator.calculate_glycan_mass(self.glycan_code)
        return peptide_mass + glycan_mass - WATER_MASS
    
    def generate_fragments(
        self,
        structure_index: Optional[int] = None,
        fragment_types: Optional[List[str]] = None,
        peptide_fragment_types: Optional[List[str]] = None,
        glycan_fragment_types: Optional[List[str]] = None,
        charges: Optional[List[int]] = None
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Generate all glycopeptide fragment ions from one or all predicted structures.
        
        Fragment types:
        - Peptide fragments: b, y, c, z ions (backbone fragments)
          Generated from series: 'by' (b+y), 'cz' (c+z)
        - Glycan fragments: B, Y, BY, YY, C, Z, CZ, ZZ ions
          Generated from series: 'BY', 'CZ'
        - Y-ions only have -PEP suffix (indicating peptide at reducing end)
        
        Args:
            structure_index: Which predicted glycan structure to use, 1-indexed (default None = all structures).
                          If None, generates fragments from ALL structures and tracks sources.
                          If int, generates fragments from that specific structure only.
            fragment_types: List of fragment types to generate (legacy parameter, deprecated)
                          Note: Y0, Y1, Y-glycan are false glycopeptide constructs and are no longer generated.
                          They are implicitly represented by peptide y-ions and glycan Y-ions.
                          Currently only 'intact' and 'peptide-b' are supported.
            peptide_fragment_types: List of peptide fragment series: ['by'] or ['cz'] or ['by', 'cz']
                          - 'by' generates: b_ions, peptide_b_ions, y_ions, peptide_y_ions
                          - 'cz' generates: c_ions, peptide_c_ions, z_ions, peptide_z_ions
            glycan_fragment_types: List of glycan fragment series: ['BY'] or ['CZ'] or ['BY', 'CZ']
                          - 'BY' generates: b_ions, y_ions, by_ions, yy_ions
                          - 'CZ' generates: c_ions, z_ions, cz_ions, zz_ions
                          A-type oxonium and custom diagnostic ions are always generated automatically.
            charges: Charge states for m/z calculation (default [1, 2, 3])
                    Note: This parameter is currently unused in generate_fragments()
                    but is stored for use in downstream table generation.
        
        Returns:
            Dict with fragment ion lists organized by type. Each fragment includes 'structure' field
            listing which structures produced it (e.g., "1", "1,2,3", "2").
        
        Example:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N')
            >>> # Generate fragments from ALL structures (default)
            >>> frags = gp.generate_fragments()
            >>> # Generate fragments from only structure 1
            >>> frags = gp.generate_fragments(structure_index=1)
            >>> # Generate BY + CZ glycan series with by + cz peptide series from all structures
            >>> frags = gp.generate_fragments(peptide_fragment_types=['by', 'cz'], glycan_fragment_types=['BY', 'CZ'])
        """
        
        if fragment_types is None:
            fragment_types = ['peptide-b', 'peptide-y', 'Y1']
        
        if peptide_fragment_types is None:
            peptide_fragment_types = ['by']  # Default to by series
        
        if glycan_fragment_types is None:
            glycan_fragment_types = ['BY']  # Default to BY series (A-type/custom ions are automatic)
        
        fragment_types_upper = [str(ft).upper() for ft in fragment_types]

        # Determine which structures to process
        if structure_index is None:
            # Use all structures
            structures_to_process = list(range(1, len(self.glycan_structures) + 1))
        else:
            # Use single structure (validate it)
            if structure_index < 1:
                raise ValueError(f"Structure index must be >= 1 (got {structure_index})")
            if structure_index > len(self.glycan_structures):
                raise ValueError(f"Structure index {structure_index} out of range (have {len(self.glycan_structures)} structures, valid range: 1-{len(self.glycan_structures)})")
            structures_to_process = [structure_index]
        
        # Initialize fragment container
        all_fragments = {
            'b_ions': [],       # Glycan B-ions with peptide (-PEP)
            'y_ions': [],       # Glycan Y-ions with peptide (-PEP)
            'by_ions': [],      # Glycan BY-ions with peptide (-PEP)
            'yy_ions': [],      # Glycan YY-ions with peptide (-PEP)
            'byy_ions': [],     # Glycan BYY-ions with peptide (-PEP)
            'yyy_ions': [],     # Glycan YYY-ions with peptide (-PEP)
            'byyy_ions': [],    # Glycan BYYY-ions with peptide (-PEP)
            'c_ions': [],       # Glycan C-ions with peptide (-PEP)
            'z_ions': [],       # Glycan Z-ions with peptide (-PEP)
            'cz_ions': [],      # Glycan CZ-ions with peptide (-PEP)
            'zz_ions': [],      # Glycan ZZ-ions with peptide (-PEP)
            'czz_ions': [],     # Glycan CZZ-ions with peptide (-PEP)
            'zzz_ions': [],     # Glycan ZZZ-ions with peptide (-PEP)
            'czzz_ions': [],    # Glycan CZZZ-ions with peptide (-PEP)
            'bzzz_ions': [],    # Glycan BZZZ-ions with peptide (-PEP)
            'y1_ions': [],      # Y1: peptide + HexNAc1 (and Fuc1 if core fucose)
            'peptide_b_ions': [],# Peptide b-ions (backbone only, no glycan possible)
            'peptide_y_ions': [],# Peptide y-ions (backbone only, no glycan possible)
            'peptide_c_ions': [],# Peptide c-ions (backbone only, no glycan possible)
            'peptide_z_ions': [],# Peptide z-ions (backbone only, no glycan possible)
            'peptide_neutral_losses': [],  # Peptide neutral loss variants (H2O, NH3, CO2, CO)
            'a_ions': [],       # A-type oxonium diagnostic ions
            'custom_ions': [],  # Custom diagnostic ions (neutral losses from mono fragments)
            'intact': []        # Complete glycopeptide
        }
        
        # Generate fragments from each structure
        for struct_num in structures_to_process:
            zero_indexed = struct_num - 1
            glycan_structure = self.glycan_structures[zero_indexed]
            
            # Generate glycan fragments with peptide (BY or CZ series with -PEP suffix)
            if glycan_fragment_types:
                glycan_frags = self._generate_comprehensive_glycan_fragments(glycan_structure, glycan_fragment_types)
                # Add all returned fragment types with structure tracking
                for key in ['b_ions', 'y_ions', 'by_ions', 'yy_ions', 'byy_ions', 'yyy_ions', 'byyy_ions',
                            'c_ions', 'z_ions', 'cz_ions', 'zz_ions', 'czz_ions', 'zzz_ions', 'czzz_ions', 'bzzz_ions',
                            'a_ions', 'custom_ions']:
                    if key in glycan_frags:
                        for frag in glycan_frags[key]:
                            frag['structure'] = str(struct_num)
                        all_fragments[key].extend(glycan_frags[key])
            
            # Generate peptide fragments (by or cz series)
            if peptide_fragment_types:
                # Use Peptide class to generate fragment series
                pep_obj = Peptide(self.peptide_sequence, use_cam=True)
                pep_frags = pep_obj.generate_all_fragments(fragment_types=peptide_fragment_types)
                
                # Add peptide b and y ions with structure tracking (same for all structures)
                if struct_num == structures_to_process[0]:  # Only add once from first structure
                    for key in ['b_ions', 'y_ions']:
                        if key in pep_frags:
                            for frag in pep_frags[key]:
                                frag['structure'] = 'all'
                            all_fragments[f'peptide_{key}'].extend(pep_frags[key])
                    
                    # Generate peptide neutral losses (H2O, NH3, CO2, CO)
                    neutral_losses = pep_obj._generate_neutral_losses()
                    for frag in neutral_losses:
                        frag['structure'] = 'all'
                    all_fragments['peptide_neutral_losses'].extend(neutral_losses)
                
                # For c-ions: only add those WITHOUT glycosylation site as regular c-ions
                if struct_num == structures_to_process[0]:  # Only add once
                    if 'c_ions' in pep_frags:
                        for c_ion in pep_frags['c_ions']:
                            if c_ion['position'] < self.glycosylation_site:
                                c_ion['structure'] = 'all'
                                all_fragments['peptide_c_ions'].append(c_ion)
                
                # Those WITH glycosylation site are added as c* with glycan (per structure)
                if 'c_ions' in pep_frags:
                    c_with_glyc = self._generate_peptide_c_with_glycan(glycan_structure, pep_frags['c_ions'])
                    for frag in c_with_glyc:
                        frag['structure'] = str(struct_num)
                    all_fragments['peptide_c_ions'].extend(c_with_glyc)
                
                # For z-ions: only add those WITHOUT glycosylation site as regular z-ions
                if struct_num == structures_to_process[0]:  # Only add once
                    if 'z_ions' in pep_frags:
                        threshold = len(self.peptide_sequence) - self.glycosylation_site + 1
                        for z_ion in pep_frags['z_ions']:
                            if z_ion['position'] < threshold:
                                z_ion['structure'] = 'all'
                                all_fragments['peptide_z_ions'].append(z_ion)
                
                # Those WITH glycosylation site are added as z* with glycan (per structure)
                if 'z_ions' in pep_frags:
                    z_with_glyc = self._generate_peptide_z_with_glycan(glycan_structure, pep_frags['z_ions'])
                    for frag in z_with_glyc:
                        frag['structure'] = str(struct_num)
                    all_fragments['peptide_z_ions'].extend(z_with_glyc)
            
            # Generate Y1 ions (peptide + HexNAc1 core) - per structure
            if 'Y1' in fragment_types_upper:
                y1_ions = self._generate_y1_ions(glycan_structure)
                for frag in y1_ions:
                    frag['structure'] = str(struct_num)
                all_fragments['y1_ions'].extend(y1_ions)

            # Generate intact glycopeptide - per structure
            if 'intact' in fragment_types:
                intact = self._generate_intact(glycan_structure)
                intact['structure'] = str(struct_num)
                all_fragments['intact'].append(intact)
        
        return all_fragments
    
    def _generate_comprehensive_glycan_fragments(self, glycan_structure: nx.DiGraph, glycan_fragment_types: List[str]) -> Dict[str, List[Dict[str, Any]]]:
        """
        Generate comprehensive glycan fragments with peptide attached (BY or CZ series with -PEP suffix).
        
        Calls glycan.generate_fragments() to get glycan BY or CZ fragments,
        adding peptide mass and replacing -redend with -PEP suffix.
        
        Args:
            glycan_structure: NetworkX graph of glycan structure
            glycan_fragment_types: List of glycan series to generate: ['BY'] or ['CZ'] or ['BY', 'CZ']
                A-type oxonium and custom diagnostic ions are always generated automatically.
        
        Returns:
            Dict with fragment ion lists (b_ions, y_ions, by_ions, yy_ions, c_ions, z_ions, cz_ions, zz_ions,
            a_ions, custom_ions)
        """
        result = {
            'b_ions': [],
            'y_ions': [],
            'by_ions': [],
            'yy_ions': [],
            'byy_ions': [],
            'yyy_ions': [],
            'byyy_ions': [],
            'c_ions': [],
            'z_ions': [],
            'cz_ions': [],
            'zz_ions': [],
            'czz_ions': [],
            'zzz_ions': [],
            'czzz_ions': [],
            'bzzz_ions': [],
            'a_ions': [],
            'custom_ions': []
        }
        
        if not glycan_fragment_types:
            return result
        
        # Normalize to uppercase (BY, CZ, A)
        glycan_fragment_types = [ft.upper() for ft in glycan_fragment_types]
        
        try:
            # Get glycan fragments using series format (BY or CZ)
            glycan_frags, _ = self.glycan.generate_fragments(glycan_structure, fragment_types=glycan_fragment_types)
            peptide_mass = self.peptide.mass
            
            # Monosaccharide keys to extract (filter out metadata)
            mono_keys = {'HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc', 'H', 'N', 'F', 'S', 'U'}
            
            # Process all fragment types returned by the glycan fragmentor
            # Pure Y, YY, YYY and Z, ZZ, ZZZ have reducing end with peptide
            # BY, CZ, BYY, CZZ, etc. are internal fragments (no reducing end peptide)
            fragment_type_map = {
                # BY Series
                'b_ions': 'B',
                'y_ions': 'Y',
                'by_ions': 'BY',
                'yy_ions': 'YY',
                'byy_ions': 'BYY',
                'yyy_ions': 'YYY',
                'byyy_ions': 'BYYY',
                # CZ Series
                'c_ions': 'C',
                'z_ions': 'Z',
                'cz_ions': 'CZ',
                'zz_ions': 'ZZ',
                'czz_ions': 'CZZ',
                'zzz_ions': 'ZZZ',
                'czzz_ions': 'CZZZ',
                'bzzz_ions': 'BZZZ'
            }
            
            for key, frag_type_label in fragment_type_map.items():
                if key in glycan_frags:
                    for frag_comp_raw in glycan_frags[key]:
                        # Extract only monosaccharide composition (filter out metadata)
                        frag_comp = {k: v for k, v in frag_comp_raw.items() if k in mono_keys}
                        comp_str = self.glycan._format_composition_string(frag_comp)
                        
                        # Only pure Y, YY, YYY and Z, ZZ, ZZZ have reducing end with peptide
                        # All other types (BY, CZ, BYY, CZZ, YYY, etc.) are internal glycan fragments (no peptide)
                        pure_reducing_end = {'Y', 'YY', 'YYY', 'Z', 'ZZ', 'ZZZ'}
                        
                        if frag_type_label in pure_reducing_end:
                            # Add peptide mass but subtract water for glycosidic linkage (consistent with Y1 and Glycopeptide.mass)
                            frag_mass = frag_comp_raw.get('mass', 0) + peptide_mass - WATER_MASS
                            name = f"{frag_type_label}-{comp_str}-PEP"  # Pure Y/Z have reducing end with peptide
                        elif frag_type_label == 'BY':
                            frag_mass = frag_comp_raw.get('mass', 0)  # No peptide
                            name = f"Y-{comp_str}-B"  # Show both cleavage types, no -PEP
                        elif frag_type_label == 'CZ':
                            frag_mass = frag_comp_raw.get('mass', 0)  # No peptide
                            name = f"Z-{comp_str}-C"  # Show both cleavage types, no -PEP
                        elif frag_type_label == 'BYY':
                            frag_mass = frag_comp_raw.get('mass', 0)
                            name = f"YY-{comp_str}-B"  # Show YY-type with B-cleavage
                        elif frag_type_label == 'BYYY':
                            frag_mass = frag_comp_raw.get('mass', 0)
                            name = f"YYY-{comp_str}-B"  # Show YYY-type with B-cleavage
                        elif frag_type_label == 'CZZ':
                            frag_mass = frag_comp_raw.get('mass', 0)
                            name = f"ZZ-{comp_str}-C"  # Show ZZ-type with C-cleavage
                        elif frag_type_label == 'CZZZ':
                            frag_mass = frag_comp_raw.get('mass', 0)
                            name = f"ZZZ-{comp_str}-C"  # Show ZZZ-type with C-cleavage
                        elif frag_type_label == 'BZZZ':
                            frag_mass = frag_comp_raw.get('mass', 0)
                            name = f"ZZZ-{comp_str}-B"  # Show ZZZ-type with B-cleavage
                        else:
                            # YYY, ZZZ (pure reducing end) - add -PEP suffix
                            frag_mass = frag_comp_raw.get('mass', 0)
                            if frag_type_label in {'YYY', 'ZZZ'}:
                                # Add peptide mass but subtract water for glycosidic linkage
                                frag_mass = frag_comp_raw.get('mass', 0) + peptide_mass - WATER_MASS
                                name = f"{frag_type_label}-{comp_str}-PEP"
                            else:
                                name = f"{frag_type_label}-{comp_str}"  # Fallback
                        
                        result[key].append({
                            'name': name,
                            'type': frag_type_label,
                            'composition': comp_str,
                            'glycan_composition': frag_comp,  # Clean composition (only monos)
                            'mass': frag_mass,
                            'peptide_mass': peptide_mass,
                            '_custom_label': name
                        })
            
            # Pass through A-type oxonium ions and custom diagnostic ions directly
            for diag_key in ('a_ions', 'custom_ions'):
                if diag_key in glycan_frags:
                    for frag in glycan_frags[diag_key]:
                        if not isinstance(frag, dict):
                            continue
                        # These ions already have mass/_direct_mz and _custom_label from the fragmentor
                        result[diag_key].append(frag)
        
        except Exception as e:
            # Return empty if fragmentation fails
            pass
        
        return result
    
    def _generate_peptide_y_no_glycan(self) -> List[Dict[str, Any]]:
        """
        Generate peptide y-ions without glycan attachment.
        
        These are C-terminal peptide backbone fragments that span the glycosylation site
        but have no glycan attached. Useful as reference for glycan loss patterns.
        
        Returns:
            List of peptide y-ion dicts without glycan
        """
        y_ions = []
        
        # Get peptide y-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        # Include y-ions that span the glycosylation site
        for y_ion in peptide_frags['y_ions']:
            # Check if this y-ion includes the glycosylation site
            if y_ion['position'] >= len(self.peptide_sequence) - self._glycosylation_index0:
                y_ions.append({
                    'name': f"y-ion[no glycan]-{y_ion['fragment_name']}",
                    'type': 'peptide_y',
                    'sequence': y_ion['sequence'],
                    'position': y_ion['position'],
                    'mass': y_ion['mass'],
                    'composition': 'Peptide',  # No glycan
                    '_custom_label': "y[no glycan]",
                    'has_glycan_site': True,
                    'glycan_attached': False
                })
        
        return y_ions
    
    def _generate_peptide_y_with_glycan(self, glycan_structure: nx.DiGraph) -> List[Dict[str, Any]]:
        """
        Generate peptide y-ions with complete intact glycan attached.
        
        These C-terminal fragments have the full glycan attached at the glycosylation site.
        Represents intact glycopeptide cleavage from C-terminus.
        
        Args:
            glycan_structure: NetworkX graph of glycan structure
        
        Returns:
            List of peptide y-ion with glycan dicts
        """
        y_ions = []
        
        # Get peptide y-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        # Calculate complete glycan mass
        glycan_comp = self.glycan.count_residues(glycan_structure)
        glycan_mass = self.mass_calculator.calculate_glycan_mass(self.glycan_code)
        
        # Include y-ions that span the glycosylation site + full glycan
        for y_ion in peptide_frags['y_ions']:
            if y_ion['position'] >= len(self.peptide_sequence) - self._glycosylation_index0:
                total_mass = y_ion['mass'] + glycan_mass
                comp_str = self.glycan._format_composition_string(glycan_comp)
                y_ions.append({
                    'name': f"y-ion[+glycan]-{y_ion['fragment_name']}",
                    'type': 'peptide_y_intact_glycan',
                    'sequence': y_ion['sequence'],
                    'position': y_ion['position'],
                    'peptide_mass': y_ion['mass'],
                    'glycan_composition': glycan_comp,
                    'glycan_mass': glycan_mass,
                    'mass': total_mass,
                    'composition': comp_str,  # Human-readable composition
                    '_custom_label': "y[+glycan]",
                    'has_glycan_site': True,
                    'glycan_attached': True,
                    'glycan_type': 'complete'
                })
        
        return y_ions
    
    def _generate_peptide_b_ions(self) -> List[Dict[str, Any]]:
        """
        Generate peptide b-ions (N-terminal, never contain glycan).
        
        B-ions from N-terminus cannot carry glycan modifications.
        
        Returns:
            List of b-ion fragment dicts
        """
        b_ions = []
        
        # Get peptide b-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        for b_ion in peptide_frags['b_ions']:
            # Check if this b-ion includes the glycosylation site
            includes_glycan_site = b_ion['position'] > self._glycosylation_index0
            
            b_ions.append({
                'name': b_ion['fragment_name'],
                'type': 'peptide_b',
                'sequence': b_ion['sequence'],
                'position': b_ion['position'],
                'mass': b_ion['mass'],
                'composition': 'Peptide',  # No glycan
                '_custom_label': "b",
                'has_glycan_site': includes_glycan_site,
                'glycan_attached': False  # B-ions never have glycan
            })
        
        return b_ions
    
    def _generate_peptide_c_ions(self) -> List[Dict[str, Any]]:
        """
        Generate peptide c-ions (N-terminal variant, never contain glycan).
        
        C-ions from N-terminus are isobaric with b-ions but have different charge state.
        Cannot carry glycan modifications.
        
        Returns:
            List of c-ion fragment dicts
        """
        c_ions = []
        
        # Get peptide c-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        for c_ion in peptide_frags['c_ions']:
            # Check if this c-ion includes the glycosylation site
            includes_glycan_site = c_ion['position'] > self._glycosylation_index0
            
            c_ions.append({
                'name': c_ion['fragment_name'],
                'type': 'peptide_c',
                'sequence': c_ion['sequence'],
                'position': c_ion['position'],
                'mass': c_ion['mass'],
                'composition': 'Peptide',  # No glycan
                '_custom_label': "c",
                'has_glycan_site': includes_glycan_site,
                'glycan_attached': False  # C-ions never have glycan
            })
        
        return c_ions

    def _generate_peptide_c_with_glycan(
        self,
        glycan_structure: nx.DiGraph,
        c_ions: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Generate peptide c-ions with complete intact glycan attached.
        
        N-terminal fragments that include the glycosylation site can carry
        the intact glycan. These are marked with * in the sequence and name.
        """
        glycan_comp = self.glycan.count_residues(glycan_structure)
        glycan_mass = self.mass_calculator.calculate_glycan_mass(self.glycan_code)
        comp_str = self.glycan._format_composition_string(glycan_comp)
        
        c_with_glycan = []
        for c_ion in c_ions:
            if c_ion['position'] >= self.glycosylation_site:
                seq = c_ion['sequence']
                marker_idx = self.glycosylation_site - 1
                if 0 <= marker_idx < len(seq):
                    seq = f"{seq[:marker_idx+1]}*{seq[marker_idx+1:]}"
                frag_name = f"c{c_ion['position']}*"
                total_mass = c_ion['mass'] + glycan_mass - WATER_MASS
                c_with_glycan.append({
                    'name': frag_name,
                    'fragment_name': frag_name,
                    'fragment_type': 'c',
                    'type': 'c',
                    'sequence': seq,
                    'position': c_ion['position'],
                    'peptide_mass': c_ion['mass'],
                    'glycan_composition': glycan_comp,
                    'glycan_mass': glycan_mass,
                    'mass': total_mass,
                    'composition': comp_str,
                    '_custom_label': "c*",
                    'has_glycan_site': True,
                    'glycan_attached': True,
                    'glycan_type': 'complete',
                    'glycosylation_site': self.glycosylation_site
                })
        

            def format_structure_tree(self, structure: nx.DiGraph, root: int = 1) -> str:
                """
                Format a glycan structure tree for display.

                Args:
                    structure: NetworkX DiGraph representing the glycan structure
                    root: Root node id for display (default 1)

                Returns:
                    String representation of the structure tree
                """
                return self.glycan.format_structure_tree(structure, root=root)
        return c_with_glycan
    
    def _generate_peptide_z_no_glycan(self) -> List[Dict[str, Any]]:
        """
        Generate peptide z-ions without glycan attachment.
        
        Z-ions are C-terminal fragments (like y-ions) that span the glycosylation site
        but have no glycan attached. Useful as reference for glycan loss patterns.
        
        Returns:
            List of peptide z-ion dicts without glycan
        """
        z_ions = []
        
        # Get peptide z-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        # Include z-ions that span the glycosylation site
        for z_ion in peptide_frags.get('z_ions', []):
            # Check if this z-ion includes the glycosylation site
            if z_ion['position'] >= len(self.peptide_sequence) - self._glycosylation_index0:
                z_ions.append({
                    'name': f"z-ion[no glycan]-{z_ion['fragment_name']}",
                    'type': 'peptide_z',
                    'sequence': z_ion['sequence'],
                    'position': z_ion['position'],
                    'mass': z_ion['mass'],
                    'composition': 'Peptide',  # No glycan
                    '_custom_label': "z[no glycan]",
                    'has_glycan_site': True,
                    'glycan_attached': False
                })
        
        return z_ions
    
    def _generate_peptide_z_with_glycan(
        self,
        glycan_structure: nx.DiGraph,
        z_ions: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Generate peptide z-ions with complete intact glycan attached.
        
        C-terminal fragments that include the glycosylation site can carry
        the intact glycan. These are marked with * in the sequence and name.
        
        Args:
            glycan_structure: NetworkX graph of glycan structure
            z_ions: List of peptide z-ions to process
        
        Returns:
            List of peptide z-ion with glycan dicts
        """
        z_with_glycan = []
        
        # Calculate complete glycan mass
        glycan_comp = self.glycan.count_residues(glycan_structure)
        glycan_mass = self.mass_calculator.calculate_glycan_mass(self.glycan_code)
        comp_str = self.glycan._format_composition_string(glycan_comp)
        
        # z-ions include glycosylation site if position >= (seq_length - glyc_site + 1)
        threshold = len(self.peptide_sequence) - self.glycosylation_site + 1
        
        # Include z-ions that span the glycosylation site + full glycan
        for z_ion in z_ions:
            if z_ion['position'] >= threshold:
                # Remove water for glycosidic linkage
                total_mass = z_ion['mass'] + glycan_mass - WATER_MASS
                frag_name = f"z{z_ion['position']}*"
                seq = z_ion['sequence']
                start_idx = len(self.peptide_sequence) - z_ion['position']
                marker_idx = self.glycosylation_site - start_idx - 1
                if 0 <= marker_idx < len(seq):
                    seq = f"{seq[:marker_idx+1]}*{seq[marker_idx+1:]}"
                z_with_glycan.append({
                    'name': frag_name,
                    'fragment_name': frag_name,
                    'fragment_type': 'z',
                    'type': 'z',
                    'sequence': seq,
                    'position': z_ion['position'],
                    'peptide_mass': z_ion['mass'],
                    'glycan_composition': glycan_comp,
                    'glycan_mass': glycan_mass,
                    'mass': total_mass,
                    'composition': comp_str,  # Human-readable composition
                    '_custom_label': "z*",
                    'has_glycan_site': True,
                    'glycan_attached': True,
                    'glycan_type': 'complete',
                    'glycosylation_site': self.glycosylation_site
                })
        
        return z_with_glycan
    
    def _generate_y_glycan_fragments(
        self,
        glycan_structure: nx.DiGraph
    ) -> List[Dict[str, Any]]:
        """
        Generate Y-glycan ions (peptide y-ions + glycan fragments).
        
        These combine peptide y-ions with B-type glycan fragments (oxonium ions).
        This creates Y + glycan fragment ions.
        
        Args:
            glycan_structure: NetworkX graph of glycan structure
        
        Returns:
            List of Y-glycan fragment dicts
        """
        y_glycan_ions = []
        
        # Get peptide y-ions
        peptide_frags = self.peptide.generate_all_fragments()
        
        # Get glycan B-ion fragments
        glycan_frags, _ = self.glycan.generate_fragments(
            glycan_structure,
            fragment_types=['BY']
        )
        
        # Combine y-ions with glycan Y-ions (which contain reducing end)
        for y_ion in peptide_frags['y_ions']:
            if y_ion['position'] >= len(self.peptide_sequence) - self._glycosylation_index0:
                # This y-ion includes the glycosylation site
                
                # Add Y + glycan Y-ions (glycan fragments with reducing end)
                for glycan_y in glycan_frags.get('y_ions', []):
                    # Skip empty fragments
                    if sum(glycan_y.values()) == 0:
                        continue
                    
                    # Calculate glycan fragment mass
                    glycan_frag_str = self.glycan._format_fragment_string(glycan_y, 'Y')
                    glycan_y_mass = self.mass_calculator.calculate_composition_mass(glycan_y)
                    # Peptide y-ion + glycan Y-ion with water loss at glycosidic bond
                    total_mass = y_ion['mass'] + glycan_y_mass - WATER_MASS
                    comp_str = self.glycan._format_composition_string(glycan_y)
                    
                    y_glycan_ions.append({
                        'name': f"{y_ion['fragment_name']}+{glycan_frag_str}",
                        'type': 'Y-glycan',
                        'peptide_name': y_ion['fragment_name'],
                        'sequence': y_ion['sequence'],
                        'position': y_ion['position'],
                        'peptide_mass': y_ion['mass'],
                        'glycan_composition': glycan_y,
                        'glycan_fragment': glycan_frag_str,
                        'glycan_mass': glycan_y_mass,
                        'mass': total_mass,  # Added total mass
                        'composition': comp_str,  # Human-readable composition
                        '_custom_label': f"Y-{comp_str}",
                        'has_glycan_site': True,
                        'glycan_attached': True,
                        'glycan_type': 'fragment'
                    })
        
        return y_glycan_ions
    
    def _generate_intact(self, glycan_structure: nx.DiGraph) -> Dict[str, Any]:
        """
        Generate intact glycopeptide ion.
        
        Args:
            glycan_structure: NetworkX graph of glycan structure
        
        Returns:
            Dict representing intact glycopeptide
        """
        # Calculate masses
        peptide_mass = self.peptide.mass
        glycan_comp = self.glycan.count_residues(glycan_structure)
        glycan_mass = self.mass_calculator.calculate_glycan_mass(self.glycan_code)
        
        # Intact glycopeptide: subtract water for glycosidic linkage (consistent with Glycopeptide.mass property)
        total_mass = peptide_mass + glycan_mass - WATER_MASS
        comp_str = self.glycan._format_composition_string(glycan_comp)
        
        return {
            'name': 'Intact',
            'type': 'intact',
            'sequence': self.peptide_sequence,
            'peptide_mass': peptide_mass,
            'glycan_mass': glycan_mass,
            'mass': total_mass,  # Added total mass
            'composition': comp_str,  # Human-readable composition
            '_custom_label': "Intact",
            'glycan_composition': glycan_comp,
            'glycan_mass': glycan_mass,
            'total_mass': peptide_mass + glycan_mass,
            'glycosylation_site': self.glycosylation_site,
            'glycan_code': self.glycan_code
        }

    def _generate_y1_ions(self, glycan_structure: nx.DiGraph) -> List[Dict[str, Any]]:
        """
        Generate Y1 ions: full peptide + core HexNAc1 (and core Fuc if present).

        Y1 is distinct from peptide y1; it represents the intact peptide with the
        core HexNAc attached at the glycosylation site.
        """
        glycan_comp = self.glycan.count_residues(glycan_structure)
        core_comp: Dict[str, int] = {'HexNAc': 1}

        # Add core fucose if present in the full glycan composition
        if glycan_comp.get('Fuc', 0) > 0:
            core_comp['Fuc'] = 1

        peptide_mass = self.peptide.mass
        glycan_mass = 0.0
        for mono, count in core_comp.items():
            glycan_mass += self.mass_calculator.MONO_MASSES.get(mono, 0.0) * count

        # Remove water for the glycosidic linkage in glycopeptide mode
        total_mass = peptide_mass + glycan_mass - WATER_MASS

        comp_str = self.glycan._format_composition_string(core_comp)
        name = "Y1"

        return [
            {
                'name': name,
                'type': 'Y1',
                'sequence': self.peptide_sequence,
                'peptide_mass': peptide_mass,
                'glycan_composition': core_comp,
                'glycan_mass': glycan_mass,
                'mass': total_mass,
                'total_mass': total_mass,
                'composition': comp_str,
                '_custom_label': name,
                'glycosylation_site': self.glycosylation_site,
                'glycan_attached': True,
                'glycan_type': 'core'
            }
        ]
    
    def __repr__(self) -> str:
        """String representation of glycopeptide."""
        glycan_aa = self.peptide_sequence[self._glycosylation_index0]
        return (f"Glycopeptide(sequence='{self.peptide_sequence}', "
                f"glycan='{self.glycan_code}' at {glycan_aa}{self.glycosylation_site})")
    
    def visualize_structure(self, structure_number: int = 1, **kwargs) -> Optional[Any]:
        """
        Visualize predicted glycan structure.
        
        Directly access glycan structure visualization without rewriting code.
        Only visualizes the glycan part (no peptide visualization).
        
        Args:
            structure_number: Which structure to visualize (1-indexed, default 1 = first structure)
            **kwargs: Additional arguments passed to GlycanVisualizer (e.g., title, figsize, output_path)
        
        Returns:
            Path to output file (if output_path provided) or None
        
        Example:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N')
            >>> gp.visualize_structure(structure_number=1)  # Visualize first predicted structure
            >>> gp.visualize_structure(structure_number=2, output_path='structure_2.png')
        
        Note:
            This is a convenience wrapper around GlycanVisualizer.visualize_structure()
            Only the glycan part is visualized; peptide is not included in visualization.
        """
        if not self.glycan_structures:
            print(f"[ERROR] No glycan structures available for {self.glycan_code}")
            return None
        
        return GlycanVisualizer.visualize_structure(
            self.glycan_structures,
            structure_number=structure_number,
            **kwargs
        )
    
    def visualize_structures(self, 
                            structures_to_draw: Union[str, int, List[int]] = "all",
                            figsize: tuple = (5, 4),
                            save_images: bool = False,
                            show_plots: bool = True,
                            output_dir: Optional[str] = None) -> List[str]:
        """
        Visualize multiple glycan structures with flexible selection options.
        
        This is a convenience method for visualizing all or selected structures
        from the predicted glycan ensemble.
        
        Args:
            structures_to_draw: Which structures to visualize
                - 'all': Visualize all predicted structures (default)
                - int: Single structure by number, e.g., 1 (visualizes structure 1)
                - list/tuple: Multiple structures, e.g., [1, 2, 5]
                - range: Range of structures, e.g., range(1, 4) for 1,2,3
                
            figsize: Figure size (width, height) in inches. Default: (5, 4)
            save_images: Save PNG files for each structure. Default: False
            show_plots: Display plots in matplotlib window. Default: True
            output_dir: Directory to save images (if save_images=True).
                       Default: current working directory
            
        Returns:
            List of output file paths (if save_images=True), otherwise empty list
        
        Examples:
            >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N')
            >>> # Visualize all structures
            >>> gp.visualize_structures()
            
            >>> # Visualize only structure 1
            >>> gp.visualize_structures(structures_to_draw=1, show_plots=False, save_images=True)
            
            >>> # Visualize structures 1-5
            >>> gp.visualize_structures(structures_to_draw=range(1, 6), save_images=True)
            
            >>> # Visualize specific structures with custom output directory
            >>> gp.visualize_structures([1, 2, 10], save_images=True, output_dir='/tmp/glycans')
        """
        import os
        from glycofrag.io.visualizer import GlycanVisualizer
        
        if not self.glycan_structures:
            print(f"[ERROR] No glycan structures available for {self.glycan_code}")
            return []
        
        total_structures = len(self.glycan_structures)
        output_files = []
        
        # Determine which structures to visualize
        if structures_to_draw == "all":
            structure_numbers = list(range(1, total_structures + 1))
        elif isinstance(structures_to_draw, int):
            structure_numbers = [structures_to_draw]
        elif isinstance(structures_to_draw, (list, tuple, range)):
            structure_numbers = list(structures_to_draw)
        else:
            print(f"[ERROR] Invalid structures_to_draw: {structures_to_draw}")
            return []
        
        # Validate structure numbers
        valid_numbers = [n for n in structure_numbers if 1 <= n <= total_structures]
        if not valid_numbers:
            print(f"[ERROR] No valid structure numbers (have {total_structures} available)")
            return []
        
        if len(valid_numbers) < len(structure_numbers):
            invalid = set(structure_numbers) - set(valid_numbers)
            print(f"[WARN] Skipping invalid structure numbers: {sorted(invalid)}")
        
        # Set output directory
        if output_dir is None and save_images:
            output_dir = os.getcwd()
        
        print(f"\n{'='*80}")
        print(f"VISUALIZING GLYCAN STRUCTURES - {self.glycan_code}")
        print(f"{'='*80}")
        print(f"Peptide sequence:   {self.peptide_sequence}")
        print(f"Glycosylation site: {self.glycosylation_site}")
        print(f"Glycan type:        {self.glycan_type}-glycan")
        print(f"Total structures:   {total_structures}")
        print(f"Structures to draw: {len(valid_numbers)} ({', '.join(map(str, valid_numbers[:10]))}{'...' if len(valid_numbers) > 10 else ''})")
        print(f"Save images:        {save_images}")
        if save_images:
            print(f"Output directory:   {output_dir}")
        print(f"{'='*80}\n")
        
        # Visualize each structure
        for struct_num in valid_numbers:
            try:
                # Create output filename if saving
                output_path = None
                if save_images and output_dir is not None:
                    filename = f"{self.glycan_code}_structure_{struct_num:03d}.png"
                    output_path = os.path.join(output_dir, filename)
                
                # Create title
                title = f"{self.glycan_code} - Structure {struct_num}/{total_structures}"
                if self.peptide_sequence:
                    title += f" ({self.peptide_sequence})"
                
                if len(valid_numbers) <= 5 or struct_num <= 3:
                    print(f"  Visualizing structure {struct_num:3d}/{total_structures}...", end=" ", flush=True)
                
                # Call the visualize method
                result = self.visualize_structure(
                    structure_number=struct_num,
                    title=title,
                    figsize=figsize,
                    show=show_plots,
                    output_path=output_path
                )
                
                if result:
                    output_files.append(result)
                    if len(valid_numbers) <= 5 or struct_num <= 3:
                        print(f"[OK] ({os.path.basename(result)})")
                elif len(valid_numbers) <= 5 or struct_num <= 3:
                    print("[OK]")
            
            except Exception as e:
                print(f"\n  [ERROR] Failed to visualize structure {struct_num}:")
                print(f"    {type(e).__name__}: {e}")
        
        if len(valid_numbers) > 3:
            print(f"  ... ({len(valid_numbers) - 3} more structures processed)")
        
        if output_files:
            print(f"\n{'='*80}")
            print(f"[SUCCESS] Successfully visualized {len(valid_numbers)} structure(s)")
            if save_images and output_files:
                print(f"[SUCCESS] Saved {len(output_files)} structure image(s) to {output_dir}")
            print(f"{'='*80}\n")
        
        return output_files
