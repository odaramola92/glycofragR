# glycofrag/glycan/facade.py
"""Glycan facade - high-level API delegating to specialized modules."""

from typing import List, Optional, Dict, Any, Tuple, Set, Union
import networkx as nx

from glycofrag.composition import GlycanComposition
from glycofrag.structure import StructureBuilder, StructureClassifier, StructurePredictor
from glycofrag.fragmentation import GlycanFragmentor
from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.modifications import get_modification_type


class Glycan:
    """
    High-level Glycan API - facade delegating to specialized modules.
    
    This class provides a clean, unified interface for glycan analysis:
    - Composition parsing
    - Structure prediction with biosynthetic rules
    - Structure classification and labeling
    - Fragment ion generation
    
    Internally delegates to:
    - GlycanComposition: Parsing and composition data
    - StructurePredictor: Structure prediction algorithms
    - StructureClassifier: Classification and labeling
    - GlycanFragmentor: Fragment ion generation
    
    This maintains backward compatibility with the original Glycan class
    while providing a cleaner, more maintainable architecture.

    Public API:
        - predict_structures(...)
        - generate_fragments(...)
        - classify_structure(...)
        - count_residues(...)
    """
    
    def __init__(self, code: str, glycan_type: str = "N", max_structures: int = 100,
                 modification_type: Union[str, int] = 0, isomer_sensitive: bool = False,
                 preferred_core: Optional[int] = None, _allow_glycopeptide: bool = False,
                 custom_reducing_end_mass: Optional[float] = None,
                 custom_mono_masses: Optional[Dict[str, float]] = None):
        """
        Initialize Glycan from composition code.
        
        Args:
            code: 4-digit composition code (e.g., '4501')
                Format: HNFS where H=HexNAc, N=Hex, F=Fuc, S=NeuAc
                For 5+ digits: HNFSG where G=NeuGc
                Examples:
                  - '4501': 4 HexNAc, 5 Hex, 0 Fuc, 1 NeuAc
                  - '5411': 5 HexNAc, 4 Hex, 1 Fuc, 1 NeuAc
                  - '2100': 2 HexNAc, 1 Hex, 0 Fuc, 0 NeuAc (O-glycan)
            
            glycan_type: 'N' for N-glycan or 'O' for O-glycan
            max_structures: Maximum structures to predict (default: 100)
            
            modification_type: Reducing end modification. Controls how the
                reducing end (node 1) mass is calculated.
                
                **String options** (recommended):
                  - 'free'              → Free reducing end (+18.01 Da, H₂O)
                  - 'reduced'           → Reduced/alditol end (+20.03 Da)
                  - 'permethylated_free' or 'pm_free'
                                        → Permethylated free reducing end
                  - 'permethylated_reduced' or 'pm_reduced'
                                        → Permethylated reduced end
                  - 'two_ab' or '2ab'   → 2-AB (2-aminobenzamide) labeled
                  - 'two_ab_permethylated' or '2ab_permethylated'
                                        → 2-AB labeled + permethylated
                  - 'custom'            → User-defined mass (see custom_reducing_end_mass
                                          and custom_mono_masses parameters below)
                
                **Integer options** (for advanced users):
                  - 0: Free reducing end (same as 'free')
                  - 1: Reduced end (same as 'reduced')
                  - 2: Permethylated free (same as 'pm_free')
                  - 3: Permethylated reduced (same as 'pm_reduced')
                  - 4: 2-AB labeled (same as '2ab')
                  - 5: 2-AB + permethylated (same as '2ab_permethylated')
                  - 7: Custom (same as 'custom')
                
                Default: 'free' (0) — standard for standalone glycan analysis.
                
                Note: Type 6 ('glycopeptide') is reserved for the Glycopeptide class
                and cannot be used here. Use Glycopeptide for peptide-attached glycans.
            
            isomer_sensitive: If True, treats mirror images as distinct structures
            preferred_core: For O-glycans only — preferred core type (0-8) for visualization.
                          Does NOT affect prediction or fragmentation, only visualization.
                          Example: preferred_core=5 will visualize as Core 5 if possible.
            
            custom_reducing_end_mass: Mass (Da) to add at the reducing end.
                Only used when modification_type='custom' (7).
                This mass is added to Y-type ions and intact glycan mass.
                Example: 18.010564684 replicates a free reducing end.
            
            custom_mono_masses: Dict of additional mass (Da) to add per monosaccharide.
                Only used when modification_type='custom' (7).
                Keys are monosaccharide names: 'HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc'.
                Values are masses to ADD to each residue of that type.
                Example: {'HexNAc': 14.016, 'Hex': 14.016} adds 14.016 Da to every
                HexNAc and Hex residue (e.g., for a custom derivatization).
        
        Examples:
            >>> # Standard free reducing end (default)
            >>> glycan = Glycan('4501', glycan_type='N')
            
            >>> # Permethylated sample
            >>> glycan = Glycan('4501', modification_type='pm_reduced')
            
            >>> # 2-AB labeled
            >>> glycan = Glycan('4501', modification_type='2ab')
            
            >>> # Custom: add 120.0 Da label at reducing end
            >>> glycan = Glycan('4501', modification_type='custom',
            ...                 custom_reducing_end_mass=120.0)
            
            >>> # Custom: add mass to specific monosaccharides
            >>> glycan = Glycan('4501', modification_type='custom',
            ...                 custom_mono_masses={'NeuAc': 28.03})
            
            >>> # Custom: reducing end label + monosaccharide modification
            >>> glycan = Glycan('4501', modification_type='custom',
            ...                 custom_reducing_end_mass=120.0,
            ...                 custom_mono_masses={'HexNAc': 14.016})
            
            >>> # O-glycan with core preference
            >>> o_glycan = Glycan('2100', glycan_type='O', preferred_core=5)
        """
        self.code = code
        self.glycan_type = glycan_type.upper()
        self.max_structures = max_structures
        self.custom_reducing_end_mass = custom_reducing_end_mass
        self.custom_mono_masses = custom_mono_masses or {}
        
        # Normalize modification_type (accept string or int)
        self.modification_type = get_modification_type(modification_type)
        self._allow_glycopeptide = _allow_glycopeptide

        if self.modification_type == get_modification_type('glycopeptide') and not _allow_glycopeptide:
            raise ValueError(
                "Glycan modification_type 'glycopeptide' is reserved for Glycopeptide. "
                "Use Glycopeptide for attached glycans."
            )
        
        # Validate custom params are only used with type 7 ('custom')
        # Note: custom_mono_masses is allowed with any modification type since it's independent
        if self.modification_type != 7:
            if custom_reducing_end_mass is not None:
                raise ValueError(
                    "custom_reducing_end_mass can only be used with modification_type='custom' (7). "
                    f"Got modification_type={modification_type}"
                )
        
        self.isomer_sensitive = isomer_sensitive
        self.preferred_core = preferred_core
        
        # Parse composition
        self.composition = GlycanComposition.from_code(code)
        
        # Create specialized components
        self.predictor = StructurePredictor(
            self.composition,
            glycan_type=self.glycan_type,
            max_structures=max_structures,
            modification_type=modification_type,
            isomer_sensitive=isomer_sensitive,
            preferred_core=preferred_core
        )
        self.classifier = StructureClassifier()
        self.fragmentor = GlycanFragmentor(
            glycan_code=code,
            glycan_type=self.glycan_type,
            modification_type=modification_type,
            custom_reducing_end_mass=custom_reducing_end_mass,
            custom_mono_masses=custom_mono_masses
        )
        
        # Results cache
        self.possible_structures: List[nx.DiGraph] = []
        self.structure_fingerprints: Set[str] = set()
    
    def predict_structures(self) -> List[nx.DiGraph]:
        """
        Generate all possible complete structures for the given glycan code.
        
        This method delegates to StructurePredictor for all structure generation logic.
        
        For O-glycans, if preferred_core is set and no structures match, raises an error
        with available core options.
        
        Returns:
            List of NetworkX DiGraphs representing possible glycan structures
            Empty list if composition is invalid or no structures can be generated
        
        Raises:
            ValueError: If preferred_core is specified but not predicted
        
        Example:
            >>> from glycofrag import Glycan
            >>> glycan = Glycan('4501', glycan_type='N', max_structures=1000)
            >>> structures = glycan.predict_structures()
            >>> print(f"Generated {len(structures)} structures")
            
            >>> # O-glycan with core preference
            >>> o_glycan = Glycan('2100', glycan_type='O', preferred_core=5)
            >>> structures = o_glycan.predict_structures()
        """
        # Delegate to predictor and cache results locally
        self.possible_structures = self.predictor.predict_structures()
        self.structure_fingerprints = self.predictor.structure_fingerprints

        # Attach 1-based structure numbers for table labeling
        for idx, structure in enumerate(self.possible_structures, start=1):
            structure.graph.setdefault('structure_number', idx)
        
        # For O-glycans with preferred_core, validate that structures were predicted
        if self.glycan_type == 'O' and self.preferred_core is not None:
            if not self.possible_structures:
                # Get available cores from builder to show user what's possible
                core_graphs = self.predictor.builder.build_o_glycan_cores(
                    self.composition.hexnac, self.composition.hex
                )
                available_cores = []
                for core_data in core_graphs:
                    cores = core_data.get('possible_cores', [])
                    available_cores.extend(cores)
                
                if available_cores:
                    raise ValueError(
                        f"Core {self.preferred_core} was not predicted for this composition.\n"
                        f"Available cores for composition {self.code}: {sorted(set(available_cores))}\n"
                        f"Please rerun with preferred_core set to one of the available cores."
                    )
        
        return self.possible_structures
    
    def get_predicted_cores(self) -> List[Dict[str, Any]]:
        """
        Get information about predicted O-glycan cores.
        
        Returns list of core information for O-glycans, empty list for N-glycans.
        Each dict contains:
        - 'possible_cores': List of core numbers this structure represents
        - 'default_core': Default core number used
        - 'selected_core': Currently selected core for visualization
        - 'core_name': Human-readable core name
        
        Example:
            >>> o_glycan = Glycan('2100', glycan_type='O')
            >>> structures = o_glycan.predict_structures()
            >>> cores = o_glycan.get_predicted_cores()
            >>> for core_info in cores:
            ...     print(f"{core_info['core_name']}: Possible cores {core_info['possible_cores']}")
        """
        if self.glycan_type != 'O' or not self.possible_structures:
            return []
        
        cores_info = []
        seen_cores = set()
        
        for structure in self.possible_structures:
            possible_cores = tuple(structure.graph.get('possible_cores', []))
            if possible_cores and possible_cores not in seen_cores:
                cores_info.append({
                    'possible_cores': list(possible_cores),
                    'default_core': structure.graph.get('default_core'),
                    'selected_core': structure.graph.get('selected_core'),
                    'core_name': structure.graph.get('core_name', 'Unknown')
                })
                seen_cores.add(possible_cores)
        
        return cores_info

    def generate_fragments(
        self,
        structure: Optional[nx.DiGraph] = None,
        fragment_types: Optional[List[str]] = None,
        modification_type: Optional[Union[str, int]] = None,
        charges: Optional[List[int]] = None
    ) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, str]]]:
        """
        Generate glycan fragment ions from a glycan structure with m/z values.
        
        If structure is not provided, automatically predicts structures and uses
        the first one.
        
        Fragment Series (use uppercase for selection):
            - BY: Generates b_ions, y_ions, by_ions, yy_ions, yyy_ions,
                  byy_ions, byyy_ions
            - CZ: Generates c_ions, z_ions, cz_ions, zz_ions, zzz_ions,
                  czz_ions, czzz_ions, bzz_ions, bzzz_ions
            - A:  Oxonium ions (monosaccharide diagnostic fragments) — always
                  generated automatically as add-ons regardless of selection.
        
        The reducing end (node 1) serves as a PLACEHOLDER for glycan-only
        analysis (use modification_type 0–5 or 7/'custom'). For glycan
        attached to a peptide, use the ``Glycopeptide`` class instead.
        
        Args:
            structure: NetworkX DiGraph representing the glycan structure.
                If ``None``, automatically predicts structures and uses the
                first one.
            fragment_types: List of fragment series to generate.
                Options: ``['BY']``, ``['CZ']``, or ``['BY', 'CZ']``.
                Default: ``['BY']`` (BY series only).
                A-type oxonium and custom diagnostic ions are always generated
                automatically as add-ons.
            modification_type: Reducing end modification type override. If
                ``None``, uses the instance ``modification_type`` set at
                ``Glycan`` initialization.
                
                Accepted values (string or int):
                    - 0 / ``'free'``:  Free reducing end (H₂O)
                    - 1 / ``'reduced'``:  Reduced end (alditol)
                    - 2 / ``'permethylated'``:  Permethylated free end
                    - 3 / ``'permethylated_reduced'``:  Permethylated reduced
                    - 4 / ``'2ab'``:  2-AB labeled
                    - 5 / ``'2ab_permethylated'``:  2-AB labeled + permethylated
                    - 7 / ``'custom'``:  Custom (uses ``custom_reducing_end_mass``
                      and ``custom_mono_masses`` set at initialization)
            charges: List of charge states to calculate m/z for.
                Default: ``[1, 2, 3]``.
                Examples: ``[1]``, ``[1, 2]``, ``[1, 2, 3, 4]``
        
        Returns:
            Tuple of ``(fragments_dict, cleavage_info_dict)``.
            
            ``fragments_dict`` maps ion type names (e.g. ``'b_ions'``,
            ``'y_ions'``) to lists of fragment dicts, each containing:
                - ``'mass'``:  Neutral mass in Da
                - ``'mz_charge_X'``:  m/z value for charge state *X*
                  (e.g. ``'mz_charge_1'``, ``'mz_charge_2'``)
                - ``'composition'``:  Fragment composition string
            
            ``cleavage_info_dict`` maps fragment labels to cleavage
            descriptions.
        
        Example:
            >>> # Generate BY series (default)
            >>> glycan = Glycan('4501', glycan_type='N')
            >>> fragments, info = glycan.generate_fragments()
            >>> print(f"B-ions: {len(fragments.get('b_ions', []))}")
            >>>
            >>> # Generate BY and CZ series
            >>> fragments, info = glycan.generate_fragments(
            ...     fragment_types=['BY', 'CZ']
            ... )
            >>>
            >>> # Custom modification with specific charges
            >>> glycan = Glycan('4501', modification_type='custom',
            ...     custom_reducing_end_mass=120.0,
            ...     custom_mono_masses={'HexNAc': 5.0})
            >>> fragments, info = glycan.generate_fragments(charges=[1, 2])
        """
        # Auto-predict structures if not provided
        if structure is None:
            if not self.possible_structures:
                self.predict_structures()
            if self.possible_structures:
                structure = self.possible_structures[0]
            else:
                raise ValueError("No structures could be predicted for the given glycan code")
        
        if modification_type is not None:
            mod_type = get_modification_type(modification_type)
            if mod_type == get_modification_type('glycopeptide') and not self._allow_glycopeptide:
                raise ValueError(
                    "Glycan modification_type 'glycopeptide' is reserved for Glycopeptide. "
                    "Use Glycopeptide for attached glycans."
                )
            modification_type = mod_type

        # Default to BY series
        if fragment_types is None:
            fragment_types = ['BY']
        
        # Normalize to uppercase and ensure valid types
        fragment_types = [ft.upper() for ft in fragment_types]
        
        # Delegate fragment generation to fragmentor
        fragments, cleavage_info = self.fragmentor.generate_fragments(
            structure, 
            fragment_types=fragment_types, 
            modification_type=modification_type,
            charges=charges
        )
        
        return fragments, cleavage_info
    
    def classify_structure(self, graph: nx.DiGraph) -> str:
        """
        Classify the glycan structure based on composition and type.
        
        Args:
            graph: NetworkX DiGraph representing the glycan structure
        
        Returns:
            String classification (e.g., 'High Mannose', 'Complex', 'O-GalNAc')
        """
        return self.classifier.classify_structure(
            graph, 
            self.glycan_type,
            self.composition.hexnac,
            self.composition.hex
        )
    
    def _count_residues(self, graph: nx.DiGraph) -> Dict[str, int]:
        """
        Count the number of each monosaccharide type in a graph.
        
        Args:
            graph: NetworkX DiGraph representing a glycan structure
        
        Returns:
            Dictionary with counts for HexNAc, Hex, Fuc, NeuAc
            Note: Man and Gal are counted as Hex
        """
        counts = {'HexNAc': 0, 'Hex': 0, 'Fuc': 0, 'NeuAc': 0}
        for node in graph.nodes():
            node_type = graph.nodes[node]['type']
            if node_type in counts:
                counts[node_type] += 1
            elif node_type in ['Man', 'Gal']:
                # Man and Gal are both Hex types
                counts['Hex'] += 1
        return counts
    
    def count_residues(self, graph: nx.DiGraph) -> Dict[str, int]:
        """
        Public API - Count the number of each monosaccharide type in a graph.
        
        This is the public version of _count_residues() used by glycopeptide analysis.
        
        Args:
            graph: NetworkX DiGraph representing a glycan structure
        
        Returns:
            Dictionary with counts for HexNAc, Hex, Fuc, NeuAc
        """
        return self._count_residues(graph)
    
    def _format_fragment_string(self, composition: Dict[str, int], ion_type: str, is_y_ion: bool = False) -> str:
        """
        Format a fragment composition as a string.
        
        This delegates to the fragmentor for consistent fragment string formatting.
        
        Args:
            composition: Dict of monosaccharide counts
            ion_type: Type of ion ('B', 'Y', 'BY', 'YY', 'YYY')
            is_y_ion: If True, append modification-specific suffix for Y-type ions
        
        Returns:
            Formatted string like "HexNAc2Hex3-B" or "HexNAc-Redend" for Y-type ions
        """
        return self.fragmentor._format_fragment_string(composition, ion_type, is_y_ion)
    
    def _format_composition_string(self, composition: Dict[str, int]) -> str:
        """
        Format a glycan composition as a compact human-readable string.
        
        Args:
            composition: Dict of monosaccharide counts (e.g., {'H': 5, 'N': 4, 'S': 2})
        
        Returns:
            Formatted string like "H5N4S2" or "H3N2" (standard IUPACfragment notation shorthand)
        """
        if not composition:
            return ""
        
        # Order: H (Hex), N (HexNAc), F (Fuc), S (Sia/NeuAc), P (Phospho), etc.
        order = ['H', 'N', 'F', 'S', 'P', 'U']  # U for unknowns
        parts = []
        
        for key in order:
            if key in composition and composition[key] > 0:
                parts.append(f"{key}{composition[key]}")
        
        # Add any other keys not in the standard order
        for key, count in sorted(composition.items()):
            if key not in order and count > 0:
                parts.append(f"{key}{count}")
        
        return "".join(parts)
    
    def _count_branches(self, structure: nx.DiGraph) -> int:
        """
        Count the number of branches (arms) in a glycan structure.
        
        For N-glycans:
        - The core is typically at node 3 (Hex)
        - Branches extend from nodes 4 and 5 (core mannoses)
        - Each arm from node 4 or 5 counts as 1 branch
        
        For O-glycans:
        - Count based on branching from the core
        
        Args:
            structure: The glycan structure graph
            
        Returns:
            Number of branches/arms in the structure
        """
        branch_count = 0
        
        if self.glycan_type == 'N':
            for core_node in [4, 5]:
                if core_node in structure.nodes():
                    successors = list(structure.successors(core_node))
                    if len(successors) > 0:
                        branch_count += 1
        else:
            for node in structure.nodes():
                out_degree = structure.out_degree(node)
                if out_degree > 1:
                    branch_count += (out_degree - 1)
            branch_count = max(1, branch_count)
        
        return branch_count
  
    def format_structure_tree(self, structure: nx.DiGraph, root: int = 1) -> str:
        """
        Format glycan structure as ASCII tree.
        
        Args:
            structure: NetworkX DiGraph
            root: Root node ID (default 1)
        
        Returns:
            ASCII tree string
        """
        lines: List[str] = []
        
        def node_label(node_id: int) -> str:
            node_type = structure.nodes[node_id]['type']
            label = structure.nodes[node_id].get('label', node_type)
            return f"{label}({node_id})"
        
        def walk(node_id: int, prefix: str = "", is_last: bool = True):
            current = "└── " if is_last else "├── "
            lines.append(prefix + current + node_label(node_id))
            children = list(structure.successors(node_id))
            for i, child in enumerate(children):
                extension = "    " if is_last else "│   "
                walk(child, prefix + extension, i == len(children) - 1)
        
        if root in structure.nodes:
            lines.append(node_label(root))
            for child in structure.successors(root):
                walk(child, "", True)
        return "\n".join(lines)
    
    def __repr__(self):
        """String representation."""
        return (f"Glycan(code='{self.code}', type='{self.glycan_type}', "
                f"structures={len(self.possible_structures)})")
