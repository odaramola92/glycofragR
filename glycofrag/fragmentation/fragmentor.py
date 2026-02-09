# glycofrag/fragmentation/fragmentor.py
"""Glycan fragmentation and fragment ion generation."""

from typing import List, Optional, Dict, Any, Tuple, cast
from itertools import combinations
import networkx as nx

from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.modifications import ReducingEndType
from glycofrag.core.constants import WATER_MASS, METHYL_MASS


class GlycanFragmentor:
    """
    Generates fragment ions from glycan structures.
    
    Responsibility:
    - Generate B, Y, BY, YY, BYY, YYY, BYYY series fragments
    - Generate C, Z, CZ, ZZ, CZZ, BZZ, ZZZ, CZZZ, BZZZ series fragments
    - Generate oxonium diagnostic ions
    - Calculate m/z values for specified charge states
    
    Does NOT handle:
    - Structure prediction
    - Structure classification
    
    Note:
    The full fragmentation logic remains in glycofrag/glycan/glycan.py.
    This class serves as an interface for future refactoring to extract
    fragmentation algorithms into independent modules.

    Public API:
        - generate_fragments(...)
    """
    
    def __init__(self, glycan_code: str, glycan_type: str = 'N', modification_type: int = 0,
                 custom_reducing_end_mass: Optional[float] = None,
                 custom_mono_masses: Optional[Dict[str, float]] = None):
        """
        Initialize fragmentor.
        
        Args:
            glycan_code: Glycan composition code (e.g., '4501')
            glycan_type: 'N' for N-glycan or 'O' for O-glycan
            modification_type: Reducing end modification type (0-7)
            custom_reducing_end_mass: Custom mass (Da) for reducing end (type 7 only)
            custom_mono_masses: Custom mass additions per monosaccharide (type 7 only)
        """
        self.glycan_code = glycan_code
        self.glycan_type = glycan_type
        self.modification_type = modification_type
        self.custom_reducing_end_mass = custom_reducing_end_mass
        self.custom_mono_masses = custom_mono_masses
        self._mass_calculator = GlycanMassCalculator(
            modification_type=modification_type,
            custom_reducing_end_mass=custom_reducing_end_mass,
            custom_mono_masses=custom_mono_masses
        )
        
        # Composition attributes (populated by generate_fragments)
        self.hexnac_total: int = 0
        self.hex_total: int = 0
        self.fuc_total: int = 0
        self.neuac_total: int = 0
        
        # Composition attributes (populated by generate_fragments)
        self.hexnac_total: int = 0
        self.hex_total: int = 0
        self.fuc_total: int = 0
        self.neuac_total: int = 0
    
    def _is_permethylated_modification(self, modification_type: int) -> bool:
        """
        Check if modification type is permethylated.
        
        Args:
            modification_type: Reducing end modification type
            
        Returns:
            True if permethylated, False otherwise
        """
        return modification_type in (
            ReducingEndType.PERMETHYLATED_REDUCED,
            ReducingEndType.PERMETHYLATED_FREE,
            ReducingEndType.TWO_AB_PERMETHYLATED
        )
    
    def _count_branches(self, structure: nx.DiGraph) -> int:
        """
        Count number of branches in structure.
        
        Args:
            structure: NetworkX DiGraph
            
        Returns:
            Number of branches
        """
        return self.count_branches(structure, self.glycan_type)
    
    def count_branches(self, structure: nx.DiGraph, glycan_type: str = 'N') -> int:
        """
        Count number of branches/arms in structure.
        
        Args:
            structure: NetworkX DiGraph of glycan structure
            glycan_type: 'N' or 'O'
        
        Returns:
            Number of branches
        """
        branch_count = 0
        
        if glycan_type == 'N':
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

    def generate_fragments(
        self,
        structure: nx.DiGraph,
        fragment_types: Optional[List[str]] = None,
        modification_type: Optional[int] = None,
        charges: Optional[List[int]] = None
    ) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, str]]]:
        """
        Generate fragment ions from a glycan structure with m/z values.
        
        This is the main entry point for fragment generation. It orchestrates
        calling the individual series generators and calculating m/z values.
        
        Args:
            structure: NetworkX DiGraph representing the glycan structure
            fragment_types: List of fragment series to generate ['BY', 'CZ'] or None for BY.
                A-type oxonium and custom diagnostic ions are always generated automatically.
            modification_type: Reducing end modification type. If None, uses instance modification_type
            charges: List of charge states to calculate m/z for (default [1, 2, 3])
        
        Returns:
            Tuple of (fragments dict, cleavage_info dict)
        """
        # Set defaults
        if charges is None:
            charges = [1, 2, 3]
        
        if modification_type is None:
            modification_type = self.modification_type
        elif modification_type != self.modification_type:
            self.modification_type = modification_type
            self._mass_calculator = GlycanMassCalculator(modification_type=modification_type)
        
        # Calculate composition from structure
        self._calculate_composition(structure)
        
        if fragment_types is None:
            fragment_types = ['BY']  # Default to BY series
        
        # Initialize fragment containers
        fragments = {}
        cleavage_info = {}
        
        # Generate BY series if requested
        if 'BY' in fragment_types:
            by_frags, by_info = self._generate_by_series(structure, modification_type)
            fragments.update(by_frags)
            cleavage_info.update(by_info)
        
        # Generate CZ series if requested
        if 'CZ' in fragment_types:
            cz_frags, cz_info = self._generate_cz_series(structure, modification_type)
            fragments.update(cz_frags)
            cleavage_info.update(cz_info)
        
        # Always generate A-type oxonium and custom diagnostic ions as add-ons
        a_frags, a_info = self._generate_oxonium_ions(structure, by_fragments=fragments)
        fragments.update(a_frags)
        cleavage_info.update(a_info)
        
        # Calculate mass for each fragment and add m/z values
        structure_label = None
        if isinstance(structure, nx.DiGraph):
            structure_label = structure.graph.get('structure_number')
        for frag_type, frag_list in fragments.items():
            for frag in frag_list:
                if structure_label is not None and '_structure_number' not in frag:
                    frag['_structure_number'] = structure_label
                # Handle oxonium ions with direct m/z values
                if '_direct_mz' in frag:
                    # _direct_mz is the z=1 m/z value, so subtract proton mass to get neutral mass
                    proton_mass = 1.007825
                    frag_mass = frag['_direct_mz'] - proton_mass
                    frag['mass'] = frag_mass
                    
                    # Calculate m/z for each charge state from the neutral mass
                    for charge in charges:
                        mz = self._mass_calculator.calculate_mz(frag_mass, charge)
                        frag[f'mz_charge_{charge}'] = mz
                else:
                    # Normal fragment mass calculation
                    frag_mass = self._mass_calculator.calculate_fragment_mass(frag, frag_type)
                    frag['mass'] = frag_mass
                    
                    # Calculate m/z for each charge state
                    for charge in charges:
                        mz = self._mass_calculator.calculate_mz(frag_mass, charge)
                        frag[f'mz_charge_{charge}'] = mz
        
        return fragments, cleavage_info

    def _generate_by_series(
        self,
        structure: nx.DiGraph,
        modification_type: int = 0
    ) -> Tuple[Dict[str, List[Dict[str, int]]], Dict[str, Dict[str, str]]]:
        """
        Generate BY-series fragments: B, Y, BY, YY, BYY, YYY, BYYY.
        
        Key Principle:
        - B-ions: Do NOT contain reducing end (node 1)
        - Y-ions: ALWAYS contain reducing end (node 1) - PLACEHOLDER
        - BY-ions: B-ions from double cleavages
        - YY-ions: Y-ions with additional losses (2 cleavages, requires 2+ branches)
        - BYY-ions: BY-ions with additional losses
        - YYY-ions: YY-ions with further losses (3 cleavages, requires 3+ branches)
        - BYYY-ions: YYY-ions from branched losses (requires 3+ branches)
        
        Returns:
            Tuple of (fragments dict, cleavage_info dict)
        """
        permethylated = self._is_permethylated_modification(modification_type)
        
        # Count branches to determine which fragment types are possible
        num_branches = self._count_branches(structure)

        fragments = {
            'b_ions': [],
            'y_ions': [],
            'by_ions': [],
            'yy_ions': [],
            'byy_ions': [],
            'yyy_ions': [],
            'byyy_ions': []
        }
        
        cleavage_info = {
            'b_ions': {},
            'y_ions': {},
            'by_ions': {},
            'yy_ions': {},
            'byy_ions': {},
            'yyy_ions': {},
            'byyy_ions': {}
        }
        
        # Track unique fragments
        unique_sets = {
            'b_ions': set(),
            'y_ions': set(),
            'by_ions': set(),
            'yy_ions': set(),
            'byy_ions': set(),
            'yyy_ions': set(),
            'byyy_ions': set()
        }
        
        # Reducing end is always node 1 (PLACEHOLDER for peptide)
        reducing_end = 1
        
        # Step 1: Generate B and Y ions from single cleavages
        for u, v in structure.edges():
            temp_graph = structure.copy()
            temp_graph.remove_edge(u, v)
            
            components = list(nx.weakly_connected_components(temp_graph))
            
            for component in components:
                comp_counts = self._count_residues(temp_graph.subgraph(component))
                
                if reducing_end in component:
                    # Y-ion (contains REDUCING END PLACEHOLDER)
                    # Mark this as a Y-ion so mass calculation can add reducing end mass
                    comp_counts['_is_y_ion'] = True
                    y_str = self._format_fragment_string(comp_counts, 'Y', is_y_ion=True)
                    if y_str not in unique_sets['y_ions'] and sum(comp_counts.values()) > 0:
                        unique_sets['y_ions'].add(y_str)
                        comp_counts['_cleavage_count'] = 1
                        fragments['y_ions'].append(comp_counts)
                        cleavage_info['y_ions'][y_str] = f"Cleavage_{u}_{v} (1_cleavage)"
                else:
                    # B-ion (non-reducing end, oxonium ion)
                    b_str = self._format_fragment_string(comp_counts, 'B')
                    if b_str not in unique_sets['b_ions'] and sum(comp_counts.values()) > 0:
                        unique_sets['b_ions'].add(b_str)
                        comp_counts['_cleavage_count'] = 1
                        fragments['b_ions'].append(comp_counts)
                        cleavage_info['b_ions'][b_str] = f"Cleavage_{u}_{v} (1_cleavage)"
                        
                        # Add neutral loss variants for special B-ions
                        neutral_variants = self._generate_neutral_loss_variants(comp_counts, permethylated)
                        for variant in neutral_variants:
                            variant_comp = {k: v for k, v in variant.items() if not k.startswith('_')}
                            variant_str = self._format_fragment_string(variant_comp, 'B') + f"(-{variant.get('_neutral_loss', 'H2O')})"
                            if variant_str not in unique_sets['b_ions']:
                                unique_sets['b_ions'].add(variant_str)
                                variant['_cleavage_count'] = 1
                                fragments['b_ions'].append(variant)
                                loss_type = variant.get('_neutral_loss', 'H2O')
                                cleavage_info['b_ions'][variant_str] = f"Cleavage_{u}_{v}_minus_{loss_type} (1_cleavage)"
        
        # Step 2: Generate BY ions (double cleavages, no reducing end)
        edges_list = list(structure.edges())
        for i, (u1, v1) in enumerate(edges_list):
            for u2, v2 in edges_list[i+1:]:
                temp_graph = structure.copy()
                temp_graph.remove_edge(u1, v1)
                temp_graph.remove_edge(u2, v2)
                
                components = list(nx.weakly_connected_components(temp_graph))
                
                for component in components:
                    if reducing_end not in component:
                        # BY fragments must be bounded by BOTH cleavages.
                        # If only one cleavage separates the component, it's a B-ion.
                        cut_edges = [(u1, v1), (u2, v2)]
                        boundary_cuts = 0
                        for a, b in cut_edges:
                            in_a = a in component
                            in_b = b in component
                            if in_a != in_b:
                                boundary_cuts += 1

                        if boundary_cuts < 2:
                            continue

                        comp_counts = self._count_residues(temp_graph.subgraph(component))
                        by_str = self._format_fragment_string(comp_counts, 'BY')
                        if by_str not in unique_sets['by_ions'] and sum(comp_counts.values()) > 0:
                            unique_sets['by_ions'].add(by_str)
                            comp_counts['_cleavage_count'] = 2
                            fragments['by_ions'].append(comp_counts)
                            cleavage_info['by_ions'][by_str] = f"Double_{u1}_{v1}_{u2}_{v2} (2_cleavages)"
        
        # Step 3: Generate YY ions (double Y cleavages) using topology-based edge pairs
        total_monos = self.hexnac_total + self.hex_total + self.fuc_total + self.neuac_total
        branch_dict = self._get_core_branches(structure)
        
        # Count antennae for 3+ antenna structures
        num_antennae = self._count_antennae(structure)
        
        for i, (u1, v1) in enumerate(edges_list):
            for u2, v2 in edges_list[i+1:]:
                # For N-glycans with 3+ antennae, use antenna-level checking
                if self.glycan_type == "N" and num_antennae >= 3:
                    antenna_dict = self._get_antenna_branches(structure)
                    ant1 = self._get_edge_branch(u1, v1, antenna_dict)
                    ant2 = self._get_edge_branch(u2, v2, antenna_dict)
                    
                    # Both edges must be in antennae AND from different antennae
                    if ant1 is None or ant2 is None or ant1 == ant2:
                        continue
                elif self.glycan_type == "N" and branch_dict:
                    # For 2-antenna structures, use branch-level checking
                    branch1 = self._get_edge_branch(u1, v1, branch_dict)
                    branch2 = self._get_edge_branch(u2, v2, branch_dict)
                    
                    # Both edges must be in branches AND from different branches
                    if branch1 is None or branch2 is None or branch1 == branch2:
                        continue
                
                temp_graph = structure.copy()
                try:
                    temp_graph.remove_edge(u1, v1)
                    temp_graph.remove_edge(u2, v2)
                except nx.NetworkXError:
                    continue

                components = list(nx.weakly_connected_components(temp_graph))
                reducing_component = None
                for component in components:
                    if reducing_end in component:
                        reducing_component = component
                        break

                if reducing_component is None:
                    continue

                # Validate that all nodes in reducing_component are reachable from root
                if not self._is_valid_yy_composition(temp_graph, reducing_component):
                    continue

                yy_counts = self._count_residues(temp_graph.subgraph(reducing_component))
                yy_str_check = self._format_fragment_string(yy_counts, 'YY', is_y_ion=True)
                
                if sum(yy_counts.values()) == 0:
                    continue

                if sum(yy_counts.values()) >= total_monos:
                    continue

                if self.glycan_type == "N":
                    # For 3+ antennae: allow smaller core fragments
                    min_hexnac = 1 if num_antennae >= 3 else 2
                    min_hex = 1 if num_antennae >= 3 else 2
                    if yy_counts.get('HexNAc', 0) < min_hexnac or yy_counts.get('Hex', 0) < min_hex:
                        continue

                yy_str = self._format_fragment_string(yy_counts, 'YY', is_y_ion=True)
                if yy_str not in unique_sets['yy_ions']:
                    unique_sets['yy_ions'].add(yy_str)
                    yy_counts['_cleavage_count'] = 2
                    fragments['yy_ions'].append(yy_counts)
                    cleavage_info['yy_ions'][yy_str] = f"Double_{u1}_{v1}_{u2}_{v2} (2_cleavages)"

        # Step 3b: Generate YY ions for 3+ antenna structures using B+Y cleavages
        # 1 B-cleavage (removes entire antenna) + 1 Y-cleavage (removes terminal part)
        if self.glycan_type == "N" and num_antennae >= 3:
            antenna_dict = self._get_antenna_branches(structure)
            
            for i, (u1, v1) in enumerate(edges_list):
                for u2, v2 in edges_list[i+1:]:
                    # Get antenna assignments for both edges
                    ant1 = self._get_edge_branch(u1, v1, antenna_dict)
                    ant2 = self._get_edge_branch(u2, v2, antenna_dict)
                    
                    # At least one edge should be in antenna, and if both are in antennae, they must be different
                    if ant1 is None and ant2 is None:
                        continue
                    if ant1 is not None and ant2 is not None and ant1 == ant2:
                        continue
                    
                    # Try both: edge1 as B-cleavage, edge2 as Y-cleavage
                    for b_edge, y_edge in [((u1, v1), (u2, v2)), ((u2, v2), (u1, v1))]:
                        temp_graph = structure.copy()
                        try:
                            temp_graph.remove_edge(b_edge[0], b_edge[1])
                            temp_graph.remove_edge(y_edge[0], y_edge[1])
                        except nx.NetworkXError:
                            continue
                        
                        components = list(nx.weakly_connected_components(temp_graph))
                        reducing_component = None
                        for component in components:
                            if reducing_end in component:
                                reducing_component = component
                                break
                        
                        if reducing_component is None:
                            continue
                        
                        # Validate reachability
                        if not self._is_valid_yy_composition(temp_graph, reducing_component):
                            continue
                        
                        yy_counts = self._count_residues(temp_graph.subgraph(reducing_component))
                        
                        if sum(yy_counts.values()) == 0:
                            continue
                        
                        if sum(yy_counts.values()) >= total_monos:
                            continue
                        
                        # For B+Y combinations: allow HexNAc >= 1, Hex >= 1
                        if yy_counts.get('HexNAc', 0) < 1 or yy_counts.get('Hex', 0) < 1:
                            continue
                        
                        yy_str = self._format_fragment_string(yy_counts, 'YY', is_y_ion=True)
                        if yy_str not in unique_sets['yy_ions']:
                            unique_sets['yy_ions'].add(yy_str)
                            yy_counts['_cleavage_count'] = 2
                            fragments['yy_ions'].append(yy_counts)
                            cleavage_info['yy_ions'][yy_str] = f"B_Y_{b_edge[0]}_{b_edge[1]}_{y_edge[0]}_{y_edge[1]} (2_cleavages)"

        # Ensure core YY fragment for N-glycans (HexNAc2Hex2)
        if self.glycan_type == "N" and self.hexnac_total >= 2 and self.hex_total >= 2:
            yy_core = {'HexNAc': 2, 'Hex': 2}
            yy_core['_cleavage_count'] = 2
            yy_core_str = self._format_fragment_string(yy_core, 'YY', is_y_ion=True)
            if yy_core_str not in unique_sets['yy_ions']:
                unique_sets['yy_ions'].add(yy_core_str)
                fragments['yy_ions'].append(yy_core)
                cleavage_info['yy_ions'][yy_core_str] = "Core_YY_HexNAc2Hex2 (2_cleavages)"

        # Ensure core YY fragment for N-glycans (HexNAc2Hex1)
        if self.glycan_type == "N" and self.hexnac_total >= 2 and self.hex_total >= 1:
            yy_core_1 = {'HexNAc': 2, 'Hex': 1}
            yy_core_1['_cleavage_count'] = 2
            yy_core_1_str = self._format_fragment_string(yy_core_1, 'YY', is_y_ion=True)
            if yy_core_1_str not in unique_sets['yy_ions']:
                unique_sets['yy_ions'].add(yy_core_1_str)
                fragments['yy_ions'].append(yy_core_1)
                cleavage_info['yy_ions'][yy_core_1_str] = "Core_YY_HexNAc2Hex1 (2_cleavages)"

        # Step 3c: For multi-antenna structures (HexNAc > 4), generate topological BYY
        # by taking YY fragments (Y-Y on arms) and adding B-cleavage at 3-4 (removing chitobiose)
        # This keeps the core Hex but removes the 2 HexNAc from chitobiose
        if num_antennae >= 3 and self.glycan_type == "N":
            for yy_ion in fragments['yy_ions'].copy():
                # Check if this YY has the chitobiose (HexNAc2)
                if yy_ion.get('HexNAc', 0) >= 2:
                    # B-cleavage at 3-4: removes chitobiose (2 HexNAc) but keeps core Hex
                    byy_topo = yy_ion.copy()
                    # Remove any stale mass/mz fields to force recalculation
                    byy_topo.pop('mass', None)
                    for key in list(byy_topo.keys()):
                        if key.startswith('mz_'):
                            byy_topo.pop(key)
                    byy_topo['HexNAc'] = byy_topo.get('HexNAc', 0) - 2
                    if byy_topo.get('HexNAc', 0) == 0:
                        byy_topo.pop('HexNAc', None)
                    
                    if sum(v for k, v in byy_topo.items() if not k.startswith('_')) > 0:
                        # Validate: for < 4 antennae, max 1 NeuAc
                        # For 4+ antennae: multiple NeuAc possible
                        if num_antennae < 4 and byy_topo.get('NeuAc', 0) > 1:
                            continue
                        byy_topo_str = self._format_fragment_string(byy_topo, 'BYY')
                        if byy_topo_str not in unique_sets['byy_ions']:
                            unique_sets['byy_ions'].add(byy_topo_str)
                            byy_topo['_cleavage_count'] = 3
                            fragments['byy_ions'].append(byy_topo)
                            cleavage_info['byy_ions'][byy_topo_str] = "BYY_topological_YY_plus_B_at_3_4 (3_cleavages)"

        # Step 4: Generate BYY ions from YY ions by removing 1 or 2 core HexNAc units
        for yy_ion in fragments['yy_ions'].copy():
            if yy_ion.get('HexNAc', 0) >= 1:
                byy_comp = yy_ion.copy()
                # Remove any stale mass/mz fields to force recalculation
                byy_comp.pop('mass', None)
                for key in list(byy_comp.keys()):
                    if key.startswith('mz_'):
                        byy_comp.pop(key)
                byy_comp['HexNAc'] = max(0, byy_comp.get('HexNAc', 0) - 1)
                if byy_comp.get('HexNAc', 0) == 0:
                    byy_comp.pop('HexNAc', None)

                if sum(v for k, v in byy_comp.items() if not k.startswith('_')) > 0:
                    # Filter out impossible NeuAc combinations
                    # For < 4 antennae: max 1 NeuAc can remain after 2 Y-cleavages
                    # For 4+ antennae: multiple NeuAc can remain
                    if num_antennae < 4 and byy_comp.get('NeuAc', 0) > 1:
                        continue
                    
                    byy_str = self._format_fragment_string(byy_comp, 'BYY')
                    if byy_str not in unique_sets['byy_ions']:
                        unique_sets['byy_ions'].add(byy_str)
                        byy_comp['_cleavage_count'] = 2
                        fragments['byy_ions'].append(byy_comp)
                        cleavage_info['byy_ions'][byy_str] = "BYY_from_YY_minus_HexNAc (derived_from_2_cleavages)"

            if self.glycan_type != "O" and yy_ion.get('HexNAc', 0) >= 2:
                byy_comp2 = yy_ion.copy()
                # Remove any stale mass/mz fields to force recalculation
                byy_comp2.pop('mass', None)
                for key in list(byy_comp2.keys()):
                    if key.startswith('mz_'):
                        byy_comp2.pop(key)
                byy_comp2['HexNAc'] = max(0, byy_comp2.get('HexNAc', 0) - 2)
                if byy_comp2.get('HexNAc', 0) == 0:
                    byy_comp2.pop('HexNAc', None)

                if sum(v for k, v in byy_comp2.items() if not k.startswith('_')) > 0:
                    # Filter out impossible NeuAc combinations
                    # For < 4 antennae: max 1 NeuAc can remain
                    # For 4+ antennae: multiple NeuAc possible
                    if num_antennae < 4 and byy_comp2.get('NeuAc', 0) > 1:
                        continue
                    
                    byy_str2 = self._format_fragment_string(byy_comp2, 'BYY')
                    if byy_str2 not in unique_sets['byy_ions']:
                        unique_sets['byy_ions'].add(byy_str2)
                        byy_comp2['_cleavage_count'] = 2
                        fragments['byy_ions'].append(byy_comp2)
                        cleavage_info['byy_ions'][byy_str2] = "BYY_from_YY_minus_2HexNAc (derived_from_2_cleavages)"

            # Additional BYY from YY by removing one Fuc (B cleavage on fucose)
            if yy_ion.get('Fuc', 0) >= 1:
                byy_fuc = yy_ion.copy()
                byy_fuc['Fuc'] = max(0, byy_fuc.get('Fuc', 0) - 1)
                if byy_fuc.get('Fuc', 0) == 0:
                    byy_fuc.pop('Fuc', None)

                if sum(v for k, v in byy_fuc.items() if not k.startswith('_')) > 0:
                    # Filter out impossible NeuAc combinations
                    # For < 4 antennae: max 1 NeuAc can remain
                    # For 4+ antennae: multiple NeuAc possible
                    if num_antennae < 4 and byy_fuc.get('NeuAc', 0) > 1:
                        continue
                    
                    byy_fuc_str = self._format_fragment_string(byy_fuc, 'BYY')
                    if byy_fuc_str not in unique_sets['byy_ions']:
                        unique_sets['byy_ions'].add(byy_fuc_str)
                        byy_fuc['_cleavage_count'] = 2
                        fragments['byy_ions'].append(byy_fuc)
                        cleavage_info['byy_ions'][byy_fuc_str] = "BYY_from_YY_minus_Fuc (derived_from_2_cleavages)"

        # Ensure BYY core fragment for multi-antenna structures (HexNAc2Hex)
        if num_antennae >= 3 and self.glycan_type == "N" and self.hexnac_total >= 2 and self.hex_total >= 1:
            byy_core = {'HexNAc': 2, 'Hex': 1}
            byy_core['_cleavage_count'] = 2
            byy_core_str = self._format_fragment_string(byy_core, 'BYY')
            if byy_core_str not in unique_sets['byy_ions']:
                unique_sets['byy_ions'].add(byy_core_str)
                fragments['byy_ions'].append(byy_core)
                cleavage_info['byy_ions'][byy_core_str] = "BYY_core_HexNAc2Hex1_ensure (derived_from_2_cleavages)"


        # Step 5: Generate YYY ions (triple Y cleavages) using topology-based edge triplets
        # Enable for: complex structures (HexNAc>4) OR with 3+ antennae OR core-fucose extra cleavage
        core_fuc_edges = set()
        if self.glycan_type == "N":
            for u, v in structure.edges():
                if structure.nodes[v].get('type') == 'Fuc':
                    parent_pos = structure.nodes[u].get('position')
                    if parent_pos == 'core_reducing':
                        core_fuc_edges.add((u, v))
        
        # Count actual antennae (HexNAc chains from core branches)
        num_antennae = self._count_antennae(structure)

        allow_yyy = (
            (self.hexnac_total > 4)
            or (self.hex_total > 4 and self.hexnac_total >= 2)
            or (num_antennae >= 3)
            or bool(core_fuc_edges)
        )

        if allow_yyy:
            # For 3+ antennae, build antenna-based branch mapping
            if num_antennae >= 3:
                antenna_dict = self._get_antenna_branches(structure)
            else:
                antenna_dict = branch_dict
            
            for (u1, v1), (u2, v2), (u3, v3) in combinations(edges_list, 3):
                branch1 = self._get_edge_branch(u1, v1, branch_dict)
                branch2 = self._get_edge_branch(u2, v2, branch_dict)
                branch3 = self._get_edge_branch(u3, v3, branch_dict)
                
                # For 3+ antennae, also check antenna assignment
                if num_antennae >= 3:
                    ant1 = self._get_edge_branch(u1, v1, antenna_dict)
                    ant2 = self._get_edge_branch(u2, v2, antenna_dict)
                    ant3 = self._get_edge_branch(u3, v3, antenna_dict)
                    
                    # Require all 3 edges from different antennae
                    if ant1 is None or ant2 is None or ant3 is None:
                        continue
                    if ant1 == ant2 or ant2 == ant3 or ant1 == ant3:
                        continue
                elif len(branch_dict) >= 3:
                    # Require all 3 edges from different branches
                    if branch1 is None or branch2 is None or branch3 is None:
                        continue
                    if branch1 == branch2 or branch2 == branch3 or branch1 == branch3:
                        continue
                else:
                    # Allow YYY if one edge is core fucose cleavage and the other two are from different branches
                    edge_triplet = [(u1, v1), (u2, v2), (u3, v3)]
                    if not any(edge in core_fuc_edges for edge in edge_triplet):
                        continue

                    branch_ids = [b for b in [branch1, branch2, branch3] if b is not None]
                    if len(set(branch_ids)) < 2:
                        continue
                
                temp_graph = structure.copy()
                try:
                    temp_graph.remove_edge(u1, v1)
                    temp_graph.remove_edge(u2, v2)
                    temp_graph.remove_edge(u3, v3)
                except nx.NetworkXError:
                    continue

                components = list(nx.weakly_connected_components(temp_graph))
                reducing_component = None
                for component in components:
                    if reducing_end in component:
                        reducing_component = component
                        break

                if reducing_component is None:
                    continue

                yyy_counts = self._count_residues(temp_graph.subgraph(reducing_component))
                if sum(yyy_counts.values()) == 0:
                    continue

                if sum(yyy_counts.values()) >= total_monos:
                    continue

                if self.glycan_type == "N":
                    # For 3+ antennae: allow Hex >= 2 (core can have partial structure)
                    # For 2 antennae: require Hex >= 3 (standard minimum)
                    min_hex = 2 if num_antennae >= 3 else 3
                    if yyy_counts.get('HexNAc', 0) < 2 or yyy_counts.get('Hex', 0) < min_hex:
                        continue

                yyy_str = self._format_fragment_string(yyy_counts, 'YYY', is_y_ion=True)
                if yyy_str not in unique_sets['yyy_ions']:
                    unique_sets['yyy_ions'].add(yyy_str)
                    yyy_counts['_cleavage_count'] = 3
                    fragments['yyy_ions'].append(yyy_counts)
                    cleavage_info['yyy_ions'][yyy_str] = f"Triple_{u1}_{v1}_{u2}_{v2}_{u3}_{v3} (3_cleavages)"

            # Ensure core YYY fragment for N-glycans with 3+ antennae (HexNAc2Hex2)
            # Only generate for 3+ antennae: 2-antenna structures reach HexNAc2Hex2 with just 2 cleavages (YY, not YYY)
            if self.glycan_type == "N" and num_antennae >= 3 and self.hexnac_total >= 2 and self.hex_total >= 2:
                yyy_core = {'HexNAc': 2, 'Hex': 2, '_cleavage_count': 3}
                yyy_core_str = self._format_fragment_string(yyy_core, 'YYY', is_y_ion=True)
                if yyy_core_str not in unique_sets['yyy_ions']:
                    unique_sets['yyy_ions'].add(yyy_core_str)
                    fragments['yyy_ions'].append(yyy_core)
                    cleavage_info['yyy_ions'][yyy_core_str] = "Core_YYY_HexNAc2Hex2 (3_cleavages)"

            # Ensure core YYY fragments for 3+ antennae structures
            if num_antennae >= 3:
                # HexNAc3Hex2 (common for tri-antennary)
                if self.hexnac_total >= 3 and self.hex_total >= 2:
                    yyy_3ant_1 = {'HexNAc': 3, 'Hex': 2, '_cleavage_count': 3}
                    yyy_3ant_1_str = self._format_fragment_string(yyy_3ant_1, 'YYY', is_y_ion=True)
                    if yyy_3ant_1_str not in unique_sets['yyy_ions']:
                        unique_sets['yyy_ions'].add(yyy_3ant_1_str)
                        fragments['yyy_ions'].append(yyy_3ant_1)
                        cleavage_info['yyy_ions'][yyy_3ant_1_str] = "Core_YYY_3antenna_HexNAc3Hex2 (3_cleavages)"
                
                # HexNAc4Hex2 (for larger structures)
                if self.hexnac_total >= 4 and self.hex_total >= 2:
                    yyy_3ant_2 = {'HexNAc': 4, 'Hex': 2, '_cleavage_count': 3}
                    yyy_3ant_2_str = self._format_fragment_string(yyy_3ant_2, 'YYY', is_y_ion=True)
                    if yyy_3ant_2_str not in unique_sets['yyy_ions']:
                        unique_sets['yyy_ions'].add(yyy_3ant_2_str)
                        fragments['yyy_ions'].append(yyy_3ant_2)
                        cleavage_info['yyy_ions'][yyy_3ant_2_str] = "Core_YYY_3antenna_HexNAc4Hex2 (3_cleavages)"
                
                # HexNAc2Hex3 (another common pattern)
                if self.hexnac_total >= 2 and self.hex_total >= 3:
                    yyy_3ant_3 = {'HexNAc': 2, 'Hex': 3, '_cleavage_count': 3}
                    yyy_3ant_3_str = self._format_fragment_string(yyy_3ant_3, 'YYY', is_y_ion=True)
                    if yyy_3ant_3_str not in unique_sets['yyy_ions']:
                        unique_sets['yyy_ions'].add(yyy_3ant_3_str)
                        fragments['yyy_ions'].append(yyy_3ant_3)
                        cleavage_info['yyy_ions'][yyy_3ant_3_str] = "Core_YYY_3antenna_HexNAc2Hex3 (3_cleavages)"

            # Ensure core YYY fragment for branched fucose cases (HexNAc3Hex2)
            if (self.glycan_type == "N" and self.fuc_total > 0 and
                self.hexnac_total >= 3 and self.hex_total >= 2):
                yyy_core_fuc = {'HexNAc': 3, 'Hex': 2, '_cleavage_count': 3}
                yyy_core_fuc_str = self._format_fragment_string(yyy_core_fuc, 'YYY', is_y_ion=True)
                if yyy_core_fuc_str not in unique_sets['yyy_ions']:
                    unique_sets['yyy_ions'].add(yyy_core_fuc_str)
                    fragments['yyy_ions'].append(yyy_core_fuc)
                    cleavage_info['yyy_ions'][yyy_core_fuc_str] = "Core_YYY_HexNAc3Hex2 (3_cleavages)"

            # Ensure core YYY fragment for core-fucose cases (HexNAc2Hex1)
            if (self.glycan_type == "N" and self.fuc_total > 0 and
                self.hexnac_total >= 2 and self.hex_total >= 1):
                yyy_core_fuc_2 = {'HexNAc': 2, 'Hex': 1, '_cleavage_count': 3}
                yyy_core_fuc_2_str = self._format_fragment_string(yyy_core_fuc_2, 'YYY', is_y_ion=True)
                if yyy_core_fuc_2_str not in unique_sets['yyy_ions']:
                    unique_sets['yyy_ions'].add(yyy_core_fuc_2_str)
                    fragments['yyy_ions'].append(yyy_core_fuc_2)
                    cleavage_info['yyy_ions'][yyy_core_fuc_2_str] = "Core_YYY_HexNAc2Hex1 (3_cleavages)"
            
            # For High Mannose structures (HexNAc=2), generate additional YYY fragments
            # These represent triple cleavages on the branched mannose tree
            # Only add fragments that REQUIRE 3 cleavages (not achievable with 2 or fewer)
            if self.glycan_type == "N" and self.hexnac_total == 2 and self.hex_total >= 4:
                # HexNAc2Hex3 is the main YYY fragment requiring 3 cleavages in high mannose
                # (need to cleave branches from both core mannose arms plus one terminal)
                if self.hex_total >= 4:
                    yyy_hm = {'HexNAc': 2, 'Hex': 3, '_cleavage_count': 3}
                    yyy_hm_str = self._format_fragment_string(yyy_hm, 'YYY', is_y_ion=True)
                    if yyy_hm_str not in unique_sets['yyy_ions']:
                        unique_sets['yyy_ions'].add(yyy_hm_str)
                        fragments['yyy_ions'].append(yyy_hm)
                        cleavage_info['yyy_ions'][yyy_hm_str] = "HighMannose_YYY_HexNAc2Hex3 (3_cleavages)"

        # Step 6: Generate BYYY ions from YYY ions by removing 1 or 2 core HexNAc units
        for yyy_ion in fragments['yyy_ions'].copy():
            # Remove 1 HexNAc
            if yyy_ion.get('HexNAc', 0) >= 1:
                byyy_comp = yyy_ion.copy()
                byyy_comp['HexNAc'] = max(0, byyy_comp.get('HexNAc', 0) - 1)
                if byyy_comp.get('HexNAc', 0) == 0:
                    byyy_comp.pop('HexNAc', None)

                if sum(v for k, v in byyy_comp.items() if not k.startswith('_')) > 0:
                    byyy_str = self._format_fragment_string(byyy_comp, 'BYYY')
                    if byyy_str not in unique_sets['byyy_ions']:
                        unique_sets['byyy_ions'].add(byyy_str)
                        byyy_comp['_cleavage_count'] = 3
                        fragments['byyy_ions'].append(byyy_comp)
                        cleavage_info['byyy_ions'][byyy_str] = "BYYY_from_YYY_minus_HexNAc (derived_from_3_cleavages)"

            # Remove 2 HexNAc
            if yyy_ion.get('HexNAc', 0) >= 2:
                byyy_comp2 = yyy_ion.copy()
                byyy_comp2['HexNAc'] = max(0, byyy_comp2.get('HexNAc', 0) - 2)
                if byyy_comp2.get('HexNAc', 0) == 0:
                    byyy_comp2.pop('HexNAc', None)

                if sum(v for k, v in byyy_comp2.items() if not k.startswith('_')) > 0:
                    byyy_str2 = self._format_fragment_string(byyy_comp2, 'BYYY')
                    if byyy_str2 not in unique_sets['byyy_ions']:
                        unique_sets['byyy_ions'].add(byyy_str2)
                        byyy_comp2['_cleavage_count'] = 3
                        fragments['byyy_ions'].append(byyy_comp2)
                        cleavage_info['byyy_ions'][byyy_str2] = "BYYY_from_YYY_minus_2HexNAc (derived_from_3_cleavages)"
        
        return fragments, cleavage_info
    
    def _generate_cz_series(
        self,
        structure: nx.DiGraph,
        modification_type: int = 0
    ) -> Tuple[Dict[str, List[Dict[str, int]]], Dict[str, Dict[str, str]]]:
        """
        Generate CZ-series fragments: C, Z, CZ, ZZ, CZZ, BZZ, ZZZ, CZZZ, BZZZ.
        
        CZ fragments are parallel to BY fragments with different mass shifts:
        - C-ions: Like B-ions but with +NH3 (add 17.026549 Da)
        - Z-ions: Like Y-ions but with -NH (add 2.015894 Da relative to Y)
        
        Implementation Strategy:
        1. Generate BY-series fragments first
        2. Create CZ equivalents with appropriate metadata flags
        3. Mass calculation happens in GlycanMassCalculator based on fragment type
        
        Returns:
            Tuple of (fragments dict, cleavage_info dict)
        """
        # First generate BY series as template
        by_fragments, by_cleavage_info = self._generate_by_series(structure, modification_type)
        
        # Initialize CZ fragments
        cz_fragments = {
            'c_ions': [],
            'z_ions': [],
            'cz_ions': [],
            'zz_ions': [],
            'czz_ions': [],
            'bzz_ions': [],
            'zzz_ions': [],
            'czzz_ions': [],
            'bzzz_ions': []
        }
        
        cz_cleavage_info = {
            'c_ions': {},
            'z_ions': {},
            'cz_ions': {},
            'zz_ions': {},
            'czz_ions': {},
            'bzz_ions': {},
            'zzz_ions': {},
            'czzz_ions': {},
            'bzzz_ions': {}
        }
        
        # Mapping from BY to CZ types
        mapping = {
            'b_ions': 'c_ions',
            'y_ions': 'z_ions',
            'by_ions': 'cz_ions',
            'yy_ions': 'zz_ions',
            'byy_ions': 'czz_ions',
            'yyy_ions': 'zzz_ions',
            'byyy_ions': 'czzz_ions'
        }
        
        # Additional: BZZ from BYY (special case)
        # BZZ ions are like BYY but in CZ nomenclature
        if 'byy_ions' in by_fragments:
            for byy_ion in by_fragments['byy_ions']:
                bzz_comp = byy_ion.copy()
                # Remove any stale mass/mz fields to force recalculation (CZ has different mass than BY)
                bzz_comp.pop('mass', None)
                for key in list(bzz_comp.keys()):
                    if key.startswith('mz_'):
                        bzz_comp.pop(key)
                bzz_str = self._format_fragment_string(bzz_comp, 'BZZ')
                cz_fragments['bzz_ions'].append(bzz_comp)
                cz_cleavage_info['bzz_ions'][bzz_str] = "BYY_to_BZZ (derived_from_2_or_3_cleavages)"
        
        # Map all BY fragments to CZ equivalents
        for by_type, cz_type in mapping.items():
            if by_type in by_fragments:
                for fragment in by_fragments[by_type]:
                    # Copy composition
                    cz_comp = fragment.copy()
                    # Remove any stale mass/mz fields to force recalculation (CZ has different mass than BY)
                    cz_comp.pop('mass', None)
                    for key in list(cz_comp.keys()):
                        if key.startswith('mz_'):
                            cz_comp.pop(key)
                    # Create CZ string
                    cz_str = self._format_fragment_string(cz_comp, cz_type.replace('_ions', '').upper())
                    cz_fragments[cz_type].append(cz_comp)
                    # Copy cleavage info
                    by_key = self._format_fragment_string(fragment, by_type.replace('_ions', '').upper())
                    if by_key in by_cleavage_info.get(by_type, {}):
                        cz_cleavage_info[cz_type][cz_str] = by_cleavage_info[by_type][by_key]
        
        # Add BZZZ from BYYY
        if 'byyy_ions' in by_fragments:
            for byyy_ion in by_fragments['byyy_ions']:
                bzzz_comp = byyy_ion.copy()
                # Remove any stale mass/mz fields to force recalculation (CZ has different mass than BY)
                bzzz_comp.pop('mass', None)
                for key in list(bzzz_comp.keys()):
                    if key.startswith('mz_'):
                        bzzz_comp.pop(key)
                bzzz_str = self._format_fragment_string(bzzz_comp, 'BZZZ')
                cz_fragments['bzzz_ions'].append(bzzz_comp)
                cz_cleavage_info['bzzz_ions'][bzzz_str] = "BYYY_to_BZZZ (derived_from_3_cleavages)"
        
        return cz_fragments, cz_cleavage_info
    
    def _get_y_ion_suffix(self) -> str:
        """
        Get the appropriate suffix for Y-type ions based on modification type.
        
        Returns:
            String suffix like '-Redend', '-FreeEnd', '-2AB', or '-PEP'
        """
        mod_type = self.modification_type
        
        if mod_type in (ReducingEndType.PERMETHYLATED_REDUCED, ReducingEndType.REDUCED):
            return "-Redend"
        elif mod_type in (ReducingEndType.PERMETHYLATED_FREE, ReducingEndType.FREE):
            return "-FreeEnd"
        elif mod_type in (ReducingEndType.TWO_AB, ReducingEndType.TWO_AB_PERMETHYLATED):
            return "-2AB"
        elif mod_type == ReducingEndType.GLYCOPEPTIDE:
            return "-PEP"
        else:
            # Default to Redend for unknown types
            return "-Redend"

    def _generate_neutral_loss_variants(
        self,
        composition: Dict[str, int],
        permethylated: bool = False
    ) -> List[Dict[str, Any]]:
        """
        Generate neutral loss variants for a fragment.
        
        Special fragments generated:
        - NeuAc-H2O: Water loss from sialic acid
        - Fuc-H2O: Water loss from fucose  
        - Fuc-CH3: Methyl loss from fucose (only for permethylated)
        - HexNAc-C2H6O3: 126 Da oxonium ion
        
        Args:
            composition: Fragment composition
            permethylated: Whether glycan is permethylated
        
        Returns:
            List of composition dicts with neutral loss metadata
        """
        variants = []
        
        # NeuAc-H2O loss (water from sialic acid)
        if composition.get('NeuAc', 0) == 1 and len(composition) == 1:
            variant: Dict[str, Any] = cast(Dict[str, Any], composition.copy())
            variant['_neutral_loss'] = 'H2O'
            variant['_mass_adjustment'] = -WATER_MASS
            variants.append(variant)
        
        # Fuc-H2O loss (water from fucose)
        if composition.get('Fuc', 0) == 1 and len(composition) == 1:
            variant: Dict[str, Any] = cast(Dict[str, Any], composition.copy())
            variant['_neutral_loss'] = 'H2O'
            variant['_mass_adjustment'] = -WATER_MASS
            variants.append(variant)
        
        # Fuc-CH3 loss (methyl from fucose, only permethylated)
        if permethylated and composition.get('Fuc', 0) == 1 and len(composition) == 1:
            variant: Dict[str, Any] = cast(Dict[str, Any], composition.copy())
            variant['_neutral_loss'] = 'CH3'
            variant['_mass_adjustment'] = -14.0157  # Methyl group
            variants.append(variant)
        
        # HexNAc-C2H6O3 (126 Da oxonium, important diagnostic)
        if composition.get('HexNAc', 0) == 1 and len(composition) == 1:
            variant: Dict[str, Any] = cast(Dict[str, Any], composition.copy())
            variant['_neutral_loss'] = 'C2H6O3'
            variant['_mass_adjustment'] = -78.0317  # C2H6O3
            variants.append(variant)
        
        return variants
    
    def _generate_oxonium_ions(
        self,
        structure: nx.DiGraph,
        by_fragments: Optional[Dict[str, List[Dict[str, Any]]]] = None
    ) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, str]]]:
        """
        Generate A-type oxonium ions (diagnostic glycan ions).
        
        Standard diagnostic ions (always generated):
        - HexNAc(-CH6O3) at 138.055 m/z (O-glycan marker)
        - HexNAc(-C2H4O2) at 144.066 m/z (N-glycan marker)
        - HexNAc-C₂H₆O₃ at 126.0166 m/z
        
        Custom diagnostic ions (generated from BY fragments):
        - Y-HexNAc1(-CH3OH)-B (from HexNAc1 in BY)
        - Y-HexNAc1(-H2O)-B (from HexNAc1 in BY)
        - Y-HexNAc1(-2H2O)-B (from HexNAc1 in BY)
        - NeuAc1(-H2O)-B (from NeuAc1 in B)
        - Fuc1(-H2O)-B (from Fuc1 in B)
        - Fuc1(-CH3)-B (from Fuc1 in B)
        
        Args:
            structure: The glycan structure graph
            by_fragments: Optional dict with 'b_ions' and 'by_ions' for custom generation
        
        Returns:
            Tuple of (fragments dict, cleavage_info dict)
        """
        fragments = {'a_ions': [], 'custom_ions': []}
        cleavage_info = {'a_ions': {}, 'custom_ions': {}}
        
        # Get composition to check what sugars are present
        total_comp = self._count_residues(structure)
        
        # Define masses for custom fragments
        WATER_MASS = 18.010565
        METHANOL_MASS = 32.026215
        METHYL_MASS = 15.023479
        ACETYL_PLUS_C2H4O2 = 78.0687
        TWO_WATER_MASS = 36.021129
        
        # Adjust for permethylation if needed
        if self.modification_type in [2, 3, 5]:
            PERMETHYLATION_MASS = 14.015650
            WATER_MASS = 18.0153 - PERMETHYLATION_MASS
            METHANOL_MASS = 32.042 - PERMETHYLATION_MASS
            ACETYL_PLUS_C2H4O2 = 78.07 - (3 * PERMETHYLATION_MASS)
            TWO_WATER_MASS = 36.02 - (2 * PERMETHYLATION_MASS)
        
        # ALWAYS generate these 3 static oxonium ions (for all glycan types)
        if total_comp.get('HexNAc', 0) > 0:
            # 1. HexNAc(-CH6O3) - O-glycan marker at 138.055 m/z (ALWAYS)
            oxonium_1 = {
                '_custom_label': 'HexNAc(-CH6O3)',
                '_direct_mz': 138.055,
                '_is_oxonium': True
            }
            fragments['a_ions'].append(oxonium_1)
            cleavage_info['a_ions']['HexNAc(-CH6O3)'] = "Diagnostic ion (all types)"
            
            # 2. HexNAc(-C2H4O2) - N-glycan marker at 144.066 m/z (ALWAYS)
            oxonium_2 = {
                '_custom_label': 'HexNAc(-C2H4O2)',
                '_direct_mz': 144.066,
                '_is_oxonium': True
            }
            fragments['a_ions'].append(oxonium_2)
            cleavage_info['a_ions']['HexNAc(-C2H4O2)'] = "Diagnostic ion (all types)"
            
            # 3. HexNAc-C2H6O3 at 126.0166 m/z (ALWAYS)
            oxonium_3 = {
                '_custom_label': 'HexNAc-C2H6O3',
                '_direct_mz': 126.0166,
                '_is_oxonium': True
            }
            fragments['a_ions'].append(oxonium_3)
            cleavage_info['a_ions']['HexNAc-C2H6O3'] = "Diagnostic ion (all types)"
        
        # Generate custom diagnostic ions from BY fragments if provided
        if by_fragments:
            mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
            
            def _mono_total(fragment: Dict[str, Any]) -> int:
                return sum(fragment.get(m, 0) for m in mono_keys)

            # Y-HexNAc1-B custom fragments (variants with neutral losses)
            by_ions = by_fragments.get('by_ions', [])
            for fragment in by_ions:
                if (
                    isinstance(fragment, dict)
                    and fragment.get('HexNAc', 0) == 1
                    and _mono_total(fragment) == 1
                ):
                    base_hexnac = {'HexNAc': 1}

                    # Y-HexNAc1(-CH3OH)-B
                    custom_hexnac_methanol = {
                        **base_hexnac,
                        '_custom_label': 'Y-HexNAc1(-CH3OH)-B',
                        '_mass_adjustment': -METHANOL_MASS,
                        '_is_oxonium': True,
                        '_custom_type': 'dehydration'
                    }
                    fragments['custom_ions'].append(custom_hexnac_methanol)
                    cleavage_info['custom_ions']['Y-HexNAc1(-CH3OH)-B'] = "HexNAc custom loss"
                    
                    # Y-HexNAc1(-H2O)-B
                    custom_hexnac_water = {
                        **base_hexnac,
                        '_custom_label': 'Y-HexNAc1(-H2O)-B',
                        '_mass_adjustment': -WATER_MASS,
                        '_is_oxonium': True,
                        '_custom_type': 'dehydration'
                    }
                    fragments['custom_ions'].append(custom_hexnac_water)
                    cleavage_info['custom_ions']['Y-HexNAc1(-H2O)-B'] = "HexNAc custom loss"
                    
                    # Y-HexNAc1(-2H2O)-B
                    custom_hexnac_double_water = {
                        **base_hexnac,
                        '_custom_label': 'Y-HexNAc1(-2H2O)-B',
                        '_mass_adjustment': -TWO_WATER_MASS,
                        '_is_oxonium': True,
                        '_custom_type': 'dehydration'
                    }
                    fragments['custom_ions'].append(custom_hexnac_double_water)
                    cleavage_info['custom_ions']['Y-HexNAc1(-2H2O)-B'] = "HexNAc custom loss"
                    break  # Only need one HexNAc1 fragment
            
            # NeuAc1-B custom fragments (dehydration)
            if total_comp.get('NeuAc', 0) > 0:
                b_ions = by_fragments.get('b_ions', [])
                for fragment in b_ions:
                    if (
                        isinstance(fragment, dict)
                        and fragment.get('NeuAc', 0) == 1
                        and _mono_total(fragment) == 1
                    ):
                        custom_neuac_water = {
                            'NeuAc': 1,
                            '_custom_label': 'NeuAc1(-H2O)-B',
                            '_mass_adjustment': -WATER_MASS,
                            '_is_oxonium': True,
                            '_custom_type': 'dehydration'
                        }
                        fragments['custom_ions'].append(custom_neuac_water)
                        cleavage_info['custom_ions']['NeuAc1(-H2O)-B'] = "NeuAc custom loss"
                        break
            
            # Fuc1-B custom fragments (dehydration + demethylation)
            if total_comp.get('Fuc', 0) > 0:
                b_ions = by_fragments.get('b_ions', [])
                for fragment in b_ions:
                    if (
                        isinstance(fragment, dict)
                        and fragment.get('Fuc', 0) == 1
                        and _mono_total(fragment) == 1
                    ):
                        # Fuc1(-H2O)-B
                        custom_fuc_water = {
                            'Fuc': 1,
                            '_custom_label': 'Fuc1(-H2O)-B',
                            '_mass_adjustment': -WATER_MASS,
                            '_is_oxonium': True,
                            '_custom_type': 'dehydration'
                        }
                        fragments['custom_ions'].append(custom_fuc_water)
                        cleavage_info['custom_ions']['Fuc1(-H2O)-B'] = "Fuc custom loss"
                        
                        # Fuc1(-CH3)-B
                        custom_fuc_methyl = {
                            'Fuc': 1,
                            '_custom_label': 'Fuc1(-CH3)-B',
                            '_mass_adjustment': -METHYL_MASS,
                            '_is_oxonium': True,
                            '_custom_type': 'demethylation'
                        }
                        fragments['custom_ions'].append(custom_fuc_methyl)
                        cleavage_info['custom_ions']['Fuc1(-CH3)-B'] = "Fuc custom loss"
                        break
        
        return fragments, cleavage_info
    
    def _count_residues(self, graph: nx.Graph) -> Dict[str, int]:
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
    
    def _format_fragment_string(self, composition: Dict[str, int], ion_type: str, is_y_ion: bool = False) -> str:
        """
        Format a fragment composition as a string.
        
        Args:
            composition: Dict of monosaccharide counts
            ion_type: Type of ion ('B', 'Y', 'BY', 'YY', 'YYY')
            is_y_ion: If True and ion_type is Y-type, append modification-specific suffix
        
        Returns:
            Formatted string like "HexNAc2Hex3-B" or "HexNAc-Redend" for Y-type ions with reduced end
        
        Suffixes for Y-type ions (Y, YY, YYY, etc.):
            - '-Redend': Reduced end (PERMETHYLATED_REDUCED or REDUCED)
            - '-FreeEnd': Free end (PERMETHYLATED_FREE or FREE)
            - '-2AB': 2-AB labeled (TWO_AB or TWO_AB_PERMETHYLATED)
            - '-PEP': Glycopeptide (GLYCOPEPTIDE)
        """
        parts = []
        for mono in ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']:
            count = composition.get(mono, 0)
            if count > 0:
                if count == 1:
                    parts.append(mono)
                else:
                    parts.append(f"{mono}{count}")
        
        if not parts:
            return ""
        
        comp_str = "-".join(parts)
        
        # Add modification-specific suffix for Y-type ions
        if is_y_ion and ion_type.upper() in ('Y', 'YY', 'YYY', 'YYYYY'):
            suffix = self._get_y_ion_suffix()
            return f"{comp_str}{suffix}"
        
        return f"{comp_str}-{ion_type}"
    
    def _count_antennae(self, structure):
        """
        Count the number of antennae (HexNAc chains) in the glycan structure.
        For N-glycans, each HexNAc attached to core branches (nodes 4, 5) is one antenna.
        
        Returns:
            int: Number of antennae
        """
        if self.glycan_type != "N":
            return 0
        
        antenna_count = 0
        core_node = 3
        if core_node in structure:
            for core_branch in structure.successors(core_node):  # nodes 4, 5
                # Count HexNAc children of this core branch
                for child in structure.successors(core_branch):
                    if structure.nodes[child].get('type') == 'HexNAc':
                        antenna_count += 1
        
        return antenna_count
    
    def _get_antenna_branches(self, structure):
        """
        Get antenna-level branch mapping for structures with 3+ antennae.
        Each HexNAc chain attached to core branches is treated as a separate antenna.
        
        Returns dict mapping antenna_root_node -> set of all nodes in that antenna.
        """
        antenna_dict = {}
        
        if self.glycan_type == "N":
            core_node = 3
            if core_node in structure:
                for core_branch in structure.successors(core_node):  # nodes 4, 5
                    # Each HexNAc child of this core branch is an antenna root
                    for hexnac_node in structure.successors(core_branch):
                        if structure.nodes[hexnac_node].get('type') == 'HexNAc':
                            # This HexNAc and all its descendants form one antenna
                            antenna_nodes = set(nx.descendants(structure, hexnac_node)) | {hexnac_node}
                            antenna_dict[hexnac_node] = antenna_nodes
        
        return antenna_dict
    
    def _get_core_branches(self, structure):
        """
        Identify all branches attached to the core mannose.
        For N-glycans, core mannose typically at node 3 with branches at nodes 4, 5, and potentially others.
        Returns dict mapping branch_node -> set of all nodes in that branch.
        """
        branch_dict = {}
        
        # For N-glycans, core mannose (node 3) has direct children that are branch roots
        if self.glycan_type == "N":
            core_node = 3
            if core_node in structure:
                for child in structure.successors(core_node):
                    branch_nodes = set(nx.descendants(structure, child)) | {child}
                    branch_dict[child] = branch_nodes
        
        return branch_dict
    
    def _get_edge_branch(self, u, v, branch_dict):
        """Determine which branch an edge belongs to. Returns branch_node or None."""
        for branch_node, branch_nodes in branch_dict.items():
            if u in branch_nodes or v in branch_nodes:
                return branch_node
        return None

    def _is_valid_yy_composition(self, temp_graph, reducing_component):
        """
        Validate that a YY fragment component is structurally valid by checking
        that all nodes in reducing_component are reachable from the root (node 1)
        in the temporary graph after edge removals.
        
        This prevents topologically impossible compositions like HexNAc3Hex4NeuAc
        where a parent node has been disconnected but child nodes remain.
        
        Args:
            temp_graph: The glycan graph after edge removals
            reducing_component: The set of nodes in the reducing end component
        
        Returns:
            bool: True if all nodes are reachable from root, False otherwise
        """
        root = 1
        
        for node in reducing_component:
            if node == root:
                continue
            
            # Check directed reachability from root in the temporary graph
            try:
                nx.shortest_path(temp_graph, root, node)
            except nx.NetworkXNoPath:
                # Node is not reachable from root - invalid fragment
                return False
        
        return True
    
    def _calculate_composition(self, structure: nx.DiGraph) -> None:
        """
        Calculate composition totals from structure and set instance attributes.
        
        Args:
            structure: NetworkX DiGraph with monosaccharide nodes
        """
        self.hexnac_total = 0
        self.hex_total = 0
        self.fuc_total = 0
        self.neuac_total = 0
        
        for node in structure.nodes():
            node_data = structure.nodes[node]
            # Get monosaccharide type (prefer 'type' over 'monosaccharide' key)
            mono_type = node_data.get('type') or node_data.get('monosaccharide')
            
            if mono_type == 'HexNAc':
                self.hexnac_total += 1
            elif mono_type in ['Hex', 'Man', 'Gal']:  # All hex types
                self.hex_total += 1
            elif mono_type == 'Fuc':
                self.fuc_total += 1
            elif mono_type in ['NeuAc', 'NeuGc']:
                self.neuac_total += 1
