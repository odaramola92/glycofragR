# glycofrag/structure/builder.py
"""Core glycan structure building - creates N-glycan and O-glycan cores."""

from typing import List, Dict, Any
import networkx as nx


class StructureBuilder:
    """
    Builds core glycan structures (N-glycan and O-glycan cores).
    
    Responsibility: Create the foundational core structures that serve as
    starting points for structure prediction algorithms.
    
    Does NOT handle:
    - Structure prediction (adding remaining residues)
    - Fragmentation
    - Classification

    Public API:
        - build_n_glycan_core()
        - build_o_glycan_cores(...)
    """
    
    @staticmethod
    def build_n_glycan_core() -> nx.DiGraph:
        """
        Build the N-glycan core structure: HexNAc2Hex3.
        
        The core structure consists of:
        - Node 1: HexNAc (reducing end, attached to peptide)
        - Node 2: HexNAc (core)
        - Node 3: Hex (central mannose)
        - Nodes 4, 5: Hex (branching mannoses, alpha-3 and alpha-6)
        
        Returns:
            NetworkX DiGraph representing the core structure
        
        Example:
            >>> builder = StructureBuilder()
            >>> core = builder.build_n_glycan_core()
            >>> print(f"Nodes: {len(core.nodes())}")
            Nodes: 5
            >>> print(f"Edges: {len(core.edges())}")
            Edges: 4
        """
        G = nx.DiGraph()
        
        # Add core structure (bottom to top)
        G.add_node(1, type='HexNAc', position='core_reducing', label='HexNAc')
        G.add_node(2, type='HexNAc', position='core', label='HexNAc')
        G.add_edge(1, 2)
        
        # First branching Hex (central mannose)
        G.add_node(3, type='Hex', position='core_central', label='Hex')
        G.add_edge(2, 3)
        
        # Left and right core Hex (alpha-3 and alpha-6 mannose)
        G.add_node(4, type='Hex', position='core_branch', label='Hex')
        G.add_node(5, type='Hex', position='core_branch', label='Hex')
        G.add_edge(3, 4)
        G.add_edge(3, 5)
        
        return G
    
    @staticmethod
    def build_o_glycan_cores(hexnac_total: int, hex_total: int) -> List[Dict[str, Any]]:
        """
        Build all possible O-glycan core structures based on composition.
        
        O-glycans can have multiple core types:
        - Core 0: Single HexNAc (Tn antigen)
        - Core 1: HexNAc-Hex (GalNAc-Gal, T antigen)
        - Core 2: HexNAc with two branching HexNAc
        - Core 3: HexNAc with single branching HexNAc
        - Single Hex core
        
        Args:
            hexnac_total: Total number of HexNAc residues in composition
            hex_total: Total number of Hex residues in composition
        
        Returns:
            List of dictionaries, each containing:
            - 'graph': NetworkX DiGraph of the core structure
            - 'core_hexnac': Number of HexNAc used in this core
            - 'core_hex': Number of Hex used in this core
            - 'remaining_hexnac': HexNAc remaining after core
            - 'remaining_hex': Hex remaining after core
        
        Example:
            >>> builder = StructureBuilder()
            >>> cores = builder.build_o_glycan_cores(hexnac_total=2, hex_total=1)
            >>> print(f"Generated {len(cores)} possible cores")
            >>> for i, core in enumerate(cores):
            ...     print(f"Core {i+1}: HexNAc={core['core_hexnac']}, Hex={core['core_hex']}")
        """
        core_structures = []
        
        # Core 1: HexNAc-Hex (GalNAc-Gal, T antigen)
        # Also represents Core 8 (same fragments, different linkage)
        if hex_total >= 1 and hexnac_total >= 1:
            G1 = nx.DiGraph()
            G1.add_node(1, type='HexNAc', position='core_reducing', label='GalNAc', specific_type='GalNAc')
            G1.add_node(2, type='Hex', position='core', label='Gal', specific_type='Gal')
            G1.add_edge(1, 2)
            core_structures.append({
                'graph': G1,
                'core_hexnac': 1,
                'core_hex': 1,
                'remaining_hexnac': hexnac_total - 1,
                'remaining_hex': hex_total - 1,
                'possible_cores': [1, 8],
                'default_core': 1,
                'core_name': 'T antigen (Core 1/8)'
            })
        
        # Core 2/4: GalNAc with branching (β1,6) GlcNAc
        # Core 2: GalNAc-Gal-GlcNAc(β1,6)
        # Core 4: GalNAc-Gal with two GlcNAc branches (β1,3 and β1,6)
        # Same fragments, different structures
        if hexnac_total >= 3:
            G2 = nx.DiGraph()
            G2.add_node(1, type='HexNAc', position='core_reducing', label='GalNAc', specific_type='GalNAc')
            G2.add_node(2, type='HexNAc', position='core_branch', label='GlcNAc', specific_type='GlcNAc')
            G2.add_node(3, type='HexNAc', position='core_branch', label='GlcNAc', specific_type='GlcNAc')
            G2.add_edge(1, 2)
            G2.add_edge(1, 3)
            core_structures.append({
                'graph': G2,
                'core_hexnac': 3,
                'core_hex': 0,
                'remaining_hexnac': hexnac_total - 3,
                'remaining_hex': hex_total,
                'possible_cores': [2, 4],
                'default_core': 2,
                'core_name': 'Branched (Core 2/4)'
            })
        
        # Core 3/5/6/7: Extended cores (structural isomers, same fragments)
        # Core 3: GlcNAc(β1,3)-GalNAc
        # Core 5: GalNAc(α1,3)-GalNAc  
        # Core 6: GlcNAc(β1,6)-GalNAc
        # Core 7: GalNAc(α1,6)-GalNAc
        if hexnac_total >= 2:
            G3 = nx.DiGraph()
            G3.add_node(1, type='HexNAc', position='core_reducing', label='GalNAc', specific_type='GalNAc')
            G3.add_node(2, type='HexNAc', position='core_branch', label='GalNAc', specific_type='GalNAc')  # Default to GalNAc
            G3.add_edge(1, 2)
            core_structures.append({
                'graph': G3,
                'core_hexnac': 2,
                'core_hex': 0,
                'remaining_hexnac': hexnac_total - 2,
                'remaining_hex': hex_total,
                'possible_cores': [3, 5, 6, 7],
                'default_core': 3,
                'core_name': 'Extended (Core 3/5/6/7)'
            })
        
        # NOTE: Removed "Gal-only" core (single Hex starting core) 
        # O-glycans MUST have GalNAc at the reducing end, so Gal-only is invalid
        
        # Core 0: Single HexNAc (Tn antigen) - ONLY when there's no Hex in composition
        # If hex_total > 0, Core 0 would just add Hex back, creating duplicates with Core 1/8
        if hexnac_total >= 1 and hex_total == 0:
            G0 = nx.DiGraph()
            G0.add_node(1, type='HexNAc', position='core_reducing', label='GalNAc', specific_type='GalNAc')
            core_structures.append({
                'graph': G0,
                'core_hexnac': 1,
                'core_hex': 0,
                'remaining_hexnac': hexnac_total - 1,
                'remaining_hex': hex_total,
                'possible_cores': [0],
                'default_core': 0,
                'core_name': 'Tn antigen (Core 0)'
            })
        
        return core_structures
