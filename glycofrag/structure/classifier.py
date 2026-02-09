# glycofrag/structure/classifier.py
"""Glycan structure classification and labeling."""

from typing import Dict
import networkx as nx


class StructureClassifier:
    """
    Classifies glycan structures and labels monosaccharide types.
    
    Responsibility:
    - Determine glycan class (High Mannose, Hybrid, Complex, O-glycan types)
    - Label Hex nodes as Man (mannose) or Gal (galactose) based on biosynthetic rules
    
    Does NOT handle:
    - Structure prediction
    - Fragmentation
    - Mass calculation

    Public API:
        - classify_structure(...)
        - label_hex_types(...)
    """
    
    @staticmethod
    def classify_structure(graph: nx.DiGraph, glycan_type: str, hexnac_total: int, 
                          hex_total: int) -> str:
        """
        Classify glycan structure based on composition and type.
        
        Args:
            graph: NetworkX DiGraph representing the glycan structure
            glycan_type: 'N' or 'O'
            hexnac_total: Total number of HexNAc in composition
            hex_total: Total number of Hex in composition
        
        Returns:
            String classification:
            - N-glycans: 'High Mannose', 'Hybrid', 'Complex', 'Unknown'
            - O-glycans: 'O-GalNAc (Tn/STn)', 'O-GalNAc-Gal (T/ST)', etc.
        
        Example:
            >>> classifier = StructureClassifier()
            >>> # Assume we have a structure graph
            >>> classification = classifier.classify_structure(graph, 'N', 4, 5)
            >>> print(classification)
            Complex
        """
        if glycan_type == "O":
            # Determine O-glycan subtypes
            if hex_total == 0 and hexnac_total >= 2:
                return "O-GalNAc (Tn/STn)"
            elif hex_total >= 1 and hexnac_total >= 1:
                return "O-GalNAc-Gal (T/ST)"
            elif hexnac_total == 3 and hex_total == 0:
                return "O-GalNAc with HexNAc branches"
            else:
                return "Complex O-glycan"
        else:
            # N-glycan classification
            if hexnac_total == 2:
                return "High Mannose"
            elif hexnac_total == 3:
                return "Hybrid"
            elif hexnac_total > 3:
                return "Complex"
            else:
                return "Unknown"
    
    @staticmethod
    def label_hex_types(graph: nx.DiGraph, glycan_type: str, classification: str) -> None:
        """
        Label Hex nodes as Man (mannose, core) or Gal (galactose, branching).
        
        This is critical for accurate glycan structure annotation. Labeling rules
        depend on glycan type and biosynthetic classification:
        
        **N-glycan High Mannose (HexNAc = 2):**
        - All Hex nodes labeled as Man (mannose)
        - No branching HexNAc present - structure remains core-only
        
        **N-glycan Hybrid (HexNAc = 3):**
        - Core Hex (nodes 3, 4, 5): Man (mannose arm backbone)
        - Branching HexNAc detected via position='branch' attribute
        - Hex attached to branching HexNAc: Gal (antenna arm)
        - NeuAc only attaches to antenna Gal, never mannose arm
        
        **N-glycan Complex (HexNAc > 3):**
        - Core Hex (nodes 3, 4, 5): Man (trimannosyl core)
        - All other Hex nodes: Gal (branching antenna structures)
        
        **O-glycan (all types):**
        - All Hex labeled as Gal (galactose)
        - No mannose backbone
        
        Args:
            graph: NetworkX DiGraph representing the glycan structure (modified in-place)
            glycan_type: 'N' or 'O'
            classification: Classification string from classify_structure()
        
        Returns:
            None (modifies graph in-place, updating node 'type' and 'label' attributes)
        
        Example:
            >>> classifier = StructureClassifier()
            >>> # Assume we have a structure graph
            >>> classifier.label_hex_types(graph, 'N', 'Complex')
            >>> # Now graph nodes have specific Man/Gal labels
        """
        core_hex_nodes = {3, 4, 5}
        
        if glycan_type == "O":
            # O-glycans: all Hex are Gal (cannot have mannose)
            for node, attrs in graph.nodes(data=True):
                if attrs.get('type') == 'Hex':
                    attrs['type'] = 'Gal'
                    attrs['label'] = 'Gal'
        
        elif "High Mannose" in classification:
            # High Mannose: all Hex are Man (no branching HexNAc)
            for node, attrs in graph.nodes(data=True):
                if attrs.get('type') == 'Hex':
                    attrs['type'] = 'Man'
                    attrs['label'] = 'Man'
        
        elif "Hybrid" in classification:
            # Hybrid: core Hex (3,4,5) are Man, others on branching HexNAc are Gal
            branching_hexnac = [n for n in graph.nodes()
                               if graph.nodes[n]['type'] == 'HexNAc' and
                                  graph.nodes[n].get('position') == 'branch']
            
            branching_hex_nodes = set()
            for parent in branching_hexnac:
                for child in graph.successors(parent):
                    if graph.nodes[child]['type'] == 'Hex':
                        branching_hex_nodes.add(child)
            
            for node, attrs in graph.nodes(data=True):
                if attrs.get('type') == 'Hex':
                    if node in core_hex_nodes:
                        attrs['type'] = 'Man'
                        attrs['label'] = 'Man'
                    elif node in branching_hex_nodes:
                        attrs['type'] = 'Gal'
                        attrs['label'] = 'Gal'
                    else:
                        # Hex not on branching HexNAc (default to Man for Hybrid)
                        attrs['type'] = 'Man'
                        attrs['label'] = 'Man'
        
        elif "Complex" in classification:
            # Complex: nodes 3,4,5 are Man, all other Hex are Gal
            for node, attrs in graph.nodes(data=True):
                if attrs.get('type') == 'Hex':
                    if node in core_hex_nodes:
                        attrs['type'] = 'Man'
                        attrs['label'] = 'Man'
                    else:
                        attrs['type'] = 'Gal'
                        attrs['label'] = 'Gal'
        
        else:
            # Unknown/other: default to Complex logic (core are Man, others are Gal)
            for node, attrs in graph.nodes(data=True):
                if attrs.get('type') == 'Hex':
                    if node in core_hex_nodes:
                        attrs['type'] = 'Man'
                        attrs['label'] = 'Man'
                    else:
                        attrs['type'] = 'Gal'
                        attrs['label'] = 'Gal'
