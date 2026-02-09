# glycofrag/structure/predictor.py
"""Glycan structure prediction algorithms."""

from typing import List, Set, Dict, Any, Tuple, Optional
import copy
import time
import networkx as nx

from glycofrag.composition.parser import GlycanComposition
from glycofrag.structure.builder import StructureBuilder
from glycofrag.structure.classifier import StructureClassifier
from glycofrag.core.mass_calculator import GlycanMassCalculator, REDUCING_END_MASSES, ADDITIONAL_MODIFICATIONS
from glycofrag.utils.caching import create_structure_cache_key, get_cached_structures, cache_structures


class StructurePredictor:
    """
    Predicts possible glycan structures from composition using biosynthetic rules.
    
    Responsibility:
    - Generate all possible structures for a given composition
    - Apply biosynthetic rules (high mannose, hybrid, complex)
    - Deduplicate structures using fingerprinting
    - Cache results for performance
    
    Does NOT handle:
    - Composition parsing (delegated to GlycanComposition)
    - Core building (delegated to StructureBuilder)
    - Classification/labeling (delegated to StructureClassifier)
    - Fragmentation

    Public API:
        - predict_structures(...)
    """
    
    def __init__(self, composition: GlycanComposition, glycan_type: str = 'N',
                 max_structures: int = 100, modification_type: int = 0,
                 isomer_sensitive: bool = False, preferred_core: Optional[int] = None):
        """
        Initialize structure predictor.
        
        Args:
            composition: GlycanComposition object
            glycan_type: 'N' or 'O'
            max_structures: Maximum structures to generate
            modification_type: Reducing end modification type
            isomer_sensitive: If True, treats mirror images as distinct
            preferred_core: For O-glycans, preferred core type for visualization (0-8)
                          Only affects visualization, not prediction or fragmentation
        """
        self.composition = composition
        self.glycan_type = glycan_type.upper()
        self.max_structures = max_structures
        self.modification_type = modification_type
        self.isomer_sensitive = isomer_sensitive
        self.preferred_core = preferred_core
        
        # Tracking variables
        self._iteration_count = 0
        self._last_log_time = time.time()
        
        # Results
        self.possible_structures: List[nx.DiGraph] = []
        self.structure_fingerprints: Set[str] = set()
        
        # Helpers
        self.builder = StructureBuilder()
        self.classifier = StructureClassifier()
        self._mass_calculator = GlycanMassCalculator(modification_type=modification_type)
        
        # Calculate core and remaining residues
        if self.glycan_type == "O":
            self.core_hexnac = 1
            self.core_hex = 0
        else:  # N-glycan
            self.core_hexnac = 2
            self.core_hex = 3
        
        self.remaining_hexnac = composition.hexnac - self.core_hexnac
        self.remaining_hex = composition.hex - self.core_hex
        self.remaining_fuc = composition.fuc
        self.remaining_neuac = composition.neuac
        self.remaining_neugc = composition.neugc
    
    def predict_structures(self) -> List[nx.DiGraph]:
        """
        Generate all possible structures for the composition.
        
        Returns:
            List of NetworkX DiGraphs representing glycan structures
        """
        # Create cache key and check cache
        glycan_code_str = f"{self.composition.hexnac}{self.composition.hex}{self.composition.fuc}{self.composition.neuac}"
        cache_key = create_structure_cache_key(
            glycan_code_str,
            self.glycan_type,
            self.max_structures,
            isomer_sensitive=self.isomer_sensitive
        )
        
        cached_result = get_cached_structures(cache_key)
        if cached_result is not None:
            structures, fingerprints = cached_result
            self.possible_structures = [copy.deepcopy(s) for s in structures]
            self.structure_fingerprints = fingerprints.copy()
            
            # Filter by preferred_core if specified (for O-glycans)
            if self.glycan_type == "O" and self.preferred_core is not None:
                filtered_structures = []
                for s in self.possible_structures:
                    possible_cores = s.graph.get('possible_cores', [])
                    if possible_cores and self.preferred_core in possible_cores:
                        # Update selected_core to match preference
                        s.graph['selected_core'] = self.preferred_core
                        filtered_structures.append(s)
                
                if not filtered_structures and self.possible_structures:
                    # No structures match the preferred core - collect available cores
                    available_cores = set()
                    for s in self.possible_structures:
                        cores = s.graph.get('possible_cores', [])
                        available_cores.update(cores)
                    
                    if available_cores:
                        raise ValueError(
                            f"Core {self.preferred_core} was not predicted for this composition.\n"
                            f"Available cores for composition {glycan_code_str}: {sorted(available_cores)}\n"
                            f"Please rerun with preferred_core set to one of the available cores."
                        )
                
                self.possible_structures = filtered_structures
            
            self._attach_mass_to_structures(self.possible_structures)
            print(f"\n{self.glycan_type}-Glycan Structure Prediction Strategy:")
            print(f"Retrieved {len(self.possible_structures)} structures from cache")
            if self.glycan_type == "O" and self.preferred_core:
                print(f"Filtered to Core {self.preferred_core}")
            print("-------------------------------------")
            return self.possible_structures
        
        # Not in cache - proceed with prediction
        print(f"\n{self.glycan_type}-Glycan Structure Prediction Strategy:")
        print(f"Limiting prediction to maximum {self.max_structures} structures")
        print("-------------------------------------")
        
        self.structure_fingerprints.clear()
        self.possible_structures.clear()
        
        # Handle special cases
        if self._handle_special_cases():
            self._attach_mass_to_structures(self.possible_structures)
            cache_structures(cache_key, self.possible_structures, self.structure_fingerprints)
            return self.possible_structures
        
        # Regular structure prediction
        if self.glycan_type == "O":
            self._predict_o_glycan_structures()
        else:
            self._predict_n_glycan_structures()
        
        if len(self.possible_structures) >= self.max_structures:
            print(f"Reached maximum limit of {self.max_structures} structures. Prediction stopped early.")
        print(f"Total unique structures generated: {len(self.possible_structures)}")
        
        # Attach mass and label structures
        self._attach_mass_to_structures(self.possible_structures)
        for structure in self.possible_structures:
            classification = self.classifier.classify_structure(
                structure, self.glycan_type, 
                self.composition.hexnac, self.composition.hex
            )
            self.classifier.label_hex_types(structure, self.glycan_type, classification)
        
        cache_structures(cache_key, self.possible_structures, self.structure_fingerprints)
        return self.possible_structures
    
    def _handle_special_cases(self) -> bool:
        """Handle special simple cases. Returns True if handled."""
        comp = self.composition
        
        # HexNAc(1) - single GlcNAc/GalNAc
        if comp.hexnac == 1 and comp.hex == 0 and comp.fuc == 0 and comp.neuac == 0:
            glycan_desc = "Tn antigen" if self.glycan_type == "O" else "Core GlcNAc fragment"
            print(f"Special case: Simple HexNAc ({glycan_desc})")
            G = nx.DiGraph()
            G.add_node(1, type='HexNAc', position='core_reducing', label='HexNAc')
            self.possible_structures.append(G)
            print(f"Generated 1 fixed structure for HexNAc(1) - {self.glycan_type}-glycan")
            return True
        
        # HexNAc(1)Fuc(1)
        if comp.hexnac == 1 and comp.hex == 0 and comp.fuc == 1 and comp.neuac == 0:
            glycan_desc = "O-glycan with Fuc" if self.glycan_type == "O" else "Core GlcNAc with fucose"
            print(f"Special case: HexNAc(1)Fuc(1) ({glycan_desc})")
            G = nx.DiGraph()
            G.add_node(1, type='HexNAc', position='core_reducing', label='HexNAc')
            G.add_node(2, type='Fuc', position='branch', label='Fuc')
            G.add_edge(1, 2)
            self.possible_structures.append(G)
            print(f"Generated 1 fixed structure for HexNAc(1)Fuc(1) - {self.glycan_type}-glycan")
            return True
        
        # O-glycan HexNAc-Hex (T antigen) - use core builder for proper metadata
        if self.glycan_type == "O" and comp.hexnac == 1 and comp.hex == 1 and comp.fuc == 0 and comp.neuac == 0:
            print("Special case: Simple HexNAc-Hex (T antigen)")
            G = nx.DiGraph()
            G.add_node(1, type='HexNAc', position='core_reducing', label='GalNAc', specific_type='GalNAc')
            G.add_node(2, type='Hex', position='core', label='Gal', specific_type='Gal')
            G.add_edge(1, 2)
            # Add core metadata
            G.graph['possible_cores'] = [1, 8]
            G.graph['default_core'] = 1
            G.graph['core_name'] = 'T antigen (Core 1/8)'
            G.graph['selected_core'] = self.preferred_core if self.preferred_core else 1
            self.possible_structures.append(G)
            print("Generated 1 fixed structure for HexNAc-Hex (1100)")
            return True
        
        return False
    
    def _predict_o_glycan_structures(self):
        """Predict O-glycan structures."""
        core_graphs = self.builder.build_o_glycan_cores(
            self.composition.hexnac, self.composition.hex
        )
        print(f"Core graphs generated: {len(core_graphs)}")
        
        for i, core_data in enumerate(core_graphs):
            if len(self.possible_structures) >= self.max_structures:
                break
            
            # Check if user selected a specific core and if this core matches
            possible_cores = core_data.get('possible_cores', [])
            default_core = core_data.get('default_core')
            core_name = core_data.get('core_name', f'Core {i+1}')
            
            if self.preferred_core is not None:
                # User specified a preferred core - check if it's valid for this structure
                if possible_cores and self.preferred_core not in possible_cores:
                    # Skip this core if user's preference doesn't match
                    continue
            
            print(f"Processing {core_name} (possible cores: {possible_cores or 'N/A'})")
            core_graph = core_data['graph']
            
            # Store core metadata in the graph for visualization
            core_graph.graph['possible_cores'] = possible_cores
            core_graph.graph['default_core'] = default_core
            core_graph.graph['core_name'] = core_name
            core_graph.graph['selected_core'] = self.preferred_core if self.preferred_core else default_core
            
            self.remaining_hexnac = core_data['remaining_hexnac']
            self.remaining_hex = core_data['remaining_hex']
            self.remaining_neuac = self.composition.neuac
            self.remaining_neugc = self.composition.neugc
            self.remaining_fuc = self.composition.fuc
            
            if self.remaining_hexnac >= 0 and self.remaining_hex >= 0:
                print(f"Starting with: HexNAc={self.remaining_hexnac}, Hex={self.remaining_hex}, "
                      f"NeuAc={self.remaining_neuac}, NeuGc={self.remaining_neugc}, Fuc={self.remaining_fuc}")
                
                self._iteration_count = 0
                self._build_complete_o_glycan_structures(copy.deepcopy(core_graph), list(core_graph.nodes()))
                print(f"Structures after {core_name}: {len(self.possible_structures)}")
    
    def _predict_n_glycan_structures(self):
        """Predict N-glycan structures."""
        core_graph = self.builder.build_n_glycan_core()
        self._build_structures(core_graph, [4, 5])
    
    def _structure_fingerprint(self, graph: nx.DiGraph) -> str:
        """
        Create canonical fingerprint for structure deduplication.
        
        Two modes:
        - Default: Deduplicates mirror images (multiset approach)
        - Isomer-sensitive: Treats mirror images as distinct (canonical labeling)
        """
        if self.isomer_sensitive:
            return self._structure_fingerprint_sensitive(graph)
        else:
            return self._structure_fingerprint_default(graph)
    
    def _structure_fingerprint_default(self, graph: nx.DiGraph) -> str:
        """Default: deduplicates mirror images using sorted multiset."""
        root = 1
        
        def get_node_fingerprint(node, parent_seen):
            node_type = graph.nodes[node]['type']
            position = graph.nodes[node]['position']
            
            child_fps = []
            for child in sorted(graph.successors(node)):
                if child not in parent_seen:
                    new_seen = parent_seen | {child}
                    child_fps.append(get_node_fingerprint(child, new_seen))
            
            child_fps.sort()
            return f"({node_type}-{position}[{','.join(child_fps)}])"
        
        return get_node_fingerprint(root, {root})
    
    def _structure_fingerprint_sensitive(self, graph: nx.DiGraph) -> str:
        """Isomer-sensitive: treats mirror images as distinct."""
        root = 1
        canonical_labels = {root: 0}
        next_label = 1
        
        queue = [root]
        visited = {root}
        
        while queue:
            node = queue.pop(0)
            
            children_by_key = {}
            for succ in graph.successors(node):
                node_type = graph.nodes[succ]['type']
                position = graph.nodes[succ]['position']
                key = (node_type, position)
                if key not in children_by_key:
                    children_by_key[key] = []
                children_by_key[key].append(succ)
            
            for (node_type, position), children in sorted(children_by_key.items()):
                children_with_hash = []
                for child in children:
                    subtree_hash = self._get_subtree_hash(graph, child, visited.copy())
                    children_with_hash.append((child, subtree_hash))
                
                children_with_hash.sort(key=lambda x: x[1])
                
                for child, _ in children_with_hash:
                    if child not in visited:
                        canonical_labels[child] = next_label
                        next_label += 1
                        queue.append(child)
                        visited.add(child)
        
        canonical_adj = []
        for orig_node in sorted(canonical_labels.keys(), key=lambda x: canonical_labels[x]):
            node_type = graph.nodes[orig_node]['type']
            position = graph.nodes[orig_node]['position']
            canon_label = canonical_labels[orig_node]
            child_labels = sorted([canonical_labels[succ] for succ in graph.successors(orig_node)])
            canonical_adj.append((canon_label, node_type, position, tuple(child_labels)))
        
        return str(canonical_adj)
    
    def _get_subtree_hash(self, graph: nx.DiGraph, node: int, already_visited: Set[int]) -> str:
        """Generate hash for subtree rooted at node."""
        if node in already_visited:
            return "visited"
        
        node_type = graph.nodes[node]['type']
        position = graph.nodes[node]['position']
        hash_parts = [f"{node_type}-{position}"]
        
        succ_by_type = {}
        for succ in graph.successors(node):
            succ_type = graph.nodes[succ]['type']
            succ_by_type.setdefault(succ_type, []).append(succ)
        
        for succ_type in sorted(succ_by_type.keys()):
            successors = succ_by_type[succ_type]
            sub_hashes = []
            for succ in successors:
                if succ not in already_visited:
                    already_visited.add(succ)
                    sub_hashes.append(self._get_subtree_hash(graph, succ, already_visited))
                    already_visited.remove(succ)
            sub_hashes.sort()
            hash_parts.append(f"{succ_type}:{len(successors)}[{','.join(sub_hashes)}]")
        
        return "(" + "-".join(hash_parts) + ")"
    
    def _build_complete_o_glycan_structures(self, graph: nx.DiGraph, available_nodes: List[int]):
        """Build complete O-glycan structures recursively."""
        self._iteration_count += 1
        
        current_time = time.time()
        if self._iteration_count % 20 == 0 or (current_time - self._last_log_time) > 3:
            print(f"Processing: Structures found: {len(self.possible_structures)}, "
                  f"Remaining: HexNAc={self.remaining_hexnac}, Hex={self.remaining_hex}, "
                  f"Fuc={self.remaining_fuc}, NeuAc={self.remaining_neuac}")
            self._last_log_time = current_time
        
        if self._iteration_count > 200 or len(self.possible_structures) >= self.max_structures:
            return
        
        if (self.remaining_hexnac == 0 and self.remaining_hex == 0 and 
            self.remaining_fuc == 0 and self.remaining_neuac == 0 and self.remaining_neugc == 0):
            fingerprint = self._structure_fingerprint(graph)
            if fingerprint not in self.structure_fingerprints:
                self.structure_fingerprints.add(fingerprint)
                self.possible_structures.append(copy.deepcopy(graph))
            return
        
        valid_parents = []
        for n in graph.nodes():
            if graph.out_degree(n) < 2:
                node_type = graph.nodes[n]['type']
                if node_type == 'NeuAc' and self.remaining_neuac > 0:
                    valid_parents.append(n)
                elif node_type != 'Fuc':
                    valid_parents.append(n)
        
        if self.remaining_hexnac > 0:
            for parent in valid_parents:
                if graph.nodes[parent]['type'] == 'NeuAc':
                    continue
                new_graph = copy.deepcopy(graph)
                new_node_id = max(new_graph.nodes()) + 1
                new_graph.add_node(new_node_id, type='HexNAc', position='branch', label='HexNAc')
                new_graph.add_edge(parent, new_node_id)
                self.remaining_hexnac -= 1
                self._build_complete_o_glycan_structures(new_graph, available_nodes + [new_node_id])
                self.remaining_hexnac += 1
        
        if self.remaining_hex > 0:
            for parent in valid_parents:
                if graph.nodes[parent]['type'] in ['NeuAc', 'Fuc']:
                    continue
                new_graph = copy.deepcopy(graph)
                new_node_id = max(new_graph.nodes()) + 1
                new_graph.add_node(new_node_id, type='Hex', position='branch', label='Hex')
                new_graph.add_edge(parent, new_node_id)
                self.remaining_hex -= 1
                self._build_complete_o_glycan_structures(new_graph, available_nodes + [new_node_id])
                self.remaining_hex += 1
        
        if self.remaining_neuac > 0:
            sialic_counts = {}
            for parent in valid_parents:
                sialic_counts[parent] = sum(1 for child in graph.successors(parent)
                                            if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
            no_sialic_parents = [p for p in valid_parents if sialic_counts[p] == 0]
            neuac_parents = no_sialic_parents if no_sialic_parents else valid_parents
            
            for parent in neuac_parents:
                if graph.nodes[parent]['type'] == 'Fuc':
                    continue
                new_graph = copy.deepcopy(graph)
                new_node_id = max(new_graph.nodes()) + 1
                new_graph.add_node(new_node_id, type='NeuAc', position='terminal', label='NeuAc')
                new_graph.add_edge(parent, new_node_id)
                self.remaining_neuac -= 1
                self._build_complete_o_glycan_structures(new_graph, available_nodes + [new_node_id])
                self.remaining_neuac += 1
        
        if self.remaining_neugc > 0:
            sialic_counts = {}
            for parent in valid_parents:
                sialic_counts[parent] = sum(1 for child in graph.successors(parent)
                                            if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
            no_sialic_parents = [p for p in valid_parents if sialic_counts[p] == 0]
            neugc_parents = no_sialic_parents if no_sialic_parents else valid_parents
            
            for parent in neugc_parents:
                if graph.nodes[parent]['type'] == 'Fuc':
                    continue
                new_graph = copy.deepcopy(graph)
                new_node_id = max(new_graph.nodes()) + 1
                new_graph.add_node(new_node_id, type='NeuGc', position='terminal', label='NeuGc')
                new_graph.add_edge(parent, new_node_id)
                self.remaining_neugc -= 1
                self._build_complete_o_glycan_structures(new_graph, available_nodes + [new_node_id])
                self.remaining_neugc += 1
        
        if self.remaining_fuc > 0:
            for parent in valid_parents:
                if graph.nodes[parent]['type'] in ['NeuAc', 'Fuc']:
                    continue
                new_graph = copy.deepcopy(graph)
                new_node_id = max(new_graph.nodes()) + 1
                new_graph.add_node(new_node_id, type='Fuc', position='terminal', label='Fuc')
                new_graph.add_edge(parent, new_node_id)
                self.remaining_fuc -= 1
                self._build_complete_o_glycan_structures(new_graph, available_nodes + [new_node_id])
                self.remaining_fuc += 1
    
    def _build_structures(self, graph: nx.DiGraph, available_nodes: List[int], 
                         node_id: int = 6, build_stage: int = 0, stage1_start_node_id: int = 6):
        """
        Build N-glycan structures using staged approach.
        
        Stages:
        0: Add all HexNAc to core Hex (4 and 5)
        1: Add all Hex
        2: Add NeuAc
        3: (Reserved)
        4: Add Fuc
        """
        if len(self.possible_structures) >= self.max_structures:
            return
        
        if (self.remaining_hexnac == 0 and self.remaining_hex == 0 and 
            self.remaining_fuc == 0 and self.remaining_neuac == 0 and self.remaining_neugc == 0):
            fingerprint = self._structure_fingerprint(graph)
            if fingerprint not in self.structure_fingerprints:
                self.structure_fingerprints.add(fingerprint)
                self.possible_structures.append(copy.deepcopy(graph))
            return
        
        # Stage 0: Add HexNAc to core Hex nodes (4 and 5)
        if build_stage == 0 and self.remaining_hexnac > 0:
            core_hex_nodes = [n for n in [4, 5] if n in graph.nodes()]
            
            if self.composition.hexnac > 4:
                hexnac_counts = {}
                for node in core_hex_nodes:
                    count = sum(1 for child in graph.successors(node)
                               if graph.nodes[child]['type'] == 'HexNAc')
                    hexnac_counts[node] = count
                
                min_count = min(hexnac_counts.values())
                valid_parents = [node for node in core_hex_nodes 
                               if hexnac_counts[node] == min_count]
            else:
                valid_parents = []
                for node in core_hex_nodes:
                    has_hexnac = any(graph.nodes[succ]['type'] == 'HexNAc' 
                                   for succ in graph.successors(node))
                    if not has_hexnac:
                        valid_parents.append(node)
            
            for parent in valid_parents:
                new_graph = copy.deepcopy(graph)
                new_graph.add_node(node_id, type='HexNAc', position='branch', label='HexNAc')
                new_graph.add_edge(parent, node_id)
                self.remaining_hexnac -= 1
                self._build_structures(new_graph, available_nodes + [node_id], node_id + 1, 0, stage1_start_node_id)
                self.remaining_hexnac += 1
            
            if self.remaining_hexnac == 0:
                self._build_structures(graph, available_nodes, node_id, 1, node_id)
        elif build_stage == 0:
            self._build_structures(graph, available_nodes, node_id, 1, node_id)
        
        # Stage 1: Add Hex
        elif build_stage == 1 and self.remaining_hex > 0:
            if self.composition.hexnac == 2:
                # High Mannose logic
                core_4_children = graph.out_degree(4) if 4 in graph.nodes() else 0
                core_5_children = graph.out_degree(5) if 5 in graph.nodes() else 0
                core_total_children = core_4_children + core_5_children
                
                available_hex = []
                core_available_hex = []
                for n in graph.nodes():
                    if graph.nodes[n]['type'] == 'Hex' and n != 3:
                        if n in [4, 5] and core_total_children < 3:
                            core_available_hex.append(n)
                        elif graph.out_degree(n) < 2:
                            available_hex.append(n)
                
                if core_available_hex:
                    available_hex = core_available_hex
                
                if available_hex:
                    depth_map = dict(nx.shortest_path_length(graph, source=3))
                    min_depth = min(depth_map.get(n, 0) for n in available_hex)
                    available_hex = [n for n in available_hex if depth_map.get(n, 0) == min_depth]
                
                if available_hex:
                    min_children = min(graph.out_degree(n) for n in available_hex)
                    available_hex = [n for n in available_hex if graph.out_degree(n) == min_children]
                
                if available_hex:
                    for parent in sorted(available_hex):
                        new_graph = copy.deepcopy(graph)
                        new_graph.add_node(node_id, type='Hex', position='branch', label='Hex')
                        new_graph.add_edge(parent, node_id)
                        self.remaining_hex -= 1
                        self._build_structures(new_graph, available_nodes + [node_id], node_id + 1, 1, stage1_start_node_id)
                        self.remaining_hex += 1
                else:
                    self._build_structures(graph, available_nodes, node_id, 2, stage1_start_node_id)
                    
            elif self.composition.hexnac == 3 or (self.composition.hex - self.composition.hexnac > 1):
                # Hybrid logic
                branch_hexnac = [n for n in graph.nodes()
                                 if graph.nodes[n]['type'] == 'HexNAc' and
                                    graph.nodes[n]['position'] == 'branch']
                
                antenna_hex = []
                for n in graph.nodes():
                    if graph.nodes[n]['type'] == 'Hex':
                        for ancestor in nx.ancestors(graph, n):
                            if (graph.nodes[ancestor]['type'] == 'HexNAc' and
                                graph.nodes[ancestor]['position'] == 'branch'):
                                antenna_hex.append(n)
                                break
                
                mannose_hex = []
                for n in graph.nodes():
                    if graph.nodes[n]['type'] == 'Hex' and n != 3:
                        is_antenna = False
                        for ancestor in nx.ancestors(graph, n):
                            if (graph.nodes[ancestor]['type'] == 'HexNAc' and
                                graph.nodes[ancestor]['position'] == 'branch'):
                                is_antenna = True
                                break
                        
                        if not is_antenna and not any(graph.nodes[succ]['type'] == 'HexNAc'
                                                      for succ in graph.successors(n)):
                            mannose_hex.append(n)
                
                mannose_hex = [n for n in mannose_hex if graph.out_degree(n) < 2]
                if mannose_hex:
                    depth_map = dict(nx.shortest_path_length(graph, source=3))
                    min_depth = min(depth_map.get(n, 0) for n in mannose_hex)
                    mannose_hex = [n for n in mannose_hex if depth_map.get(n, 0) == min_depth]
                if mannose_hex:
                    min_children = min(graph.out_degree(n) for n in mannose_hex)
                    mannose_hex = [n for n in mannose_hex if graph.out_degree(n) == min_children]
                
                if self.composition.neuac > 0:
                    if len(antenna_hex) < self.composition.neuac:
                        all_targets = sorted(branch_hexnac, key=lambda x: (graph.nodes[x]['position'], x))
                    else:
                        all_targets = sorted(mannose_hex, key=lambda x: (graph.nodes[x]['position'], x))
                else:
                    branch_hexnac = [n for n in branch_hexnac
                                     if not any(graph.nodes[succ]['type'] == 'Hex'
                                                for succ in graph.successors(n))]
                    all_targets = sorted(mannose_hex + branch_hexnac, key=lambda x: (graph.nodes[x]['position'], x))
                
                if all_targets:
                    for parent in all_targets:
                        new_graph = copy.deepcopy(graph)
                        new_graph.add_node(node_id, type='Hex', position='branch', label='Hex')
                        new_graph.add_edge(parent, node_id)
                        self.remaining_hex -= 1
                        self._build_structures(new_graph, available_nodes + [node_id], node_id + 1, 1, stage1_start_node_id)
                        self.remaining_hex += 1
                else:
                    self._build_structures(graph, available_nodes, node_id, 2, stage1_start_node_id)
                    
            else:
                # Complex logic
                branch_hexnac = [n for n in graph.nodes() 
                               if graph.nodes[n]['type'] == 'HexNAc' and 
                                  graph.nodes[n]['position'] == 'branch' and
                                  not any(graph.nodes[succ]['type'] == 'Hex' 
                                         for succ in graph.successors(n))]
                
                if branch_hexnac:
                    for parent in sorted(branch_hexnac):
                        new_graph = copy.deepcopy(graph)
                        new_graph.add_node(node_id, type='Hex', position='branch', label='Hex')
                        new_graph.add_edge(parent, node_id)
                        self.remaining_hex -= 1
                        self._build_structures(new_graph, available_nodes + [node_id], node_id + 1, 1, stage1_start_node_id)
                        self.remaining_hex += 1
                else:
                    self._build_structures(graph, available_nodes, node_id, 2, stage1_start_node_id)
                    
        elif build_stage == 1:
            self._build_structures(graph, available_nodes, node_id, 2, stage1_start_node_id)
        
        # Stage 2: Add NeuAc
        elif build_stage == 2 and self.remaining_neuac > 0:
            branch_hex = [n for n in graph.nodes()
                         if graph.nodes[n]['type'] == 'Hex' and
                            graph.nodes[n]['position'] == 'branch']
            
            if self.composition.hexnac == 3 or (self.composition.hex - self.composition.hexnac > 1):
                antenna_hex = []
                for hex_node in branch_hex:
                    is_antenna = False
                    for ancestor in nx.ancestors(graph, hex_node):
                        if (graph.nodes[ancestor]['type'] == 'HexNAc' and
                            graph.nodes[ancestor]['position'] == 'branch'):
                            is_antenna = True
                            break
                    if is_antenna:
                        antenna_hex.append(hex_node)
                branch_hex = antenna_hex
            
            hex_sialic_counts = {}
            for hex_node in branch_hex:
                sialic_count = sum(1 for child in graph.successors(hex_node)
                                   if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
                hex_sialic_counts[hex_node] = sialic_count
            
            hex_without_sialic = [hex_node for hex_node in branch_hex
                                  if hex_sialic_counts[hex_node] == 0]
            
            if hex_without_sialic:
                min_sialic_count = 0
                valid_parents = [hex_node for hex_node in branch_hex
                                 if hex_sialic_counts[hex_node] == min_sialic_count]
                
                for parent in valid_parents:
                    new_graph = copy.deepcopy(graph)
                    new_graph.add_node(node_id, type='NeuAc', position='terminal', label='NeuAc')
                    new_graph.add_edge(parent, node_id)
                    self.remaining_neuac -= 1
                    self._build_structures(new_graph, available_nodes, node_id + 1, 2, stage1_start_node_id)
                    self.remaining_neuac += 1
            else:
                if not (self.composition.hexnac == 3 or (self.composition.hex - self.composition.hexnac > 1)):
                    branching_hexnac = [n for n in graph.nodes()
                                       if graph.nodes[n]['type'] == 'HexNAc' and
                                          graph.nodes[n].get('position') == 'branch']
                    
                    if branching_hexnac:
                        hexnac_sialic_counts = {}
                        for hexnac_node in branching_hexnac:
                            sialic_count = sum(1 for child in graph.successors(hexnac_node)
                                               if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
                            hexnac_sialic_counts[hexnac_node] = sialic_count
                        hexnac_without_sialic = [n for n in branching_hexnac
                                                 if hexnac_sialic_counts[n] == 0]
                        valid_hexnac = hexnac_without_sialic if hexnac_without_sialic else branching_hexnac
                        
                        for parent in valid_hexnac:
                            new_graph = copy.deepcopy(graph)
                            new_graph.add_node(node_id, type='NeuAc', position='terminal', label='NeuAc')
                            new_graph.add_edge(parent, node_id)
                            self.remaining_neuac -= 1
                            self._build_structures(new_graph, available_nodes, node_id + 1, 2, stage1_start_node_id)
                            self.remaining_neuac += 1
            
            if self.remaining_neuac == 0:
                self._build_structures(graph, available_nodes, node_id, 3, stage1_start_node_id)
        elif build_stage == 2:
            self._build_structures(graph, available_nodes, node_id, 3, stage1_start_node_id)
        
        # Stage 3: Add NeuGc (treat exactly like NeuAc)
        elif build_stage == 3 and self.remaining_neugc > 0:
            branch_hex = [n for n in graph.nodes()
                         if graph.nodes[n]['type'] == 'Hex' and
                            graph.nodes[n]['position'] == 'branch']
            
            if self.composition.hexnac == 3 or (self.composition.hex - self.composition.hexnac > 1):
                antenna_hex = []
                for hex_node in branch_hex:
                    is_antenna = False
                    for ancestor in nx.ancestors(graph, hex_node):
                        if (graph.nodes[ancestor]['type'] == 'HexNAc' and
                            graph.nodes[ancestor]['position'] == 'branch'):
                            is_antenna = True
                            break
                    if is_antenna:
                        antenna_hex.append(hex_node)
                branch_hex = antenna_hex
            
            hex_sialic_counts = {}
            for hex_node in branch_hex:
                sialic_count = sum(1 for child in graph.successors(hex_node)
                                   if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
                hex_sialic_counts[hex_node] = sialic_count
            
            hex_without_sialic = [hex_node for hex_node in branch_hex
                                  if hex_sialic_counts[hex_node] == 0]
            
            if hex_without_sialic:
                min_sialic_count = 0
                valid_parents = [hex_node for hex_node in branch_hex
                                 if hex_sialic_counts[hex_node] == min_sialic_count]
                
                for parent in valid_parents:
                    new_graph = copy.deepcopy(graph)
                    new_graph.add_node(node_id, type='NeuGc', position='terminal', label='NeuGc')
                    new_graph.add_edge(parent, node_id)
                    self.remaining_neugc -= 1
                    self._build_structures(new_graph, available_nodes, node_id + 1, 3, stage1_start_node_id)
                    self.remaining_neugc += 1
            else:
                if not (self.composition.hexnac == 3 or (self.composition.hex - self.composition.hexnac > 1)):
                    branching_hexnac = [n for n in graph.nodes()
                                       if graph.nodes[n]['type'] == 'HexNAc' and
                                          graph.nodes[n].get('position') == 'branch']
                    
                    if branching_hexnac:
                        hexnac_sialic_counts = {}
                        for hexnac_node in branching_hexnac:
                            sialic_count = sum(1 for child in graph.successors(hexnac_node)
                                               if graph.nodes[child]['type'] in ('NeuAc', 'NeuGc'))
                            hexnac_sialic_counts[hexnac_node] = sialic_count
                        hexnac_without_sialic = [n for n in branching_hexnac
                                                 if hexnac_sialic_counts[n] == 0]
                        valid_hexnac = hexnac_without_sialic if hexnac_without_sialic else branching_hexnac
                        
                        for parent in valid_hexnac:
                            new_graph = copy.deepcopy(graph)
                            new_graph.add_node(node_id, type='NeuGc', position='terminal', label='NeuGc')
                            new_graph.add_edge(parent, node_id)
                            self.remaining_neugc -= 1
                            self._build_structures(new_graph, available_nodes, node_id + 1, 3, stage1_start_node_id)
                            self.remaining_neugc += 1
            
            if self.remaining_neugc == 0:
                self._build_structures(graph, available_nodes, node_id, 4, stage1_start_node_id)
        elif build_stage == 3:
            self._build_structures(graph, available_nodes, node_id, 4, stage1_start_node_id)
        
        # Stage 4: Add Fuc
        elif build_stage == 4 and self.remaining_fuc > 0:
            hexnac_nodes = [n for n in graph.nodes()
                          if graph.nodes[n]['type'] == 'HexNAc' and
                             graph.nodes[n].get('position') in ['core_reducing', 'branch']]
            
            # Only attach Fuc to HexNAc nodes that don't already have a Fuc child
            for parent in hexnac_nodes:
                # Check if this HexNAc already has a Fucose attached
                fuc_children = [child for child in graph.successors(parent)
                               if graph.nodes[child]['type'] == 'Fuc']
                
                # Only attach Fuc if this HexNAc doesn't already have one
                if not fuc_children:
                    new_graph = copy.deepcopy(graph)
                    new_graph.add_node(node_id, type='Fuc', position='terminal', label='Fuc')
                    new_graph.add_edge(parent, node_id)
                    self.remaining_fuc -= 1
                    self._build_structures(new_graph, available_nodes, node_id + 1, 4, stage1_start_node_id)
                    self.remaining_fuc += 1
        elif build_stage == 4:
            if (self.remaining_hexnac == 0 and self.remaining_hex == 0 and 
                self.remaining_neuac == 0 and self.remaining_neugc == 0 and self.remaining_fuc == 0):
                fingerprint = self._structure_fingerprint(graph)
                if fingerprint not in self.structure_fingerprints:
                    self.structure_fingerprints.add(fingerprint)
                    self.possible_structures.append(copy.deepcopy(graph))
    
    def _get_tree_structure_debug(self, structure: nx.DiGraph) -> str:
        """Generate a debug string showing the tree structure of a glycan."""
        lines = []
        
        # Find root node
        root = None
        for node in structure.nodes():
            if structure.in_degree(node) == 0:
                root = node
                break
        
        if root is None:
            return "[No root found]"
        
        def format_node(node, depth=0):
            node_data = structure.nodes[node]
            node_type = node_data.get('type', 'Unknown')
            prefix = "  " * depth + "|- "
            fuc_children = [child for child in structure.successors(node)
                           if structure.nodes[child]['type'] == 'Fuc']
            fuc_info = f" [HAS {len(fuc_children)} Fuc]" if fuc_children else ""
            result = [f"{prefix}{node_type}(id={node}){fuc_info}"]
            
            for child in structure.successors(node):
                result.extend(format_node(child, depth + 1))
            return result
        
        return "\n".join(format_node(root))
    
    def _attach_mass_to_structures(self, structures: List[nx.DiGraph]) -> None:
        """Attach mass and breakdown metadata to each structure."""
        for structure in structures:
            composition = self._count_residues(structure)
            breakdown, mass = self._format_mass_breakdown(composition)
            structure.graph['mass'] = mass
            structure.graph['mass_breakdown'] = breakdown
    
    def _count_residues(self, graph: nx.DiGraph) -> Dict[str, int]:
        """Count residues in a structure."""
        counts = {'HexNAc': 0, 'Hex': 0, 'Fuc': 0, 'NeuAc': 0}
        for node in graph.nodes():
            node_type = graph.nodes[node]['type']
            if node_type in counts:
                counts[node_type] += 1
            elif node_type in ['Man', 'Gal']:
                counts['Hex'] += 1
        return counts
    
    def _format_mass_breakdown(self, composition: Dict[str, int]) -> Tuple[str, float]:
        """Format mass breakdown and calculate total mass."""
        hexnac_count = composition.get('HexNAc', 0)
        hex_count = composition.get('Hex', 0)
        fuc_count = composition.get('Fuc', 0)
        neuac_count = composition.get('NeuAc', 0)
        
        hexnac_mass = self._mass_calculator.MONO_MASSES['HexNAc']
        hex_mass = self._mass_calculator.MONO_MASSES['Hex']
        fuc_mass = self._mass_calculator.MONO_MASSES['Fuc']
        neuac_mass = self._mass_calculator.MONO_MASSES['NeuAc']
        reducing_end_mass = REDUCING_END_MASSES.get(self.modification_type, 0.0)
        additional_mass = ADDITIONAL_MODIFICATIONS.get(self.modification_type, 0.0)
        
        total_mass = (
            hexnac_count * hexnac_mass +
            hex_count * hex_mass +
            fuc_count * fuc_mass +
            neuac_count * neuac_mass +
            reducing_end_mass +
            additional_mass
        )
        
        breakdown = (
            f"{hexnac_count}*{hexnac_mass:.4f} + "
            f"{hex_count}*{hex_mass:.4f} + "
            f"{fuc_count}*{fuc_mass:.4f} + "
            f"{neuac_count}*{neuac_mass:.4f} + "
            f"RE({reducing_end_mass:.4f}) + "
            f"ADD({additional_mass:.4f})"
        )
        
        return breakdown, total_mass
