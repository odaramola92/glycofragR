"""
Glycan structure visualization module.

Provides Matplotlib-based visualization of glycan structures with hierarchical tree layout.
"""

import os
import tempfile
from typing import Optional, List, Union, Dict
import matplotlib
import matplotlib.figure
import networkx as nx


class GlycanVisualizer:
    """
    Visualize glycan structures with SNFG-compliant symbols and colors.
    
    **SNFG Symbol Guide:**
    
    ========== ======== =========
    Symbol     Color    Residue
    ========== ======== =========
    Square     Blue     GlcNAc
    Square     Yellow   GalNAc
    Circle     Green    Mannose
    Circle     Yellow   Galactose
    Diamond    Magenta  NeuAc
    Diamond    Silver   NeuGc
    Triangle   Red      Fucose
    ========== ======== =========
    
    **Quick Start:**
    
        >>> from glycofrag import Glycan, Glycopeptide, GlycanVisualizer
        >>>
        >>> # From a Glycan object
        >>> glycan = Glycan('4501', glycan_type='N')
        >>> structures = glycan.predict_structures()
        >>> GlycanVisualizer.draw(structures, 'all', glycan_code='4501')
        >>>
        >>> # From a Glycopeptide object
        >>> gp = Glycopeptide('EEQYNSTYR', '4501', 5, 'N')
        >>> GlycanVisualizer.draw(gp.glycan_structures, 'all', glycan_code='4501')
    
    **Methods:**
    
    - ``draw(structures, structures_to_draw=1, ...)`` — Main entry point (int, list, or 'all')
    - ``visualize(structure, ...)`` — Draw a single structure graph
    - ``visualize_structure(structures, structure_number=1, ...)`` — Draw by 1-based index
    
    **Advanced Visualization Control:**
    
    All visualization methods accept optional kwargs for fine-tuning layout:
    
    - ``vertical_gap``: Vertical spacing between levels (float or dict per level)
    - ``horizontal_spacing``: Horizontal spacing between siblings (float or dict per level)
    - ``node_size``: Monosaccharide marker sizes (float or dict per type)
    - ``show_node_numbers``: Display node connectivity numbers (bool, default False)
    - ``output_path``: Path to save PNG file (str)
    - ``show``: Display plot interactively (bool, default True)
    - ``title``: Figure title (str, default 'Glycan Structure')
    - ``figsize``: Figure size in inches (tuple, default (5, 4))
    """
    
    # Color scheme for monosaccharides
    NODE_COLORS = {
        'HexNAc': '#0066FF',      # Blue (generic HexNAc)
        'GlcNAc': '#0066FF',      # Blue (GlcNAc specifically)
        'GalNAc': '#FFFF00',      # Yellow (GalNAc specifically for O-glycans)
        'Man': '#00A651',          # Green (mannose - core Hex)
        'Gal': '#FFFF00',          # Yellow (galactose - branching Hex)
        'Hex': '#00A651',          # Green (fallback for unlabeled Hex)
        'NeuAc': '#FF00FF',        # Magenta (diamond)
        'NeuGc': '#C0C0C0',        # Silver (diamond)
        'Fuc': '#FF0000',          # Red (triangle)
    }
    
    NODE_SHAPES = {
        'HexNAc': 'box',
        'GlcNAc': 'box',
        'GalNAc': 'box',
        'Man': 'circle',
        'Gal': 'circle',
        'Hex': 'circle',
        'NeuAc': 'diamond',
        'NeuGc': 'diamond',
        'Fuc': 'triangle',
    }

    # Visualization sizing controls (edit these to change geometry)
    NODE_MARKER_SIZE = 200  # scatter size in points^2 (constant on screen)
    NODE_MARKER_SIZE_NONCIRCLE = 170  # slightly smaller for box/diamond/triangle
    NODE_MARKER_SIZE_NEUAC = 120  # smaller size specifically for NeuAc diamonds
    NODE_EDGE_LINEWIDTH = 0.5  # thinner outline for shapes
    EDGE_LINEWIDTH = 2
    LEGEND_MARKER_SIZE = 10  # in points

    # Layout spacing controls (data units)
    VERTICAL_GAP = 1.2
    HORIZONTAL_SPACING = 1.5
    FUC_HORIZONTAL_OFFSET = 1.0  # horizontal distance for fucose on edge branches
    FUC_MIDDLE_HORIZONTAL_OFFSET = 0.4  # horizontal distance for fucose on middle branches (shorter)
    FUC_VERTICAL_OFFSET = -0.9  # upward bend for fucose (negative = upward)
    CORE_BRANCH_BASE_SPACING = 0.9  # base spacing between nodes 4 and 5
    CORE_BRANCH_SPACING_INCREMENT = 1.0  # additional spacing per branch (increased for Hybrid clarity)
    
    @staticmethod
    def _get_vertical_gap(depth: int, custom_gaps: Optional['Union[float, Dict[int, float]]']) -> float:
        """
        Get vertical gap for a specific depth level.
        
        Args:
            depth: Current depth (0=root, 1=first level, etc.)
            custom_gaps: Custom vertical gaps (float for global, dict for per-level)
        
        Returns:
            Vertical gap to apply at this depth
        """
        if custom_gaps is None:
            return GlycanVisualizer.VERTICAL_GAP
        if isinstance(custom_gaps, dict):
            return custom_gaps.get(depth, GlycanVisualizer.VERTICAL_GAP)
        return custom_gaps  # Single value override
    
    @staticmethod
    def _get_horizontal_spacing(depth: int, custom_spacings: Optional['Union[float, Dict[int, float]]']) -> float:
        """
        Get horizontal spacing for a specific depth level.
        
        Args:
            depth: Current depth (0=root, 1=first level, etc.)
            custom_spacings: Custom horizontal spacings (float for global, dict for per-level)
        
        Returns:
            Horizontal spacing to apply at this depth
        """
        if custom_spacings is None:
            return GlycanVisualizer.HORIZONTAL_SPACING
        if isinstance(custom_spacings, dict):
            return custom_spacings.get(depth, GlycanVisualizer.HORIZONTAL_SPACING)
        return custom_spacings
    
    @staticmethod
    def _get_node_size(node_type: str, custom_sizes: Optional['Union[float, Dict[str, float]]']) -> float:
        """
        Get node size for a specific monosaccharide type.
        
        Args:
            node_type: Monosaccharide type (HexNAc, Gal, Man, NeuAc, Fuc, etc.)
            custom_sizes: Custom node sizes (float for global, dict for per-type)
        
        Returns:
            Node size (marker size in points^2) to apply
        """
        if custom_sizes is None:
            # Return current size based on type
            if node_type in ('NeuAc', 'NeuGc'):
                return GlycanVisualizer.NODE_MARKER_SIZE_NEUAC
            shape = GlycanVisualizer.NODE_SHAPES.get(node_type, 'circle')
            if shape == 'circle':
                return GlycanVisualizer.NODE_MARKER_SIZE
            else:
                return GlycanVisualizer.NODE_MARKER_SIZE_NONCIRCLE
        if isinstance(custom_sizes, dict):
            return custom_sizes.get(node_type, GlycanVisualizer._get_node_size(node_type, None))
        return custom_sizes  # Single value override
    
    @staticmethod
    def visualize_with_matplotlib(structure: 'nx.DiGraph',
                                 title: str = "Glycan Structure",
                                 figsize: tuple = (5, 4),
                                 show: bool = True,
                                 output_path: Optional[str] = None,
                                 show_node_numbers: bool = False,
                                 vertical_gap: Optional['Union[float, Dict[int, float]]'] = None,
                                 horizontal_spacing: Optional['Union[float, Dict[int, float]]'] = None,
                                 node_size: Optional['Union[float, Dict[str, float]]'] = None) -> Optional[str]:
        """
        Visualize glycan structure using Matplotlib with hierarchical tree layout.
        
        Creates publication-quality PNG images with SNFG-compliant monosaccharide symbols
        and intelligent layout for clear structure visualization.
        
        **Visualization Control:**
        
        This method provides fine-grained control over the visual appearance through
        three optional parameters. Each parameter accepts either:
        
        1. **Single float** - Apply uniformly to all levels/nodes
        2. **Dictionary** - Apply specific values per level (depth) or node type
        3. **None** (default) - Use built-in defaults (no change to current behavior)
        
        Args:
            structure: NetworkX DiGraph representing glycan structure
            title: Figure title (appears above visualization)
            figsize: Tuple of (width, height) in inches (default: 5x4)
                    **IMPORTANT:** If you provide a custom figsize, automatic figure scaling
                    based on vertical_gap/horizontal_spacing is DISABLED. Your exact dimensions
                    will be used. This allows you to override the auto-scaling behavior.
                    
                    Examples:
                    - figsize=(5, 4): Default, auto-scaling applies if spacing parameters provided
                    - figsize=(8, 10): Your exact size, auto-scaling disabled
            show: Whether to display in matplotlib window (default: True)
            output_path: Optional path to save PNG file. If provided, saves without showing.
            show_node_numbers: Show node connectivity numbers inside monosaccharides.
                             Default: False (clean, publication-ready).
                             Set to True to display numbering (1, 2, 3...) for debugging
                             or structure analysis.
            vertical_gap: Control vertical spacing between levels.
                         - None (default): Use default spacing (1.2)
                         - Float (e.g., 1.5): Apply uniformly to all levels
                         - Dict (e.g., {0: 1.0, 1: 1.5, 4: 2.0}): Specify per depth level
                         
                         **Level Definition:** Depth from root (reducing end):
                         Level 0 = Core reducing HexNAc
                         Level 1 = Core HexNAc
                         Level 2 = Central Man
                         Level 3 = Branch Man
                         Level 4 = Branch HexNAc
                         Level 5 = Gal/terminal sugars
                         Level 6 = NeuAc/NeuGc
            
            horizontal_spacing: Control horizontal spacing between siblings at same level.
                               - None (default): Use default spacing (1.5)
                               - Float (e.g., 2.0): Apply uniformly to all levels
                               - Dict (e.g., {3: 2.0, 4: 2.5}): Specify per depth level
                               
                               **Use Case:** Increase spacing at branching levels to prevent overlap
            
            node_size: Control size of monosaccharide markers (in points^2).
                      - None (default): Use default sizes (200 for circles, 170 for shapes, 120 for sialic acids)
                      - Float (e.g., 250): Apply uniformly to all nodes
                      - Dict (e.g., {'HexNAc': 300, 'Gal': 200, 'Fuc': 150}): Specify per monosaccharide type
                      
                      **Monosaccharide Types:** HexNAc, GlcNAc, GalNAc, Man, Gal, Hex, NeuAc, NeuGc, Fuc
        
        Returns:
            Matplotlib figure object if successful, None if matplotlib unavailable
        
        Example - Default Behavior:
            >>> from glycofrag import Glycan
            >>> from glycofrag.io.visualizer import GlycanVisualizer
            >>> 
            >>> glycan = Glycan('4501', glycan_type='N')
            >>> structures = glycan.predict_structures()
            >>> 
            >>> # Use defaults (no change to current behavior)
            >>> GlycanVisualizer.visualize_with_matplotlib(structures[0])
        
        Example - Global Overrides:
            >>> # Increase all spacing and node sizes uniformly
            >>> GlycanVisualizer.visualize_with_matplotlib(
            ...     structures[0],
            ...     vertical_gap=1.8,          # More vertical space everywhere
            ...     horizontal_spacing=2.5,    # More horizontal space everywhere
            ...     node_size=300              # Larger nodes everywhere
            ... )
        
        Example - Per-Level Control:
            >>> # Fine-tune specific levels
            >>> GlycanVisualizer.visualize_with_matplotlib(
            ...     structures[0],
            ...     vertical_gap={0: 1.0, 1: 1.2, 2: 1.5, 3: 2.0},  # Gradual increase
            ...     horizontal_spacing={4: 3.0, 5: 3.5},            # Wide branching area
            ...     node_size={'HexNAc': 350, 'Gal': 250, 'Fuc': 180}  # Emphasize HexNAc
            ... )
        
        Example - Mixed Control:
            >>> # Global default + specific overrides
            >>> GlycanVisualizer.visualize_with_matplotlib(
            ...     structures[0],
            ...     vertical_gap=1.5,                    # Global vertical spacing
            ...     horizontal_spacing={3: 2.5, 4: 3.0}, # Only widen branch areas
            ...     node_size=250                        # Global node size
            ... )
        
        Example - With GlycanAnalysis:
            >>> from glycofrag import GlycanAnalysis
            >>> 
            >>> analysis = GlycanAnalysis(
            ...     glycan_code='6623',
            ...     glycan_type='N',
            ...     modification_type='permethylated_reduced'
            ... )
            >>> 
            >>> # Visualize with custom spacing (passed through to visualizer)
            >>> analysis.visualize_structure(
            ...     1,  # Structure index
            ...     vertical_gap={4: 2.0, 5: 2.5},
            ...     horizontal_spacing=2.5,
            ...     node_size={'NeuAc': 180}
            ... )
        
        Notes:
            - All positioning is in data space units (not pixels) for consistent scaling
            - Colors follow SNFG standard (Standard Notation for Glycans)
            - Legend automatically positioned on left side
            - High DPI (200) for publication-quality output
            - Custom parameters do NOT affect structure prediction, only visualization
        """
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Circle, Rectangle, RegularPolygon
            import math
        except ImportError:
            print("[WARN] Matplotlib not installed. Install with: pip install matplotlib")
            return None
        
        # Auto-scale figure size based on custom spacing parameters
        # BUT: if user explicitly provides custom figsize, respect it and skip auto-scaling
        adjusted_figsize = figsize
        user_provided_custom_figsize = (figsize != (5, 4))  # Check if user changed default
        
        if not user_provided_custom_figsize:
            # User is using default figsize, so apply auto-scaling
            adjusted_figsize = list(figsize)
            
            # Height scaling based on vertical_gap (very conservative: cube root scaling)
            if vertical_gap is not None:
                if isinstance(vertical_gap, dict):
                    avg_gap = sum(vertical_gap.values()) / len(vertical_gap) if vertical_gap else 1.2
                else:
                    avg_gap = vertical_gap
                
                # Very conservative scaling: use power of 0.35 for gentler growth
                # Examples: gap=1.2 (default) → scale 1.0x
                #           gap=2.4 (2x)      → scale 1.26x (was 1.41x)
                #           gap=4.8 (4x)      → scale 1.57x (was 2.0x)
                #           gap=9.6 (8x)      → scale 1.96x (was 2.5x capped)
                import math
                ratio = avg_gap / 1.2
                scale_factor = math.pow(max(1.0, ratio), 0.35)
                # Cap maximum scaling at 1.8x to keep images reasonable
                scale_factor = min(scale_factor, 1.8)
                adjusted_figsize[1] = figsize[1] * scale_factor
            
            # Width scaling based on horizontal_spacing (very conservative)
            if horizontal_spacing is not None:
                if isinstance(horizontal_spacing, dict):
                    avg_spacing = sum(horizontal_spacing.values()) / len(horizontal_spacing) if horizontal_spacing else 1.5
                else:
                    avg_spacing = horizontal_spacing
                
                import math
                ratio = avg_spacing / 1.5
                scale_factor = math.pow(max(1.0, ratio), 0.35)
                scale_factor = min(scale_factor, 1.8)
                adjusted_figsize[0] = figsize[0] * scale_factor
            
            adjusted_figsize = tuple(adjusted_figsize)
        
        # Create figure with two subplots: legend (left, smaller) and structure (right, larger)
        fig = plt.figure(figsize=adjusted_figsize)
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 3], wspace=0.3)
        
        ax_legend = fig.add_subplot(gs[0])
        ax_structure = fig.add_subplot(gs[1])
        
        # Build a tree layout (root = 1) for clear glycan structure visualization
        root = 1
        positions = {}
        
        # Capture custom parameters in closure so assign_positions can access them
        custom_vertical_gaps = vertical_gap
        custom_horizontal_spacings = horizontal_spacing
        # Track if user provided custom spacing (to respect their preferences)
        user_provided_horizontal = horizontal_spacing is not None

        def assign_positions(node, x_pos, y_pos, depth_level, parent_x_pos=None, sibling_positions=None):
            """
            Assign positions to nodes in tree layout.
            
            Args:
                node: Current node ID
                x_pos: X position of current node
                y_pos: Y position (vertical coordinate)
                depth_level: Depth level (0=root, 1=first level, etc.)
                parent_x_pos: X position of parent (for context)
                sibling_positions: Dict of {sibling_node: x_position} for context
            """
            children = list(structure.successors(node))
            
            # Separate Fuc from other children
            fuc_children = [c for c in children if structure.nodes[c].get('type') == 'Fuc']
            other_children = [c for c in children if structure.nodes[c].get('type') != 'Fuc']
            num_spread = len(other_children)  # Define here so it's available for Fuc positioning
            
            # Get custom vertical gap for this depth level
            v_gap = GlycanVisualizer._get_vertical_gap(depth_level, custom_vertical_gaps)
            # Get custom horizontal spacing for this depth level
            h_spacing_default = GlycanVisualizer._get_horizontal_spacing(depth_level, custom_horizontal_spacings)
            
            if not children:
                # Leaf node
                positions[node] = (x_pos, y_pos)
            elif len(other_children) == 1 and not fuc_children:
                # Single non-Fuc child: align vertically
                child = other_children[0]
                assign_positions(child, x_pos, y_pos + v_gap, depth_level + 1, parent_x_pos=x_pos)
                positions[node] = (x_pos, y_pos)
            else:
                # Multiple children: spread them out
                if other_children:
                    spacing = h_spacing_default
                    
                    # Only apply automatic spacing adjustments if user didn't provide custom horizontal spacing
                    if not user_provided_horizontal:
                        # Increase spacing for branching HexNAc to prevent visual overlap
                        if num_spread > 2:
                            spacing = h_spacing_default * 1.5
                        else:
                            spacing = h_spacing_default

                        # Extra spacing for HexNAc branching into multiple Hex (Gal) children
                        if structure.nodes[node].get('type') == 'HexNAc' and num_spread >= 2:
                            hex_children = [c for c in other_children if structure.nodes[c].get('type') == 'Hex']
                            if len(hex_children) >= 2:
                                spacing = max(spacing, h_spacing_default * 1.5)
                        
                        # Widen spacing specifically between core branches (node 4 and 5) to avoid crowding
                        if node == 3 and set(other_children) == {4, 5}:
                            # Count total grandchildren to determine spacing
                            total_grandchildren = 0
                            for child in other_children:
                                grandchildren = list(structure.successors(child))
                                non_fuc_grandchildren = [g for g in grandchildren if structure.nodes[g].get('type') != 'Fuc']
                                total_grandchildren += len(non_fuc_grandchildren)

                                # If a grandchild is a branch HexNAc, count its Hex children as additional branches
                                for gc in non_fuc_grandchildren:
                                    if structure.nodes[gc].get('type') == 'HexNAc':
                                        hex_children = [h for h in structure.successors(gc)
                                                        if structure.nodes[h].get('type') == 'Hex']
                                        total_grandchildren += len(hex_children)
                            
                            # Calculate spacing: base + (increment * num_branches)
                            dynamic_spacing = GlycanVisualizer.CORE_BRANCH_BASE_SPACING + (total_grandchildren - 1) * GlycanVisualizer.CORE_BRANCH_SPACING_INCREMENT
                            spacing = max(spacing, dynamic_spacing)
                    
                    # For better branch separation, use deeper recursion positioning
                    start_x = x_pos - (spacing * (num_spread - 1) / 2)
                    child_pos_map = {}
                    for i, child in enumerate(sorted(other_children)):
                        child_x = start_x + (i * spacing)
                        child_pos_map[child] = child_x
                        # Don't inherit spacing - let children calculate their own
                        assign_positions(child, child_x, y_pos + v_gap, depth_level + 1,
                                       parent_x_pos=x_pos, sibling_positions=child_pos_map)
                else:
                    positions[node] = (x_pos, y_pos)
                
                # Center node between children
                child_xs = [positions[c][0] for c in other_children if c in positions]
                if child_xs:
                    positions[node] = ((min(child_xs) + max(child_xs)) / 2, y_pos)
                else:
                    positions[node] = (x_pos, y_pos)
        
        # First pass: position all non-Fuc nodes
        assign_positions(root, 0, 0, 0, parent_x_pos=None, sibling_positions={root: 0})
        
        # Second pass: position all Fuc nodes after HexNAcs are known
        def position_fuc_nodes():
            """Position fucose nodes using complete HexNAc position information."""
            for node in structure.nodes():
                if structure.nodes[node].get('type') != 'Fuc':
                    continue
                
                # Get parent node
                parents = list(structure.predecessors(node))
                if not parents:
                    continue
                parent = parents[0]
                parent_x, parent_y = positions.get(parent, (0, 0))
                parent_depth = parent_y
                
                # Get all Fuc siblings (other Fuc attached to same parent)
                fuc_siblings = [n for n in structure.successors(parent) 
                               if structure.nodes[n].get('type') == 'Fuc']
                fuc_idx = fuc_siblings.index(node) if node in fuc_siblings else 0
                
                # Spread multiple Fuc vertically
                fuc_vertical_spacing = 0.5
                fuc_y = parent_depth + (fuc_idx * fuc_vertical_spacing)
                horizontal_offset = GlycanVisualizer.FUC_HORIZONTAL_OFFSET
                
                # Determine parent's branch side
                side = "right" if parent_x > 0.01 else "left" if parent_x < -0.01 else "center"
                
                # Find ALL HexNAc positions on parent's side (now fully positioned)
                side_hexnac_x = []
                for n, attrs in structure.nodes(data=True):
                    if attrs.get("type") == "HexNAc" and n in positions:
                        node_x = positions[n][0]
                        if side == "right" and node_x > 0.01:
                            side_hexnac_x.append(node_x)
                        elif side == "left" and node_x < -0.01:
                            side_hexnac_x.append(node_x)
                
                if side in ("right", "left") and side_hexnac_x:
                    if side == "right":
                        inner_x = min(side_hexnac_x)  # Closest to center (x=0)
                        outer_x = max(side_hexnac_x)  # Furthest from center
                        is_inner = abs(parent_x - inner_x) <= abs(parent_x - outer_x)
                        fuc_x = parent_x - horizontal_offset if is_inner else parent_x + horizontal_offset
                    else:  # left
                        inner_x = max(side_hexnac_x)  # Closest to center (x=0)
                        outer_x = min(side_hexnac_x)  # Furthest from center
                        is_inner = abs(parent_x - inner_x) <= abs(parent_x - outer_x)
                        fuc_x = parent_x + horizontal_offset if is_inner else parent_x - horizontal_offset
                else:
                    # Fallback
                    fuc_x = parent_x + horizontal_offset
                
                positions[node] = (fuc_x, fuc_y)
        
        position_fuc_nodes()

        # Helper function to calculate edge endpoints accounting for node sizes
        def get_edge_endpoints(parent_pos, child_pos, parent_node_attrs, child_node_attrs):
            """Calculate edge start/end points accounting for node shapes and sizes."""
            import math

            x0, y0 = parent_pos
            x1, y1 = child_pos

            # Convert data coords to display coords (pixels)
            p0_disp = ax_structure.transData.transform((x0, y0))
            p1_disp = ax_structure.transData.transform((x1, y1))

            dx = p1_disp[0] - p0_disp[0]
            dy = p1_disp[1] - p0_disp[1]
            dist = math.hypot(dx, dy)
            if dist == 0:
                return (x0, y0), (x1, y1)

            # Unit direction vector in display coords
            nx = dx / dist
            ny = dy / dist

            def marker_shrink_points(node_attrs, ux, uy):
                node_type = node_attrs.get('type', 'Unknown')
                shape = GlycanVisualizer.NODE_SHAPES.get(node_type, 'circle')
                # Use custom node size if provided
                size = GlycanVisualizer._get_node_size(node_type, node_size)

                # Direction-aware shrink distance in points based on marker geometry
                abs_ux = abs(ux)
                abs_uy = abs(uy)
                eps = 1e-9

                if shape == 'circle':
                    # area = size (points^2): r = sqrt(area / pi)
                    return math.sqrt(size / math.pi)

                if shape == 'box':
                    # area = size => side = sqrt(size), half-width = sqrt(size)/2
                    half_width = math.sqrt(size) / 2.0
                    return half_width / max(abs_ux, abs_uy, eps)

                if shape == 'diamond':
                    # diamond is rotated square: vertices at distance s/sqrt(2)
                    r_vertex = math.sqrt(size) / math.sqrt(2.0)
                    return r_vertex / max(abs_ux + abs_uy, eps)

                if shape == 'triangle':
                    # equilateral triangle: area = (sqrt(3)/4) * s^2
                    side = math.sqrt((4.0 * size) / math.sqrt(3.0))
                    circumradius = side / math.sqrt(3.0)
                    # Slightly reduce to reach edge more often
                    return circumradius * 0.85

                # Fallback
                return math.sqrt(size / math.pi)

            # Convert points to pixels (72 pts per inch)
            dpi = ax_structure.figure.dpi
            parent_r_px = marker_shrink_points(parent_node_attrs, nx, ny) * dpi / 72.0
            child_r_px = marker_shrink_points(child_node_attrs, nx, ny) * dpi / 72.0

            # Start/end points in display coords adjusted by marker radius
            start_disp = (p0_disp[0] + nx * parent_r_px, p0_disp[1] + ny * parent_r_px)
            end_disp = (p1_disp[0] - nx * child_r_px, p1_disp[1] - ny * child_r_px)

            # Convert back to data coords
            start_data = ax_structure.transData.inverted().transform(start_disp)
            end_data = ax_structure.transData.inverted().transform(end_disp)

            return (start_data[0], start_data[1]), (end_data[0], end_data[1])

        # Draw edges on structure axis
        for parent, child in structure.edges():
            parent_pos = positions[parent]
            child_pos = positions[child]
            parent_attrs = structure.nodes[parent]
            child_attrs = structure.nodes[child]
            
            (x0, y0), (x1, y1) = get_edge_endpoints(parent_pos, child_pos, parent_attrs, child_attrs)
            ax_structure.plot([x0, x1], [y0, y1], color='gray', linewidth=GlycanVisualizer.EDGE_LINEWIDTH, zorder=1)

        # Draw nodes with SNFG-like shapes (constant size in points)
        marker_map = {
            'circle': 'o',
            'box': 's',
            'diamond': 'D',
            'triangle': '<'  # Default: tip points left (toward parent connection)
        }
        for node, attrs in structure.nodes(data=True):
            x, y = positions[node]
            # Use specific_type if available (GalNAc/GlcNAc), otherwise use general type
            node_type = attrs.get('specific_type', attrs.get('type', 'Unknown'))
            color = GlycanVisualizer.NODE_COLORS.get(node_type, '#CCCCCC')
            # For shapes, use general type (HexNAc, Hex, etc.)
            shape_type = attrs.get('type', 'Unknown')
            shape = GlycanVisualizer.NODE_SHAPES.get(shape_type, 'circle')
            marker = marker_map.get(shape, 'o')
            
            # For Fuc triangles, always point toward the attachment line (parent)
            if shape == 'triangle':
                parents = list(structure.predecessors(node))
                if parents:
                    parent = parents[0]
                    parent_x = positions[parent][0]
                    node_x = x
                    marker = '<' if node_x > parent_x else '>'
                else:
                    marker = '<'  # Default if no parent found
            
            # Use custom node size if provided
            size = GlycanVisualizer._get_node_size(shape_type, node_size)

            ax_structure.scatter(
                [x], [y],
                s=size,
                marker=marker,
                c=color,
                edgecolors='black',
                linewidths=GlycanVisualizer.NODE_EDGE_LINEWIDTH,
                zorder=2
            )
            
            # Conditionally display node numbers based on show_node_numbers parameter
            if show_node_numbers:
                ax_structure.text(x, y, str(node), ha='center', va='center', fontsize=4, fontweight='bold', zorder=3)

        # Configure structure axis with fixed data limits
        ax_structure.set_title(title, fontsize=14, fontweight='bold')
        ax_structure.set_aspect('equal', adjustable='box')
        
        # Calculate data bounds
        all_x = [pos[0] for pos in positions.values()]
        all_y = [pos[1] for pos in positions.values()]
        
        # Add padding around the structure (20% of range or minimum 2 units)
        x_range = max(all_x) - min(all_x)
        y_range = max(all_y) - min(all_y)
        x_padding = max(x_range * 0.2, 2.0)
        y_padding = max(y_range * 0.2, 2.0)
        
        # Set fixed data limits so distances remain constant
        ax_structure.set_xlim(min(all_x) - x_padding, max(all_x) + x_padding)
        ax_structure.set_ylim(min(all_y) - y_padding, max(all_y) + y_padding)
        ax_structure.axis('off')

        # Legend on left axis
        from matplotlib.lines import Line2D
        legend_handles = []
        legend_labels = []
        marker_map = {
            'circle': 'o',
            'box': 's',
            'diamond': 'D',
            'triangle': '<'  # Tip points toward connection
        }
        # Legend order: GlcNAc, GalNAc, Man, Gal, NeuAc, NeuGc, Fuc
        legend_order = ['GlcNAc', 'GalNAc', 'Man', 'Gal', 'NeuAc', 'NeuGc', 'Fuc']
        for label in legend_order:
            if label in GlycanVisualizer.NODE_SHAPES:
                shape = GlycanVisualizer.NODE_SHAPES[label]
                color = GlycanVisualizer.NODE_COLORS.get(label, '#CCCCCC')
                marker = marker_map.get(shape, 'o')
                handle = Line2D([0], [0], marker=marker, color='none',
                                markerfacecolor=color, markeredgecolor='black',
                                markersize=GlycanVisualizer.LEGEND_MARKER_SIZE, linestyle='None')
                legend_handles.append(handle)
                legend_labels.append(label)

        ax_legend.legend(legend_handles, legend_labels, loc='center', fontsize=10,
                        frameon=True, edgecolor='black', fancybox=False)
        ax_legend.axis('off')
        
        plt.tight_layout()
        
        if output_path:
            try:
                fig.savefig(output_path, dpi=200, bbox_inches='tight')
                print(f"[PASS] Matplotlib image saved to: {output_path}")
                return output_path
            except Exception as e:
                print(f"[WARN] Could not save matplotlib image: {e}")
                return None

        if show:
            plt.show()
        
        return None
    

    @staticmethod
    def visualize(structure: 'nx.DiGraph',
                 **kwargs) -> Optional[str]:
        """
        Visualize a glycan structure using Matplotlib.
        
        Args:
            structure: NetworkX DiGraph representing glycan structure
            **kwargs: Additional arguments passed to visualize_with_matplotlib:
                - title (str): Figure title. Default: 'Glycan Structure'
                - figsize (tuple): Figure size (width, height). Default: (5, 4)
                - show (bool): Whether to display the plot. Default: True
                - output_path (str): Optional path to save PNG file
        
        Returns:
            Path to output file if output_path provided, None otherwise
        
        Example:
            >>> from glycofrag import Glycan
            >>> glycan = Glycan('5604', glycan_type='N')
            >>> structures = glycan.predict_structures()
            >>> GlycanVisualizer.visualize(structures[0], output_path='glycan.png')
        """
        return GlycanVisualizer.visualize_with_matplotlib(structure, **kwargs)

    @staticmethod
    def visualize_structure(structures: List['nx.DiGraph'],
                           structure_number: int = 1,
                           **kwargs) -> Optional[str]:
        """
        Visualize a specific structure from a list of structures using 1-based indexing.
        
        Args:
            structures: List of NetworkX DiGraph objects
            structure_number: Structure number to visualize (1, 2, 3, ...). Default: 1
            **kwargs: Advanced visualization options:
                - vertical_gap (float | dict): Vertical spacing between levels
                  * float: Apply uniformly (e.g., 1.8)
                  * dict: Per-level control (e.g., {3: 2.0, 4: 2.5})
                - horizontal_spacing (float | dict): Horizontal spacing between siblings
                  * float: Apply uniformly (e.g., 2.5)
                  * dict: Per-level control (e.g., {4: 3.0, 5: 3.5})
                - node_size (float | dict): Monosaccharide marker sizes
                  * float: Apply uniformly (e.g., 300)
                  * dict: Per-type control (e.g., {'HexNAc': 350, 'NeuAc': 180})
                - show_node_numbers (bool): Display node numbers (default: False)
                - output_path (str): Path to save PNG file
                - show (bool): Display plot interactively (default: True)
                - title (str): Figure title (default: 'Glycan Structure')
                - figsize (tuple): Figure size in inches (default: (5, 4))
        
        Returns:
            Path to output file, or None if visualization failed
            
        Example:
            >>> GlycanVisualizer.visualize_structure(
            ...     structures,
            ...     structure_number=1,
            ...     vertical_gap={3: 2.0, 4: 2.5},
            ...     horizontal_spacing=2.5,
            ...     node_size={'HexNAc': 350, 'NeuAc': 180},
            ...     output_path='glycan.png'
            ... )
        """
        structure_index = structure_number - 1
        
        if structure_index < 0:
            print(f"[ERROR] Structure number must be >= 1. Got: {structure_number}")
            return None
        
        if structure_index >= len(structures):
            print(f"[ERROR] Structure {structure_number} not found. Valid range: 1 to {len(structures)}")
            return None
        
        structure = structures[structure_index]
        print(f"Visualizing Structure {structure_number} of {len(structures)}...")
        
        return GlycanVisualizer.visualize(structure, **kwargs)

    @staticmethod
    def draw(structures: List['nx.DiGraph'],
             structures_to_draw=1,
             glycan_code: str = "",
             save_images: bool = True,
             show: bool = False,
             output_dir: str = "",
             **kwargs) -> List[str]:
        """
        Draw glycan structures with a simple interface.
        
        Args:
            structures: List of NetworkX DiGraph objects (e.g. gp.glycan_structures)
            structures_to_draw: Which structures to visualize:
                - 1       → draw structure 1 only (default)
                - 'all'   → draw all structures
                - [1,3,5] → draw structures 1, 3, and 5
            glycan_code: Glycan code for filenames (optional)
            save_images: Save PNG files (default True)
            show: Display plots (default False)
            output_dir: Directory for saved images (default: current directory)
            **kwargs: Additional args passed to visualize:
                - show_node_numbers (bool): Show node IDs inside symbols (default: False).
                  Set True for debugging/analysis.
                - figsize (tuple): Figure size in inches, e.g. (6, 5)
                - title (str): Custom figure title
        
        Returns:
            List of saved file paths
        
        Examples:
            >>> GlycanVisualizer.draw(gp.glycan_structures)                           # structure 1
            >>> GlycanVisualizer.draw(gp.glycan_structures, 'all')                    # all structures
            >>> GlycanVisualizer.draw(gp.glycan_structures, [1, 3, 5])                # specific ones
            >>> GlycanVisualizer.draw(gp.glycan_structures, 'all', glycan_code='4501')
            >>> GlycanVisualizer.draw(gp.glycan_structures, 'all', show_node_numbers=True)  # with numbering
        """
        import os
        
        if not structures:
            print("[ERROR] No structures to draw")
            return []
        
        total = len(structures)
        
        # Normalize to list of 1-based ints
        if isinstance(structures_to_draw, str) and structures_to_draw.lower() == 'all':
            nums = list(range(1, total + 1))
        elif isinstance(structures_to_draw, int):
            nums = [structures_to_draw]
        elif isinstance(structures_to_draw, (list, tuple, range)):
            nums = [int(n) for n in structures_to_draw]
        else:
            print(f"[ERROR] Invalid structures_to_draw: {structures_to_draw}")
            return []
        
        # Filter valid
        nums = [n for n in nums if 1 <= n <= total]
        if not nums:
            print(f"[ERROR] No valid structure numbers (total: {total})")
            return []
        
        output_files = []
        for n in nums:
            out_path = None
            if save_images:
                label = glycan_code or 'glycan'
                filename = f"{label}_structure_{n:03d}.png"
                out_path = os.path.join(output_dir or os.getcwd(), filename)
            
            title = kwargs.pop('title', None) or f"{glycan_code or 'Glycan'} - Structure {n}/{total}"
            
            result = GlycanVisualizer.visualize_structure(
                structures, structure_number=n,
                title=title, show=show, output_path=out_path, **kwargs
            )
            if result:
                output_files.append(result)
            
            # Restore title for next iteration
            kwargs.pop('title', None)
        
        if output_files:
            print(f"[SAVED] {len(output_files)} image(s)")
        
        return output_files


def visualize_structures(structures: List['nx.DiGraph'],
                        indices: Optional[List[int]] = None) -> None:
    """
    Convenience function to visualize multiple structures.
    
    Args:
        structures: List of NetworkX DiGraph objects
        indices: List of structure indices to visualize. If None, visualizes all.
    
    Example:
        >>> from glycofrag import Glycan
        >>> from glycofrag.io.visualizer import visualize_structures
        >>> glycan = Glycan('5604', glycan_type='N')
        >>> structures = glycan.predict_structures()
        >>> visualize_structures(structures)
    """
    if indices is None:
        indices = list(range(len(structures)))
    
    for idx in indices:
        if 0 <= idx < len(structures):
            print(f"\n{'='*60}")
            print(f"Visualizing Structure {idx + 1}")
            print(f"{'='*60}")
            GlycanVisualizer.visualize(structures[idx])
        else:
            print(f"[WARN] Structure index {idx} out of range")
