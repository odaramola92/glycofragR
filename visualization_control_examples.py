"""
Visualization Control Examples for GlycoFrag

This script demonstrates how to control glycan structure visualization at different levels:
1. Using GlycanVisualizer directly (most control)
2. Using GlycanAnalysis (released glycan workflow)
3. Using GlycoPeptideAnalysis (glycopeptide workflow)
4. Using Glycan facade (intermediate API)

Author: Oluwatosin Daramola
Date: February 8, 2026
"""

from glycofrag import Glycan, GlycanAnalysis, GlycoPeptideAnalysis
from glycofrag.io.visualizer import GlycanVisualizer
import matplotlib.pyplot as plt


def example1_direct_visualizer():
    """
    Example 1: Direct control using GlycanVisualizer
    
    This provides the most fine-grained control over every aspect of visualization.
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Direct GlycanVisualizer Control")
    print("="*80)
    
    # Generate structures
    glycan = Glycan('6623', glycan_type='N')
    structures = glycan.predict_structures()
    print(f"Generated {len(structures)} structures for glycan code 6623")
    
    # Example 1a: Default visualization (no changes)
    print("\n1a. Default visualization (baseline):")
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Default Visualization",
        output_path="viz_example1a_default.png",
        show=False
    )
    print("   Saved: viz_example1a_default.png")
    
    # Example 1b: Global spacing increase
    print("\n1b. Global spacing increase:")
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Global Spacing Increase",
        vertical_gap=1.8,           # More vertical space
        horizontal_spacing=2.5,     # More horizontal space
        node_size=300,              # Larger nodes
        output_path="viz_example1b_global.png",
        show=False
    )
    print("   Saved: viz_example1b_global.png")
    
    # Example 1c: Per-level control (professional publication style)
    print("\n1c. Per-level control for publication:")
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Publication-Quality Layout",
        vertical_gap={
            0: 1.0,   # Level 0: Core reducing HexNAc
            1: 1.2,   # Level 1: Core HexNAc
            2: 1.5,   # Level 2: Central Man
            3: 2.0,   # Level 3: Branch Man (wider)
            4: 1.8,   # Level 4: Branch HexNAc
            5: 1.5    # Level 5: Terminal sugars
        },
        horizontal_spacing={
            3: 2.5,   # Widen at branch Man level
            4: 3.0    # Extra wide at branching HexNAc level
        },
        node_size={'HexNAc': 350, 'Gal': 280, 'NeuAc': 180, 'Fuc': 200},
        figsize=(8, 6),
        output_path="viz_example1c_publication.png",
        show=False
    )
    print("   Saved: viz_example1c_publication.png")
    
    # Example 1d: Compact layout for presentations
    print("\n1d. Compact layout for presentations:")
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Compact Presentation Style",
        vertical_gap=0.8,           # Tighter vertical spacing
        horizontal_spacing=1.2,     # Tighter horizontal spacing
        node_size=180,              # Smaller nodes
        figsize=(4, 3),
        output_path="viz_example1d_compact.png",
        show=False
    )
    print("   Saved: viz_example1d_compact.png")


def example2_glycan_analysis():
    """
    Example 2: Using GlycanAnalysis (released glycan workflow)
    
    GlycanAnalysis provides high-level MS/MS analysis. Visualization parameters
    pass through to the underlying GlycanVisualizer.
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: GlycanAnalysis with Custom Visualization")
    print("="*80)
    
    # Initialize analysis (this generates structures + fragments automatically)
    analysis = GlycanAnalysis(
        glycan_code='4501',
        glycan_type='N',
        modification_type='permethylated_reduced',
        fragment_types=['BY']
    )
    print(f"Analysis complete: {len(analysis.structures)} structures predicted")
    print(f"Theoretical fragments: {len(analysis.theoretical_df)} rows")
    
    # Example 2a: Visualize with custom spacing
    print("\n2a. Visualize structure 1 with custom spacing:")
    analysis.visualize_structure(
        1,  # Structure number (1-indexed)
        vertical_gap=1.5,
        horizontal_spacing=2.0,
        node_size=250,
        output_image='viz_example2a_analysis.png'
    )
    print("   Saved: viz_example2a_analysis.png")
    
    # Example 2b: Emphasize specific monosaccharides
    print("\n2b. Emphasize HexNAc and NeuAc:")
    analysis.visualize_structure(
        1,
        node_size={
            'HexNAc': 400,  # Large
            'Man': 250,      # Medium
            'Gal': 250,      # Medium
            'NeuAc': 200,    # Larger than default
            'Fuc': 180       # Medium-small
        },
        vertical_gap=1.6,
        output_image='viz_example2b_emphasis.png'
    )
    print("   Saved: viz_example2b_emphasis.png")


def example3_glycopeptide_analysis():
    """
    Example 3: Using GlycoPeptideAnalysis (glycopeptide workflow)
    
    GlycoPeptideAnalysis handles glycopeptide fragmentation. Visualization parameters
    also pass through to GlycanVisualizer.
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: GlycoPeptideAnalysis with Custom Visualization")
    print("="*80)
    
    # Initialize glycopeptide analysis
    gp_analysis = GlycoPeptideAnalysis(
        peptide_sequence='EEQYNSTYR',
        glycan_code='4501',
        glycosylation_site=5,
        glycan_type='N',
        modification_type='free',
        peptide_fragment_types=['by'],
        glycan_fragment_types=['BY']
    )
    print(f"Glycopeptide analysis complete: {len(gp_analysis.glycan_structures)} glycan structures")
    print(f"Theoretical fragments: {len(gp_analysis.theoretical_df)} rows")
    
    # Example 3a: Visualize with wide horizontal spacing (for complex glycans)
    print("\n3a. Wide horizontal spacing for clarity:")
    gp_analysis.visualize_glycan_structure(
        structure_number=1,
        vertical_gap=1.4,
        horizontal_spacing={4: 3.5, 5: 4.0},  # Very wide at branching levels
        node_size=280,
        output_image='viz_example3a_glycopeptide.png'
    )
    print("   Saved: viz_example3a_glycopeptide.png")
    
    # Example 3b: Compact for supplementary figures
    print("\n3b. Compact style for supplementary figures:")
    gp_analysis.visualize_glycan_structure(
        structure_number=1,
        vertical_gap=0.9,
        horizontal_spacing=1.3,
        node_size=160,
        output_image='viz_example3b_supplementary.png'
    )
    print("   Saved: viz_example3b_supplementary.png")


def example4_glycan_facade():
    """
    Example 4: Using Glycan facade (intermediate API)
    
    The Glycan class provides a mid-level API. You can predict structures and
    then visualize them with custom parameters.
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Glycan Facade with Custom Visualization")
    print("="*80)
    
    # Create glycan and predict structures
    glycan = Glycan('5502', glycan_type='N')
    structures = glycan.predict_structures()
    print(f"Generated {len(structures)} structures for glycan code 5502")
    
    # Example 4a: Loop through all structures with consistent styling
    print("\n4a. Batch visualization with consistent styling:")
    for i, structure in enumerate(structures[:3], 1):  # First 3 structures
        GlycanVisualizer.visualize_with_matplotlib(
            structure,
            title=f"Structure {i} of {len(structures)}",
            vertical_gap=1.5,
            horizontal_spacing=2.2,
            node_size=270,
            output_path=f"viz_example4a_structure_{i}.png",
            show=False
        )
        print(f"   Saved: viz_example4a_structure_{i}.png")
    
    # Example 4b: Compare two visualization styles side-by-side
    print("\n4b. Comparison: default vs custom:")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    plt.suptitle("Visualization Style Comparison")
    
    # Default on left
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Default Style",
        show=False,
        output_path="viz_example4b_default_only.png"
    )
    
    # Custom on right
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Custom Style",
        vertical_gap=1.8,
        horizontal_spacing=2.5,
        node_size={'HexNAc': 320, 'Man': 280, 'Gal': 280, 'NeuAc': 180},
        show=False,
        output_path="viz_example4b_custom_only.png"
    )
    print("   Saved: viz_example4b_default_only.png")
    print("   Saved: viz_example4b_custom_only.png")


def example5_level_definitions():
    """
    Example 5: Demonstrate level (depth) definitions
    
    Shows what each level number represents in the glycan tree.
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: Understanding Level Definitions")
    print("="*80)
    
    print("""
Level Definition (depth from root, the reducing end):
------------------------------------------------------
Level 0: Core reducing HexNAc (node 1)
Level 1: Core HexNAc (node 2)
Level 2: Central Man (node 3)
Level 3: Branch Man (nodes 4, 5)
Level 4: Branch HexNAc (nodes 6, 7, 8, 9)
Level 5: Gal/terminal sugars (nodes 10, 11, 12)
Level 6: NeuAc/NeuGc (nodes 13, 14, 15)
Level 7+: Additional terminal modifications

Example glycan code 6623 (tri-antennary sialylated):
    """)
    
    glycan = Glycan('6623', glycan_type='N')
    structures = glycan.predict_structures()
    
    # Visualize with color-coded levels (by varying spacing)
    print("\nVisualizing with progressive spacing increase by level:")
    GlycanVisualizer.visualize_with_matplotlib(
        structures[0],
        title="Progressive Level Spacing (6623)",
        vertical_gap={
            0: 0.8,   # Level 0: Tight
            1: 1.0,   # Level 1: Slightly wider
            2: 1.3,   # Level 2: Mid
            3: 1.6,   # Level 3: Wider (branches start)
            4: 2.0,   # Level 4: Wide (branching HexNAc)
            5: 1.8,   # Level 5: Still wide
            6: 1.5    # Level 6: Terminal
        },
        horizontal_spacing={
            4: 3.0,   # Very wide at branching HexNAc level
            5: 2.5    # Wide at terminal level
        },
        output_path="viz_example5_level_definitions.png",
        show=False
    )
    print("   Saved: viz_example5_level_definitions.png")
    print("   (Notice how spacing increases as we move up the tree)")


def main():
    """
    Run all examples
    """
    print("\n" + "="*80)
    print("GLYCOFRAG VISUALIZATION CONTROL EXAMPLES")
    print("="*80)
    print("\nThis script demonstrates comprehensive control over glycan visualization.")
    print("All output images will be saved to the current directory.")
    print("\nPress Ctrl+C to cancel at any time...")
    
    try:
        example1_direct_visualizer()
        example2_glycan_analysis()
        example3_glycopeptide_analysis()
        example4_glycan_facade()
        example5_level_definitions()
        
        print("\n" + "="*80)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY")
        print("="*80)
        print(f"\n✓ Generated visualization examples for all 4 entry points:")
        print("  1. GlycanVisualizer (direct)")
        print("  2. GlycanAnalysis (released glycan)")
        print("  3. GlycoPeptideAnalysis (glycopeptide)")
        print("  4. Glycan facade (intermediate)")
        print("\nCheck the current directory for all PNG files (viz_example*.png)")
        
    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user.")
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
