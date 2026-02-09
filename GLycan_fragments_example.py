"""
Glycan Structure Prediction and Fragmentation Example (Simplified with GlycanAnalysis)

This example demonstrates the high-level GlycanAnalysis API:
1. Structure prediction with isomer_sensitive mode
2. Automatic fragment generation and table creation
3. Visualization with/without node numbers
4. Excel export with all data
5. O-glycan core selection and visualization

Configuration Options:
- GLYCAN_TYPE: 'N' for N-glycans or 'O' for O-glycans
- GLYCAN_CODE: 4-digit composition code (e.g., '1101' for HexNAc(1)Hex(1)Fuc(0)NeuAc(1))
- MODIFICATION_TYPE: Reducing end modification
  * 'free': Free reducing end
  * 'reduced': Reduced/alditol end
  * 'permethylated_free' or 'pm_reduced': Permethylated free
  * 'permethylated_reduced' or 'pm_reduced': Permethylated reduced
  * '2ab': 2-AB labeled
  * '2ab_permethylated': 2-AB + permethylated
  * 'custom': Custom (requires custom_reducing_end_mass parameter)
  
- PREFERRED_CORE: For O-glycans only - select specific core (0-8) for visualization
  * None (default): Uses default core for each group
  * 1-8: Select specific core (must be predicted for composition)
  
- ISOMER_SENSITIVE: Controls how mirror images are handled
  * False (default): Deduplicates topologically identical mirror images (fewer structures)
  * True: Treats mirror images as distinct structures (exhaustive enumeration)
  
- FRAGMENT_TYPES: Fragment series to generate
  * ['BY']: B and Y ions (default)
  * ['CZ']: C and Z ions
  * ['BY', 'CZ']: Both series
  
- CHARGES: Charge states for m/z calculation
  * [1, 2, 3]: Default charge states
  * [1, 2, 3, 4, 5]: Extended range
  
- VISUALIZE_STRUCTURE: Which structures to visualize
  * 'all': All predicted structures
  * 'best': First structure (or best after MS/MS matching)
  * int (1, 2, 3...): Specific structure number
  * [1, 2, 3]: List of structures
    * None: No visualization
  
- CUSTOM_VERTICAL_GAP: Advanced control over vertical spacing between levels (for manual visualization)
    * None (default): Use built-in spacing (1.2)
    * Single float (e.g., 1.8): Apply uniformly to all levels
    * Dict (e.g., {3: 2.0, 4: 2.5}): Specify per level (0=root, 1=level 1, etc.)

- CUSTOM_HORIZONTAL_SPACING: Advanced control over horizontal spacing at each level (for manual visualization)
    * None (default): Use built-in spacing (1.5)
    * Single float (e.g., 2.5): Apply uniformly
    * Dict (e.g., {4: 3.0, 5: 3.5}): Widen specific levels (useful for branching areas)

- CUSTOM_NODE_SIZE: Advanced control over monosaccharide marker sizes (for manual visualization)
    * None (default): Use built-in sizes (200 for circles, 120 for sialic acids)
    * Single float (e.g., 300): Apply uniformly to all nodes
    * Dict (e.g., {'HexNAc': 350, 'NeuAc': 180}): Specify per monosaccharide type

NOTES FOR O-GLYCANS:
- GalNAc (yellow boxes) vs GlcNAc (blue boxes) properly distinguished
- Cores within same group generate IDENTICAL fragments
- Core selection only affects visualization, not fragmentation or masses

VISUALIZATION CONTROL EXAMPLES:
- Publication quality: CUSTOM_VERTICAL_GAP={3: 2.0, 4: 2.5}, CUSTOM_HORIZONTAL_SPACING={4: 3.5}
- Compact layout: CUSTOM_VERTICAL_GAP=0.8, CUSTOM_HORIZONTAL_SPACING=1.2, CUSTOM_NODE_SIZE=180
- Emphasize HexNAc: CUSTOM_NODE_SIZE={'HexNAc': 400, 'Man': 250, 'Gal': 250, 'NeuAc': 180}

NEW SIMPLIFIED API:
- Only need to import GlycanAnalysis
- Fragment generation, table creation, and export handled automatically
- Access data via analysis.theoretical_df, analysis.clean_theoretical_df, analysis.summary_df
"""

import sys
import io
import os

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    # Change console code page to UTF-8
    os.system('chcp 65001 > nul')
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)


from glycofrag import GlycanAnalysis

def main():
    # ========== CONFIGURATION ==========
    
    # Glycan selection
    GLYCAN_TYPE = 'O'  # 'N' for N-glycans, 'O' for O-glycans
    GLYCAN_CODE = '1103'  # 5-digit composition code
    
    # Reducing end modification
    MODIFICATION_TYPE = 'permethylated_reduced'  # or 'pm_reduced', 'free', 'reduced', '2ab', etc.
    
    # Structure prediction settings
    ISOMER_SENSITIVE = False  # Set to True to treat mirror images as distinct structures
    MAX_STRUCTURES = 100      # Maximum number of structures to predict
    
    # O-glycan specific settings (only used when GLYCAN_TYPE='O')
    PREFERRED_CORE = None  # Set to 1-8 to select specific core for visualization
                          # Examples: PREFERRED_CORE=5 to force Core 5 interpretation
                          # None (default): use default core for each group
    
    # Fragment generation settings
    FRAGMENT_TYPES = ['BY']  # ['BY'], ['CZ'], or ['BY', 'CZ']
    CHARGES = [1, 2, 3]     # Charge states for m/z calculation
    
    # Visualization settings
    VISUALIZE_STRUCTURE = 'all' # 'all', 'best', int (1,2,3...), [1,2,3], or None
    
    # Excel export
    OUTPUT_FILE = f"glycan_analysis_{GLYCAN_CODE}_{GLYCAN_TYPE}.xlsx"
    
    # =================================================
    
    print("="*80)
    print("GLYCAN FRAGMENTATION ANALYSIS")
    print("="*80)
    print(f"Glycan Type: {GLYCAN_TYPE}, Code: {GLYCAN_CODE}")
    print(f"Modification: {MODIFICATION_TYPE}")
    if GLYCAN_TYPE == 'O' and PREFERRED_CORE is not None:
        print(f"Preferred core: Core {PREFERRED_CORE}")
    print(f"Configuration: isomer_sensitive={ISOMER_SENSITIVE}, fragment_types={FRAGMENT_TYPES}")
    print("="*80)
    
    # Create GlycanAnalysis - handles everything automatically!
    try:
        analysis = GlycanAnalysis(
            glycan_code=GLYCAN_CODE,
            glycan_type=GLYCAN_TYPE,
            modification_type=MODIFICATION_TYPE,
            max_structures=MAX_STRUCTURES,
            isomer_sensitive=ISOMER_SENSITIVE,
            fragment_types=FRAGMENT_TYPES,
            charges=CHARGES,
            visualize_structure=VISUALIZE_STRUCTURE,
            preferred_core=PREFERRED_CORE,  # O-glycan only
        )
        print(f"\n[OK] GlycanAnalysis created successfully")
        print(f"     {len(analysis.structures)} structures predicted")
        print(f"     {len(analysis.theoretical_df)} theoretical fragments generated")
        print(f"     {len(analysis.clean_theoretical_df)} unique fragments (deduplicated)")
    except Exception as e:
        print(f"[ERROR] Failed to create analysis: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Display summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(analysis.summary_df.to_string(index=False))
    
    # For O-glycans, show core information
    if GLYCAN_TYPE == 'O':
        cores_info = analysis.glycan.get_predicted_cores()
        if cores_info:
            print(f"\n" + "="*80)
            print("O-GLYCAN CORE INFORMATION")
            print("="*80)
            print("(Cores within the same group generate identical fragments and masses)")
            for i, core in enumerate(cores_info, 1):
                possible = core['possible_cores']
                default = core['default_core']
                selected = core['selected_core']
                name = core['core_name']
                
                core_str = f"Cores {possible}" if len(possible) > 1 else f"Core {possible[0]}"
                print(f"  {i}. {name}")
                print(f"     {core_str}, Default: {default}, Selected: {selected}")
    
    # Show structures with masses and classification
    print("\n" + "="*80)
    print("PREDICTED STRUCTURES")
    print("="*80)
    
    for i, struct in enumerate(analysis.structures):
        composition = analysis.glycan.count_residues(struct)
        mass = struct.graph['mass']
        breakdown = struct.graph['mass_breakdown']
        classification = analysis.glycan.classify_structure(struct)
        
        # For O-glycans, show core info
        core_info_str = ""
        if GLYCAN_TYPE == 'O':
            possible_cores = struct.graph.get('possible_cores', [])
            selected_core = struct.graph.get('selected_core', 'N/A')
            if possible_cores:
                core_info_str = f" [Cores: {possible_cores}, Selected: {selected_core}]"
        
        print(f"\nStructure {i+1}{core_info_str}:")
        print(f"Type: {classification}")
        print(f"Mass: {mass:.4f} Da")
       
        print(f"Composition: {composition}")
        print("Tree:")
        tree_str = analysis.glycan.format_structure_tree(struct)
        for line in tree_str.split('\n'):
            print(f"  {line}")
    
    # Display fragment summary
    print("\n" + "="*80)
    print("FRAGMENT SUMMARY")
    print("="*80)
    print(f"Total theoretical fragments: {len(analysis.theoretical_df)}")
    print(f"Unique fragments (deduplicated): {len(analysis.clean_theoretical_df)}")
    
    if not analysis.clean_theoretical_df.empty:
        # Show fragment type distribution
        if 'Fragment Type' in analysis.clean_theoretical_df.columns:
            frag_counts = analysis.clean_theoretical_df['Fragment Type'].value_counts()
            print("\nFragments by type:")
            for frag_type, count in frag_counts.items():
                print(f"  {frag_type}: {count}")
    
    # ========== OPTIONAL: Custom visualization with advanced parameters ==========
    # Uncomment the code below to use custom visualization parameters
    # Note: You can now pass these parameters directly to analysis.visualize_structure()
    #       without needing to import GlycanVisualizer!
    
    # Example 1: Custom visualization with specific parameters (using GlycanAnalysis)
    # try:
    #     print("\n" + "="*80)
    #     print("CUSTOM VISUALIZATION - Example 1 (Tight/Compact Layout)")
    #     print("="*80)
    #     print("This compresses the structure vertically and horizontally")
    #     print("Notice how Nodes 4-5 (branch mannoses) are VERY CLOSE together")
    #     analysis.visualize_structure(
    #         structure_number=1,
    #         output_path='glycan_compact_viz.png',
    #         vertical_gap=0.6,                   # Very small vertical gaps everywhere
    #         horizontal_spacing=0.8,             # Very small horizontal gaps everywhere
    #         node_size=150                       # Small nodes
    #     )
    #     print("[OK] Compact visualization saved as glycan_compact_viz.png")
    # except Exception as e:
    #     print(f"[INFO] Compact visualization example not available: {e}")
    
    
    # Example 2: Per-level/per-type custom visualization (using GlycanAnalysis)
    # try:
    #     print("\n" + "="*80)
    #     print("CUSTOM VISUALIZATION - TEST 1 (Small vertical gap at level 2)")
    #     print("="*80)
    #     print("Parameters: vertical_gap={2: 1.0}")
    #     print("This will make Node 3 very CLOSE to Node 2")
        
    #     analysis.visualize_structure(
    #         structure_number=1,
    #         output_path='glycan_test_small_gap.png',
    #         vertical_gap={2: 1.0}
    #     )
    #     print("[OK] Test 1 saved as glycan_test_small_gap.png")
    # except Exception as e:
    #     print(f"[ERROR] Test 1 failed: {e}")
    #     import traceback
    #     traceback.print_exc()
    
    # try:
    #     print("\n" + "="*80)
    #     print("CUSTOM VISUALIZATION - TEST 2 (Large vertical gap at level 2)")
    #     print("="*80)
    #     print("Parameters: vertical_gap={2: 6.0}")
    #     print("This will make Node 3 very FAR from Node 2")
        
    #     analysis.visualize_structure(
    #         structure_number=1,
    #         output_path='glycan_test_large_gap.png',
    #         vertical_gap={2: 6.0}
    #     )
    #     print("[OK] Test 2 saved as glycan_test_large_gap.png")
    #     print("\n" + "="*80)
    #     print("COMPARE: glycan_test_small_gap.png vs glycan_test_large_gap.png")
    #     print("You should see a BIG difference in spacing between Node 2 and Node 3!")
    #     print("="*80)
    # except Exception as e:
    #     print(f"[ERROR] Test 2 failed: {e}")
    #     import traceback
    #     traceback.print_exc()
    
    # Original example
    try:
        print("\n" + "="*80)
        print("CUSTOM VISUALIZATION - Example (Your Settings)")
        print("="*80)
        print("Depth Level Guide:")
        print("  Depth 0: Node 1 (HexNAc reducing end)")
        print("  Depth 1: Node 2 (HexNAc core)")
        print("  Depth 2: Node 3 (Mannose core)")
        print("  Depth 3: Nodes 4-5 (Branch Mannose - SIBLINGS, spread by horizontal_spacing)")
        print("  Depth 4: Nodes 6-9 (Branch HexNAc)")
        print("  Depth 5: Nodes 10-12 (Gal - terminal)")
        print("  Depth 6: Nodes 13-14 (NeuAc - terminal sialic acids)")
        print("\nParameters:")
        print("  vertical_gap={2: 3.0, 3: 2.0}")
        print("    → Node 3 is 3.0 units below Node 2")
        print("    → Nodes 4-5 are 2.0 units below Node 3")
        print("  horizontal_spacing={3: 5.0, 4: 2.5}")
        print("    → Nodes 4-5 (siblings) are SPREAD FAR APART (5.0 units)")
        print("    → Nodes 6-9 stay closer (2.5 units)")
        print("  node_size for emphasis: HexNAc=400, NeuAc=150")
        print("\n" + "="*80)
        
        analysis.visualize_structure(
            structure_number=1,
            output_path='glycan_advanced_viz.png',
            vertical_gap={2: 1.0, 3: 1.0},      # Spread vertically
            horizontal_spacing={2: 4.0, 3: 2}   # Nodes 4-5 FAR APART, nodes 6-9 closer
            #node_size={'HexNAc': 400, 'NeuAc': 150}  # Large HexNAc, small NeuAc
        )
        print("[OK] Spacious visualization saved as glycan_advanced_viz.png")
    except Exception as e:
        print(f"[INFO] Spacious visualization example not available: {e}")
    
    # # Example 3: Extreme example to show the difference
    # try:
    #     print("\n" + "="*80)
    #     print("CUSTOM VISUALIZATION - Example 3 (Extreme: Very Wide Branches)")
    #     print("="*80)
    #     print("This MASSIVELY spreads Nodes 4-5 apart to show the effect")
    #     analysis.visualize_structure(
    #         structure_number=1,
    #         output_path='glycan_extreme_viz.png',
    #         horizontal_spacing={3: 10.0},  # Nodes 4-5: EXTREMELY FAR apart
    #         vertical_gap={3: 0.8},            # But compressed vertically
    #         node_size={'HexNAc': 300}
    #     )
    #     print("[OK] Extreme visualization (very wide branches) saved as glycan_extreme_viz.png")
    # except Exception as e:
    #     print(f"[INFO] Extreme visualization example not available: {e}")
    
    # # Export to Excel
    # print("\n" + "="*80)
    # print("EXPORTING RESULTS")
    # print("="*80)
    # try:
    #     analysis.export_to_excel(OUTPUT_FILE, include_prediction=False)
    #     print(f"[OK] All analysis exported successfully!")
    #     print(f"     Sheets: Summary, Theoretical Fragments, Clean Theoretical")
    # except Exception as e:
    #     print(f"[ERROR] Failed to export: {e}")
    #     import traceback
    #     traceback.print_exc()
    
    # Final summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"Summary: {len(analysis.structures)} structures, "
          f"{len(analysis.clean_theoretical_df)} unique fragments")
    if VISUALIZE_STRUCTURE:
        print(f"Visualizations saved as: glycan_{GLYCAN_CODE}_structure_*.png")
    print(f"Excel file: {OUTPUT_FILE}")
    print("="*80)

if __name__ == "__main__":
    main()
