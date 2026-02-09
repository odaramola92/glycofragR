"""
Glycopeptide Structure Prediction and Fragmentation Example

This example demonstrates:
1. Creating a glycopeptide with specific sequence and glycan composition
2. Generating fragments across multiple fragment types
3. Displaying fragments with neutral masses and m/z for charge states +1, +2, +3
4. Saving fragments to Excel for analysis

Configuration:
- PEPTIDE_SEQUENCE: Target peptide sequence
- GLYCAN_CODE: 4-digit composition code (e.g., '4501' for Hex5HexNAc4NeuAc1)
- GLYCOSYLATION_SITE: Position of N/O-glycosylation (1-indexed)
- GLYCAN_TYPE: 'N' for N-glycans or 'O' for O-glycans
"""

import sys
import io
import os
from typing import Dict, List, Any

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    os.system('chcp 65001 > nul')
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)

from glycofrag import GlycoPeptideAnalysis

def main():
    # ========== CONFIGURATION ==========
    
    PEPTIDE_SEQUENCE = "LCPDCPLLAPLNDSR"
    GLYCAN_CODE = "4501"  # Hex5HexNAc4NeuAc1
    GLYCOSYLATION_SITE = 12  # N at position 12
    GLYCAN_TYPE = 'N'
    
    OUTPUT_FILE = f"glycopeptide_fragments_{PEPTIDE_SEQUENCE}_{''.join(filter(str.isalnum, GLYCAN_CODE))}.xlsx"
    
    # =================================================
    
    print("="*80, flush=True)
    print("GLYCOPEPTIDE FRAGMENTATION EXAMPLE", flush=True)
    print("="*80, flush=True)
    print(f"Peptide: {PEPTIDE_SEQUENCE}", flush=True)
    print(f"Glycan: {GLYCAN_CODE} ({GLYCAN_TYPE}-glycan)", flush=True)
    print(f"Glycosylation Site: Position {GLYCOSYLATION_SITE} ({PEPTIDE_SEQUENCE[GLYCOSYLATION_SITE-1]})", flush=True)
    print("="*80, flush=True)
    
    # Create glycopeptide analysis
    try:
        analysis = GlycoPeptideAnalysis(
            PEPTIDE_SEQUENCE,
            GLYCAN_CODE,
            GLYCOSYLATION_SITE,
            GLYCAN_TYPE,
            use_cam=True,  # Carbamidomethylation on Cys
            max_structures=100,
            mod_string="",  # N-terminal acetylation
            peptide_fragment_types=['by'],
            glycan_fragment_types=['BY']
        )
        print(f"\n[OK] Glycopeptide analysis created successfully", flush=True)
    except Exception as e:
        print(f"[ERROR] Failed to create glycopeptide analysis: {e}", flush=True)
        return
    
    # Generate fragments using updated API
    # structure_index parameter: None (default) = all structures, or specify 1-indexed structure number
    # peptide_fragment_types: ['by'] for b/y ions, ['cz'] for c/z ions, or ['by', 'cz'] for both
    # glycan_fragment_types: ['BY'] for B/Y ions, ['CZ'] for C/Z ions, or ['BY', 'CZ'] for both
    theoretical_df = analysis.theoretical_df
    total_frags = len(theoretical_df.index) if theoretical_df is not None else 0
    print(f"[OK] Generated {total_frags} theoretical fragments from all structures", flush=True)
    print(f"     Parameters: structure_index=None, peptide_fragment_types=['by'], glycan_fragment_types=['BY']", flush=True)
    print(f"     Note: structure_index=None aggregates unique fragments across all predicted glycan structures", flush=True)
    
    # Calculate masses
    peptide_mass = analysis.gp.peptide.mass
    glycan_mass = analysis.gp.mass_calculator.calculate_glycan_mass(GLYCAN_CODE)
    glycopeptide_mass = peptide_mass + glycan_mass
    
    print(f"\n[OK] Calculated masses:", flush=True)
    print(f"     Peptide mass: {peptide_mass:.4f} Da", flush=True)
    print(f"     Glycan mass: {glycan_mass:.4f} Da", flush=True)
    print(f"     Glycopeptide mass: {glycopeptide_mass:.4f} Da", flush=True)
    
    # Display fragments organized by type
    print("\n" + "="*80, flush=True)
    print("FRAGMENTS WITH MASSES (Neutral mass and m/z for charge states +1, +2, +3)", flush=True)
    print("="*80, flush=True)
    print("(Unique fragments aggregated from structure prediction)", flush=True)
    print("="*80, flush=True)
    
    # Organize by display groups from table
    display_groups = {
        'GLYCOPEPTIDE': ['intact'],
        'Y1': ['y1'],
        'Y0': ['y0'],
        'Y1F': ['y_glycan'],
        'B IONS': ['B'],
        'Y IONS': ['Y'],
        'BY IONS': ['BY'],
        'YY IONS': ['YY'],
        'C IONS': ['C'],
        'Z IONS': ['Z'],
        'CZ IONS': ['CZ'],
        'ZZ IONS': ['ZZ'],
        'PEPTIDE B IONS': ['b'],
        'PEPTIDE Y IONS': ['y'],
        'PEPTIDE C IONS': ['c'],
        'PEPTIDE Z IONS': ['z'],
        'A IONS': ['A'],
        'CUSTOM IONS': ['custom']
    }

    if theoretical_df is not None and not theoretical_df.empty:
        for group_name, frag_types in display_groups.items():
            group_df = theoretical_df[theoretical_df['Fragment Type'].isin(frag_types)]
            if group_df.empty:
                continue

            print(f"\n{group_name}:", flush=True)
            print("-" * 80, flush=True)

            group_sorted = group_df.sort_values(by='Mass(Da)')
            for _, row in group_sorted.iterrows():
                name = str(row['Fragment'])
                mass = float(row['Mass(Da)'])
                mz_1 = row.get('m/z(z=1)', '')
                mz_2 = row.get('m/z(z=2)', '')
                mz_3 = row.get('m/z(z=3)', '')

                print(f"  {name:35s} {mass:10.4f} Da", flush=True)
                print(f"    +1: {mz_1:10}  +2: {mz_2:10}  +3: {mz_3:10}", flush=True)
    
    print("\n" + "="*80, flush=True)
    
    # Save to Excel
    print(f"\n[Step] Saving fragments to Excel...", flush=True)
    analysis.export_to_excel(OUTPUT_FILE)


if __name__ == "__main__":
    main()
