"""
Batch processing utilities for glycofrag.

This module provides efficient batch processing of multiple glycopeptides,
enabling high-throughput analysis of glycoproteomics data.
"""

from typing import List, Dict, Tuple, Optional, Any
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

from glycofrag import Glycopeptide, Glycan
from glycofrag.io.tables import generate_fragment_table


def batch_process_glycopeptides(
    glycopeptide_list: List[Dict[str, Any]],
    charge_states: Tuple[int, ...] = (1, 2, 3),
    use_cam: bool = True,
    max_workers: Optional[int] = None,
    show_progress: bool = False
) -> pd.DataFrame:
    """
    Process multiple glycopeptides in batch and return combined fragment table.
    
    This function efficiently processes a list of glycopeptide specifications
    and generates a combined DataFrame containing all fragment ions.
    
    Args:
        glycopeptide_list: List of dicts with keys:
            - 'peptide_sequence': str (amino acid sequence)
            - 'glycan_code': str (glycan composition code)
            - 'glycosylation_site': int (1-indexed glycan attachment site)
            - 'glycan_type': str (optional, default 'N')
            - 'id': str (optional, unique identifier)
        charge_states: Charge states to consider for fragments
        use_cam: Whether to apply carbamidomethylation
        max_workers: Maximum parallel workers (None = auto)
        show_progress: Print progress messages
    
    Returns:
        DataFrame with all fragments from all glycopeptides, with added
        'GlycopeptideID' column to track source.
    
    Example:
        >>> glycopeptides = [
        ...     {
        ...         'peptide_sequence': 'EEQYNSTYR',
        ...         'glycan_code': '4501',
        ...         'glycosylation_site': 5,
        ...         'glycan_type': 'N',
        ...         'id': 'IgG_G0F'
        ...     },
        ...     {
        ...         'peptide_sequence': 'GTTPSPVPTR',
        ...         'glycan_code': '1100',
        ...         'glycosylation_site': 2,
        ...         'glycan_type': 'O',
        ...         'id': 'Mucin_T'
        ...     }
        ... ]
        >>> df = batch_process_glycopeptides(glycopeptides)
        >>> print(f"Generated {len(df)} total fragments")
    """
    if not glycopeptide_list:
        return pd.DataFrame()
    
    start_time = time.time()
    all_dataframes = []
    
    if show_progress:
        print(f"Processing {len(glycopeptide_list)} glycopeptides...")
    
    # Process each glycopeptide
    for idx, gp_spec in enumerate(glycopeptide_list):
        gp_id = gp_spec.get('id', f"GP_{idx+1}")
        
        try:
            # Create glycopeptide
            glycopeptide = Glycopeptide(
                peptide_sequence=gp_spec['peptide_sequence'],
                glycan_code=gp_spec['glycan_code'],
                glycosylation_site=gp_spec['glycosylation_site'],
                glycan_type=gp_spec.get('glycan_type', 'N'),
                use_cam=use_cam
            )
            
            # Check if glycan structures were predicted
            # Glycopeptide stores structures in glycan_structures
            if len(glycopeptide.glycan_structures) == 0:
                if show_progress:
                    print(f"  Warning: No structures predicted for {gp_id}")
                continue
            
            # Generate glycan-only fragments for the table
            # Use the glycan directly to get pure glycan fragments
            glycan_structure = glycopeptide.glycan_structures[0]
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                structure=glycan_structure,
                fragment_types=['BY', 'A'],
                charges=[1, 2, 3]
            )
            
            # Extract glycan fragment compositions (not peptide fragments)
            glycan_fragments_dict = {
                'b_ions': glycan_fragments.get('b_ions', []),
                'y_ions': glycan_fragments.get('y_ions', []),
                'by_ions': glycan_fragments.get('by_ions', []),
                'yy_ions': glycan_fragments.get('yy_ions', []),
                'a_ions': glycan_fragments.get('a_ions', [])
            }
            
            # Generate fragment table
            df = generate_fragment_table(
                fragments=glycan_fragments_dict,
                glycan_code=gp_spec['glycan_code'],
                peptide=gp_spec['peptide_sequence'],
                glycan_type=gp_spec.get('glycan_type', 'N'),
                modification_type=6,  # Glycopeptide modification type
                charge_states=charge_states,
                use_cam=use_cam
            )
            
            # Add glycopeptide ID
            df['GlycopeptideID'] = gp_id
            
            all_dataframes.append(df)
            
            if show_progress:
                print(f"  ✓ {gp_id}: {len(df)} fragments")
        
        except Exception as e:
            if show_progress:
                print(f"  ✗ Error processing {gp_id}: {e}")
    
    # Combine all dataframes
    if not all_dataframes:
        return pd.DataFrame()
    
    result_df = pd.concat(all_dataframes, ignore_index=True)
    
    if show_progress:
        elapsed = time.time() - start_time
        print(f"\nCompleted: {len(result_df)} total fragments from "
              f"{len(all_dataframes)}/{len(glycopeptide_list)} glycopeptides "
              f"in {elapsed:.2f}s")
    
    return result_df


def batch_process_with_multiple_structures(
    glycopeptide_list: List[Dict[str, Any]],
    max_structures_per_glycan: int = 3,
    charge_states: Tuple[int, ...] = (1, 2, 3),
    use_cam: bool = True,
    show_progress: bool = False
) -> pd.DataFrame:
    """
    Process glycopeptides with multiple glycan structures per composition.
    
    For each glycopeptide, generates fragments for up to N glycan structures,
    useful for exploring structural heterogeneity.
    
    Args:
        glycopeptide_list: List of glycopeptide specifications
        max_structures_per_glycan: Maximum structures to process per glycan
        charge_states: Charge states for fragments
        use_cam: Apply carbamidomethylation
        show_progress: Print progress messages
    
    Returns:
        DataFrame with fragments from all structures, with 'StructureID' column
    """
    if not glycopeptide_list:
        return pd.DataFrame()
    
    start_time = time.time()
    all_dataframes = []
    total_structures = 0
    
    if show_progress:
        print(f"Processing {len(glycopeptide_list)} glycopeptides "
              f"(up to {max_structures_per_glycan} structures each)...")
    
    for idx, gp_spec in enumerate(glycopeptide_list):
        gp_id = gp_spec.get('id', f"GP_{idx+1}")
        
        try:
            # Create glycopeptide
            glycopeptide = Glycopeptide(
                peptide_sequence=gp_spec['peptide_sequence'],
                glycan_code=gp_spec['glycan_code'],
                glycosylation_site=gp_spec['glycosylation_site'],
                glycan_type=gp_spec.get('glycan_type', 'N'),
                use_cam=use_cam
            )
            
            # Get predicted structures (stored in glycan_structures)
            structures = glycopeptide.glycan_structures[:max_structures_per_glycan]
            
            if not structures:
                if show_progress:
                    print(f"  Warning: No structures predicted for {gp_id}")
                continue
            
            # Process each structure
            for struct_idx in range(len(structures)):
                # Generate glycan-only fragments for the table
                # Use the glycan directly to get pure glycan fragments
                glycan_structure = structures[struct_idx]
                glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                    structure=glycan_structure,
                    fragment_types=['BY', 'A'],
                    charges=[1, 2, 3]
                )
                
                # Extract glycan fragment compositions (not peptide fragments)
                glycan_fragments_dict = {
                    'b_ions': glycan_fragments.get('b_ions', []),
                    'y_ions': glycan_fragments.get('y_ions', []),
                    'by_ions': glycan_fragments.get('by_ions', []),
                    'yy_ions': glycan_fragments.get('yy_ions', []),
                    'a_ions': glycan_fragments.get('a_ions', [])
                }
                
                df = generate_fragment_table(
                    fragments=glycan_fragments_dict,
                    glycan_code=gp_spec['glycan_code'],
                    peptide=gp_spec['peptide_sequence'],
                    glycan_type=gp_spec.get('glycan_type', 'N'),
                    modification_type=6,
                    charge_states=charge_states,
                    use_cam=use_cam
                )
                
                # Add identifiers
                df['GlycopeptideID'] = gp_id
                df['StructureID'] = f"{gp_id}_S{struct_idx+1}"
                
                all_dataframes.append(df)
                total_structures += 1
            
            if show_progress:
                print(f"  ✓ {gp_id}: {len(structures)} structures")
        
        except Exception as e:
            if show_progress:
                print(f"  ✗ Error processing {gp_id}: {e}")
    
    # Combine all dataframes
    if not all_dataframes:
        return pd.DataFrame()
    
    result_df = pd.concat(all_dataframes, ignore_index=True)
    
    if show_progress:
        elapsed = time.time() - start_time
        print(f"\nCompleted: {len(result_df)} total fragments from "
              f"{total_structures} structures in {elapsed:.2f}s")
    
    return result_df


def generate_glycan_library(
    peptide_sequence: str,
    glycan_codes: List[str],
    glycosylation_site: int,
    glycan_type: str = 'N',
    use_cam: bool = True,
    charge_states: Tuple[int, ...] = (1, 2, 3),
    show_progress: bool = False
) -> pd.DataFrame:
    """
    Generate fragment library for a single peptide with multiple glycans.
    
    Useful for creating spectral libraries for peptides with known
    glycan microheterogeneity.
    
    Args:
        peptide_sequence: Amino acid sequence
        glycan_codes: List of glycan composition codes
        glycosylation_site: Glycan attachment site (1-indexed)
        glycan_type: 'N' or 'O'
        use_cam: Apply carbamidomethylation
        charge_states: Charge states for fragments
        show_progress: Print progress messages
    
    Returns:
        DataFrame with all fragments across glycan variants
    """
    glycopeptide_list = [
        {
            'peptide_sequence': peptide_sequence,
            'glycan_code': code,
            'glycosylation_site': glycosylation_site,
            'glycan_type': glycan_type,
            'id': f"{peptide_sequence}_{code}"
        }
        for code in glycan_codes
    ]
    
    return batch_process_glycopeptides(
        glycopeptide_list=glycopeptide_list,
        charge_states=charge_states,
        use_cam=use_cam,
        show_progress=show_progress
    )


def summarize_batch_results(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Generate summary statistics for batch processing results.
    
    Args:
        df: DataFrame from batch_process_glycopeptides
    
    Returns:
        Dictionary with summary statistics
    """
    if df.empty:
        return {
            'total_fragments': 0,
            'total_glycopeptides': 0,
            'fragment_types': []
        }
    
    summary = {
        'total_fragments': len(df),
        'total_glycopeptides': df['GlycopeptideID'].nunique() if 'GlycopeptideID' in df.columns else 0,
        'fragment_types': df['FragmentType'].value_counts().to_dict() if 'FragmentType' in df.columns else {},
        'charge_states': df['Ions'].value_counts().to_dict() if 'Ions' in df.columns else {},
        'mass_range': {
            'min': float(df['Fragment_mz'].min()) if 'Fragment_mz' in df.columns else 0,
            'max': float(df['Fragment_mz'].max()) if 'Fragment_mz' in df.columns else 0
        }
    }
    
    if 'GlycanType' in df.columns:
        summary['glycan_types'] = df['GlycanType'].value_counts().to_dict()
    
    return summary
