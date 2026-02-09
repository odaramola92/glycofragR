"""
Real Glycopeptide Analysis from Thermo RAW Files

Workflow:
1. Load RAW file and search for precursor at target RT ± RT tolerance
2. Collect all MS/MS spectra within RT range around target
3. Generate theoretical fragments using glycofrag
4. Match experimental peaks with theoretical fragments (±PPM tolerance)
5. Export matched fragments to Excel with PPM deviations and intensities

Example:
    python Real_glycopeptide_Test.py
    Then enter:
        - Raw file path: Fetuin_PRM_1_6mz_1.raw
        - Target m/z: 1218.8447
        - Charge: 3
        - Target RT (min): 62.1
        - RT tolerance (min): 1.0
        - Fragment m/z tolerance (ppm): 10.0
        - Peptide sequence: LCPDCPLLAPLNDSR
        - Glycan code: 4501 (for N-glycan example)
"""

import os
import fisher_py
from glycofrag import GlycoPeptideAnalysis
import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional, Any
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================================
# RAW File I/O Functions
# ============================================================================

def load_msms_spectrum_from_raw(raw_file: str, scan_number: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load MS/MS spectrum from Thermo RAW file using fisher_py.
    
    Args:
        raw_file: Path to .raw file
        scan_number: Scan number to load
    
    Returns:
        Tuple of (mz_array, intensity_array)
    """
    try:
        raw_reader = fisher_py.RawFile(raw_file)
        
        # fisher_py.get_scan_from_scan_number returns: (mz_array, intensity_array, charge_array, ?)
        mz_array, intensity_array, charge_array, _ = raw_reader.get_scan_from_scan_number(scan_number)
        
        return mz_array, intensity_array
        
    except Exception as e:
        raise ValueError(f"Failed to load scan {scan_number} from {raw_file}: {e}")

def get_scan_count(raw_file: str) -> int:
    """Get total number of scans in RAW file."""
    try:
        raw_reader = fisher_py.RawFile(raw_file)
        return raw_reader.number_of_scans
    except Exception as e:
        raise ValueError(f"Failed to read scan count from {raw_file}: {e}")

def get_scan_retention_time(raw_file: str, scan_number: int) -> float:
    """
    Get retention time (in minutes) for a specific scan.
    
    Returns:
        Retention time in minutes
    """
    try:
        raw_reader = fisher_py.RawFile(raw_file)
        rt = raw_reader.get_retention_time_from_scan_number(scan_number)
        return rt
    except Exception as e:
        raise ValueError(f"Failed to get RT for scan {scan_number}: {e}")

def get_scan_precursor_mz(raw_file: str, scan_number: int) -> Optional[float]:
    """
    Get precursor m/z for a specific scan (MS/MS).
    
    For PRM experiments, extracts from scan event string.
    
    Returns:
        Precursor m/z or None if not available
    """
    try:
        raw_reader = fisher_py.RawFile(raw_file)
        # Try to get precursor from scan event string
        # Format example: "FTMS + p NSI SRM ms2 m/z=1218.84 X 3 [1050.00-2000.00]"
        scan_event = raw_reader.get_scan_event_str_from_scan_number(scan_number)
        
        # Parse precursor m/z from scan event
        if "m/z=" in scan_event:
            parts = scan_event.split("m/z=")
            if len(parts) > 1:
                mz_str = parts[1].split()[0]  # Get first token after m/z=
                return float(mz_str)
        
        return None
    except:
        return None

# ============================================================================
# Precursor Targeting Functions
# ============================================================================

def find_precursor_scans(
    raw_reader,
    target_mz: float,
    target_rt: float,
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance: float = 1.0,
    charge_hint: Optional[int] = None
) -> List[Tuple[int, float, float]]:
    """
    Find all MS/MS scans matching target m/z and RT.
    
    Args:
        raw_reader: Already-opened fisher_py.RawFile object
        target_mz: Target precursor m/z (full precision, e.g., 1218.8447)
        target_rt: Target retention time (minutes)
        mz_tolerance_ppm: m/z tolerance (ppm)
        rt_tolerance: RT tolerance (minutes)
        charge_hint: Expected charge state (for validation, optional)
    
    Returns:
        List of (scan_number, mz, rt) tuples matching criteria
    """
    matching_scans = []
    all_precursors_in_rt_window = []  # Track all precursors found in RT window
    
    try:
        total_scans = raw_reader.number_of_scans
        
        print(f"Scanning {total_scans} scans for target m/z={target_mz:.2f} (±{mz_tolerance_ppm} ppm) at RT={target_rt:.2f}±{rt_tolerance:.2f} min...")
        
        # Calculate ppm-based tolerance
        mz_tolerance_da = target_mz * mz_tolerance_ppm / 1e6
        
        # Progress reporting interval
        report_interval = max(1000, total_scans // 20)
        
        # Diagnostic counters
        ms1_count = 0
        ms2_count = 0
        ms2_in_rt_window = 0
        
        for scan_num in range(1, total_scans + 1):
            try:
                # Progress indicator
                if scan_num % report_interval == 0:
                    print(f"  Progress: {scan_num}/{total_scans} scans ({100*scan_num//total_scans}%) - MS1:{ms1_count}, MS2:{ms2_count}")
                
                # Get retention time directly from raw_reader (no file reopen)
                rt = raw_reader.get_retention_time_from_scan_number(scan_num)
                
                # Get scan event to determine MS level and precursor
                scan_event = raw_reader.get_scan_event_str_from_scan_number(scan_num)
                
                # Check MS level (MS1 scans contain "Full ms", MS2 contain "ms2" or "SRM")
                is_ms2 = ("ms2" in scan_event.lower() or "srm" in scan_event.lower() or "prm" in scan_event.lower())
                
                if is_ms2:
                    ms2_count += 1
                else:
                    ms1_count += 1
                    continue  # Skip MS1 scans
                
                # Quick RT filter first (most efficient)
                rt_diff = abs(rt - target_rt)
                if rt_diff > rt_tolerance:
                    continue
                
                ms2_in_rt_window += 1
                
                # Parse precursor m/z from scan event
                # Format: "FTMS + p NSI Full ms2 711.0330@hcd26.67 [146.3333-2195.0000]"
                # or:     "FTMS + p NSI SRM ms2 m/z=1218.84 X 3 [1050.00-2000.00]"
                precursor_mz = None
                
                if "m/z=" in scan_event:
                    # Format 1: m/z=XXX
                    parts = scan_event.split("m/z=")
                    if len(parts) > 1:
                        mz_str = parts[1].split()[0]
                        precursor_mz = float(mz_str)
                elif "ms2 " in scan_event.lower():
                    # Format 2: ms2 XXX@hcd or ms2 XXX [
                    parts = scan_event.lower().split("ms2 ")
                    if len(parts) > 1:
                        # Extract m/z (number before @ or [ or space)
                        mz_part = parts[1].strip()
                        # Get first token, remove @ and everything after
                        mz_str = mz_part.split('@')[0].split('[')[0].split()[0]
                        try:
                            precursor_mz = float(mz_str)
                        except:
                            pass
                
                if precursor_mz is None:
                    continue
                
                # Store all precursors in RT window for diagnostics
                all_precursors_in_rt_window.append((scan_num, precursor_mz, rt))
                
                # Check if m/z matches criteria
                mz_diff = abs(precursor_mz - target_mz)
                
                if mz_diff <= mz_tolerance_da:
                    matching_scans.append((scan_num, precursor_mz, rt))
                    print(f"  Found: Scan {scan_num}, m/z={precursor_mz:.4f}, RT={rt:.3f} min")
            
            except:
                continue
        
        if not matching_scans:
            print(f"\n{'='*80}")
            print(f"WARNING: No matching scans found!")
            print(f"{'='*80}")
            print(f"\nFile statistics:")
            print(f"  Total scans: {total_scans}")
            print(f"  MS1 scans: {ms1_count}")
            print(f"  MS2 scans: {ms2_count}")
            print(f"  MS2 scans in RT window {target_rt-rt_tolerance:.2f}-{target_rt+rt_tolerance:.2f} min: {ms2_in_rt_window}")
            
            # Report what precursors ARE present in the RT window
            if all_precursors_in_rt_window:
                print(f"\nFound {len(all_precursors_in_rt_window)} total precursors in RT window {target_rt-rt_tolerance:.2f}-{target_rt+rt_tolerance:.2f} min")
                print(f"\nShowing 10 unique precursor examples (clustered by 20 ppm):\n")
                
                # Cluster precursors by 20 ppm to find unique ones
                unique_precursors = []
                cluster_tolerance_ppm = 20.0
                
                # Sort by m/z
                sorted_precursors = sorted(all_precursors_in_rt_window, key=lambda x: x[1])
                
                for scan_num, mz, rt in sorted_precursors:
                    # Check if this m/z is unique (not within 20 ppm of existing)
                    is_unique = True
                    for _, existing_mz, _ in unique_precursors:
                        ppm_diff = abs(mz - existing_mz) / existing_mz * 1e6
                        if ppm_diff <= cluster_tolerance_ppm:
                            is_unique = False
                            break
                    
                    if is_unique:
                        unique_precursors.append((scan_num, mz, rt))
                        if len(unique_precursors) >= 10:
                            break
                
                # Display the unique precursors
                for i, (scan_num, mz, rt) in enumerate(unique_precursors, 1):
                    ppm_from_target = abs(mz - target_mz) / target_mz * 1e6
                    print(f"  {i:2d}. m/z={mz:10.4f}  RT={rt:6.2f} min  (Scan {scan_num:5d})  [d {ppm_from_target:6.1f} ppm from target]")
                
                print(f"\n{'='*80}")
                print(f"HINT: Your target m/z={target_mz:.4f} was not found in the RT window.")
                print(f"      Check if your target m/z or RT values are correct.")
                print(f"{'='*80}\n")
            else:
                print(f"\nNo precursors found in MS2 scans within RT window {target_rt-rt_tolerance:.2f}-{target_rt+rt_tolerance:.2f} min")
                
                if ms2_in_rt_window == 0:
                    print(f"  → No MS2 scans found in this RT window at all!")
                else:
                    print(f"  → Found {ms2_in_rt_window} MS2 scans in RT window, but none matched target m/z")
                
                print(f"\nSearching entire file for precursors to help diagnose the issue...")
                
                # Sample the entire file to find where precursors actually exist
                sample_precursors = []
                sample_scan_events = []  # Store example scan events for debugging
                sample_interval = max(100, total_scans // 100)  # Sample ~100 scans across file
                
                for scan_num in range(1, total_scans + 1, sample_interval):
                    try:
                        rt = raw_reader.get_retention_time_from_scan_number(scan_num)
                        scan_event = raw_reader.get_scan_event_str_from_scan_number(scan_num)
                        
                        # Check if MS2
                        is_ms2 = ("ms2" in scan_event.lower() or "srm" in scan_event.lower() or "prm" in scan_event.lower())
                        
                        # Store first 5 MS2 scan events for debugging
                        if is_ms2 and len(sample_scan_events) < 5:
                            sample_scan_events.append((scan_num, rt, scan_event))
                        
                        if is_ms2:
                            # Parse precursor m/z (support both formats)
                            precursor_mz = None
                            if "m/z=" in scan_event:
                                parts = scan_event.split("m/z=")
                                if len(parts) > 1:
                                    mz_str = parts[1].split()[0]
                                    precursor_mz = float(mz_str)
                            elif "ms2 " in scan_event.lower():
                                parts = scan_event.lower().split("ms2 ")
                                if len(parts) > 1:
                                    mz_part = parts[1].strip()
                                    mz_str = mz_part.split('@')[0].split('[')[0].split()[0]
                                    try:
                                        precursor_mz = float(mz_str)
                                    except:
                                        pass
                            
                            if precursor_mz is not None:
                                sample_precursors.append((scan_num, precursor_mz, rt))
                    except:
                        continue
                
                if sample_precursors:
                    # Cluster by 20 ppm and show 10 examples
                    unique_precursors = []
                    cluster_tolerance_ppm = 20.0
                    sorted_precursors = sorted(sample_precursors, key=lambda x: x[1])
                    
                    for scan_num, mz, rt in sorted_precursors:
                        is_unique = True
                        for _, existing_mz, _ in unique_precursors:
                            ppm_diff = abs(mz - existing_mz) / existing_mz * 1e6
                            if ppm_diff <= cluster_tolerance_ppm:
                                is_unique = False
                                break
                        
                        if is_unique:
                            unique_precursors.append((scan_num, mz, rt))
                            if len(unique_precursors) >= 10:
                                break
                    
                    print(f"\nFound {len(sample_precursors)} precursors in sampled file")
                    print(f"\nShowing 10 unique precursor examples from entire file (clustered by 20 ppm):\n")
                    
                    for i, (scan_num, mz, rt) in enumerate(unique_precursors, 1):
                        ppm_from_target = abs(mz - target_mz) / target_mz * 1e6
                        print(f"  {i:2d}. m/z={mz:10.4f}  RT={rt:6.2f} min  (Scan {scan_num:5d})  [Δ {ppm_from_target:6.1f} ppm from target]")
                    
                    # Find RT range of data
                    min_rt = min(p[2] for p in sample_precursors)
                    max_rt = max(p[2] for p in sample_precursors)
                    
                    print(f"\n{'='*80}")
                    print(f"HINT: Your target RT={target_rt:.2f} min may be outside the data range.")
                    print(f"      Observed RT range in file: {min_rt:.2f} - {max_rt:.2f} min")
                    print(f"      Your target RT window: {target_rt-rt_tolerance:.2f} - {target_rt+rt_tolerance:.2f} min")
                    print(f"{'='*80}\n")
                else:
                    print(f"\nNo MS2 precursors with 'm/z=' found in scan events.")
                    print(f"\nExample MS2 scan event strings from file:")
                    if sample_scan_events:
                        for scan_num, rt, event_str in sample_scan_events:
                            print(f"  Scan {scan_num:5d} (RT={rt:6.2f}): {event_str[:120]}")
                    else:
                        print(f"  → No MS2 scans found in entire file!")
                    print(f"\nHINT: Check if this is a PRM/SRM RAW file with MS2 scans.")
                    print(f"      Expected MS2 scan events to contain 'ms2' or 'SRM' and 'm/z=XXX'.\n")
        else:
            print(f"[OK] Search complete: Found {len(matching_scans)} matching scans")
        
        return matching_scans
        
    except Exception as e:
        raise ValueError(f"Error finding precursor scans: {e}")

# ============================================================================
# Fragment Matching Functions
# ============================================================================

def calculate_ppm_error(theoretical_mz: float, experimental_mz: float) -> float:
    """Calculate PPM error between theoretical and experimental m/z."""
    return abs(theoretical_mz - experimental_mz) / theoretical_mz * 1e6

def format_fragment_name(frag_dict: Dict[str, Any], fragment_type: str) -> str:
    """
    Create a GlypPRM-style fragment name from composition and type.
    
    Examples:
        {'NeuAc': 1, ...} -> "NeuAc1"
        {'HexNAc': 1, 'Hex': 1, ...} -> "HexNAc1-Hex1"
        {'HexNAc': 2, 'Hex': 3, ...} -> "HexNAc2-Hex3"
    """
    # Map of abbreviation keys to proper names (in order of priority)
    name_map = [
        ('NeuAc', 'NeuAc'),
        ('Fuc', 'Fuc'),
        ('HexNAc', 'HexNAc'),
        ('Hex', 'Hex'),
    ]
    
    parts = []
    for key, display_name in name_map:
        count = frag_dict.get(key, 0)
        if count > 0:
            parts.append(f"{display_name}{count}")
    
    # Create base name
    base_name = '-'.join(parts) if parts else "Unknown"
    
    # Add water loss notation for certain fragments (B ions with water loss)
    if frag_dict.get('_water_loss', False):
        base_name += "(-H2O)"
    
    return base_name

def add_fragment_names_to_theoretical(theoretical_fragments: Dict[str, List[Dict]]) -> None:
    """
    Add proper 'name' field to all fragments for display in Excel.
    Modifies fragments in-place.
    """
    for frag_type, frag_list in theoretical_fragments.items():
        if not isinstance(frag_list, list):
            continue
        
        for frag in frag_list:
            if not isinstance(frag, dict):
                continue
            
            # Determine the name based on fragment type
            if frag_type.startswith('oxonium'):
                # For oxonium ions: use composition-based name
                # composition has already been set to the proper format
                name = frag.get('composition', 'Unknown')
            else:
                # For peptide/glycopeptide fragments: use existing name or _custom_label
                name = frag.get('name', frag.get('_custom_label', 'Unknown'))
            
            frag['name'] = name

# ============================================================================
# Main Workflow
# ============================================================================

def analyze_glycopeptide_from_raw(
    raw_file: str,
    precursor_mz: float,
    charge: int,
    target_rt: float,
    rt_tolerance: float,
    peptide_sequence: str,
    glycan_code: str,
    glycosylation_site: int,
    glycan_type: str = 'N',
    ppm_tolerance: float = 20.0,
    output_file: Optional[str] = None,
    mod_string: Optional[str] = None,
    best_structure_report: bool = True,
    structure_output_dir: Optional[str] = None,
    peptide_fragment_types: Optional[List[str]] = None,
    glycan_fragment_types: Optional[List[str]] = None,
    figsize: tuple = (5, 4),
):
    """
    Complete workflow: Load RAW, find precursor, match fragments, export results.
    
    Args:
        raw_file: Path to Thermo RAW file
        precursor_mz: Target precursor m/z
        charge: Charge state
        target_rt: Target retention time (minutes)
        rt_tolerance: RT search window (minutes)
        peptide_sequence: Peptide sequence (e.g., "LCPDCPLLAPLNDSR")
        glycan_code: Glycan composition code (e.g., "4501")
        glycosylation_site: 1-indexed position where glycan attaches (e.g., 7 for 'N' in LCPDCPLLAPLNDSR)
        glycan_type: 'N' or 'O' for N/O-glycan
        ppm_tolerance: Fragment m/z tolerance (ppm)
        output_file: Output Excel file path (auto-generated if None)
        mod_string: Modification string (optional)
        best_structure_report: Whether to predict and visualize structures (default: True)
        structure_output_dir: Directory to save structure images (e.g., 'output/'). If None, displays interactively.
        figsize: Tuple of (width, height) in inches for structure figures (default: (5, 4))
    """
    
    print("=" * 80)
    print("GLYCOPEPTIDE MS/MS ANALYSIS FROM THERMO RAW FILE")
    print("=" * 80)
    
    # Validate input
    raw_path = Path(raw_file)
    if not raw_path.exists():
        raise FileNotFoundError(f"RAW file not found: {raw_file}")
    
    # Open RAW file ONCE for the entire analysis
    print(f"\nOpening RAW file: {raw_path.name}...")
    raw_reader = fisher_py.RawFile(raw_file)
    print(f"[OK] RAW file opened successfully")
    
    # Step 1: Find precursor scans
    print(f"\n[Step 1] Finding precursor scans...")
    print(f"  Searching for m/z={precursor_mz:.4f} (±10.0 ppm) at RT={target_rt:.2f}±{rt_tolerance:.2f} min")
    matching_scans = find_precursor_scans(
        raw_reader,
        precursor_mz,
        target_rt,
        mz_tolerance_ppm=10.0,
        rt_tolerance=rt_tolerance,
        charge_hint=charge
    )
    
    if not matching_scans:
        raise ValueError("No matching precursor scans found!")
    
    print(f"[OK] Found {len(matching_scans)} matching scans")
    
    # Step 2: Combine MS/MS spectra (using cached raw_reader)
    print(f"\n[Step 2] Loading MS/MS spectra from {len(matching_scans)} scans...")
    all_mz = []
    all_intensity = []
    all_scan_rt = []  # Track RT for each fragment
    
    progress_interval = max(50, len(matching_scans) // 10)
    for idx, (scan_num, mz, rt) in enumerate(matching_scans, 1):
        try:
            if idx % progress_interval == 0:
                print(f"  Progress: {idx}/{len(matching_scans)} scans ({100*idx//len(matching_scans)}%)")
            
            # Load spectrum directly from cached raw_reader (NO FILE REOPEN)
            mz_array, intensity_array, charge_array, _ = raw_reader.get_scan_from_scan_number(scan_num)
            
            all_mz.extend(mz_array)
            all_intensity.extend(intensity_array)
            # Store RT for each peak from this scan
            all_scan_rt.extend([rt] * len(mz_array))
        except Exception as e:
            print(f"  Warning: Could not load scan {scan_num}: {e}")
            continue
    
    if not all_mz:
        raise ValueError("No MS/MS spectra could be loaded!")
    
    # Combine raw peaks into DataFrame
    combined_spectrum = pd.DataFrame({
        'm/z': all_mz,
        'Intensity': all_intensity,
        'RT': all_scan_rt
    })

    # Deduplicate experimental peaks by clustering peaks within 20 ppm
    cluster_tolerance_ppm = 20.0

    # Sort by m/z for clustering
    combined_spectrum = combined_spectrum.sort_values('m/z').reset_index(drop=True)

    clusters = []  # list of lists of row indices
    if len(combined_spectrum) == 0:
        clusters = []
    else:
        current_cluster = [0]
        current_ref = combined_spectrum.at[0, 'm/z']

        for idx in range(1, len(combined_spectrum)):
            mz = combined_spectrum.at[idx, 'm/z']
            ppm_diff = abs(mz - current_ref) / current_ref * 1e6 if current_ref != 0 else float('inf')
            if ppm_diff <= cluster_tolerance_ppm:
                current_cluster.append(idx)
                # update cluster reference as mean m/z
                current_ref = combined_spectrum.loc[current_cluster, 'm/z'].mean()
            else:
                clusters.append(current_cluster)
                current_cluster = [idx]
                current_ref = mz

        if current_cluster:
            clusters.append(current_cluster)

    # Build deduplicated peaks: keep max-intensity representative and compute area (integration over RT)
    dedup_rows = []
    for cluster in clusters:
        dfc = combined_spectrum.loc[cluster]
        # Representative is row with max intensity
        idx_max = dfc['Intensity'].idxmax()
        rep_mz = float(dfc.loc[idx_max, 'm/z'])
        rep_int = float(dfc.loc[idx_max, 'Intensity'])
        rep_rt = float(dfc.loc[idx_max, 'RT'])

        # Compute chromatographic area: first try Gaussian fit (scipy.curve_fit),
        # fall back to trapezoidal integration when scipy unavailable or fit fails.
        try:
            rt_vals = dfc['RT'].values
            int_vals = dfc['Intensity'].values
            # sort by RT
            order = rt_vals.argsort()
            rt_sorted = rt_vals[order]
            int_sorted = int_vals[order]

            peak_area = float(rep_int)

            # Only attempt fit when we have at least 3 points with RT spread
            if len(rt_sorted) >= 3 and (rt_sorted.max() - rt_sorted.min()) > 1e-6:
                try:
                    from scipy.optimize import curve_fit

                    def gaussian(x, A, mu, sigma, baseline):
                        return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + baseline

                    # initial parameters
                    A0 = float(int_sorted.max() - np.median(int_sorted))
                    if A0 <= 0:
                        A0 = float(int_sorted.max())
                    mu0 = float(rt_sorted[np.argmax(int_sorted)])
                    sigma0 = max(0.01, float((rt_sorted.max() - rt_sorted.min()) / 4.0))
                    baseline0 = float(np.median(int_sorted))

                    p0 = [A0, mu0, sigma0, baseline0]

                    # bounds: A>0, sigma>0
                    lower = [0.0, rt_sorted.min(), 1e-6, 0.0]
                    upper = [np.inf, rt_sorted.max(), np.inf, np.inf]

                    popt, pcov = curve_fit(gaussian, rt_sorted, int_sorted, p0=p0, bounds=(lower, upper), maxfev=2000)

                    A_fit, mu_fit, sigma_fit, baseline_fit = popt

                    # Validate fit parameters
                    if sigma_fit > 0 and A_fit > 0:
                        # Area of fitted Gaussian (analytic): A * sigma * sqrt(2*pi)
                        peak_area = float(A_fit * sigma_fit * np.sqrt(2.0 * np.pi))
                    else:
                        # fallback to trapezoid if fit is nonsensical
                        if len(rt_sorted) > 1:
                            peak_area = float(np.trapz(int_sorted, rt_sorted))
                        else:
                            peak_area = float(rep_int)

                except Exception:
                    # scipy not available or fit failed — use trapezoidal integration as fallback
                    if len(rt_sorted) > 1:
                        peak_area = float(np.trapz(int_sorted, rt_sorted))
                    else:
                        peak_area = float(rep_int)
            else:
                # Not enough points for fitting — use trapezoid or single-point fallback
                if len(rt_sorted) > 1:
                    peak_area = float(np.trapz(int_sorted, rt_sorted))
                else:
                    peak_area = float(rep_int)

        except Exception:
            peak_area = float(rep_int)

        dedup_rows.append({
            'm/z': rep_mz,
            'Intensity': rep_int,
            'RT': rep_rt,
            'Area': peak_area,
            'Count': int(len(dfc))
        })

    dedup_df = pd.DataFrame(dedup_rows)

    # If no peaks found, keep empty arrays
    if dedup_df.empty:
        exp_spectrum = (np.array([]), np.array([]))
        exp_rt_map = {}
        exp_area_map = {}
    else:
        # Sort deduped peaks by intensity descending for matching convenience
        dedup_df = dedup_df.sort_values('Intensity', ascending=False).reset_index(drop=True)
        exp_spectrum = (dedup_df['m/z'].values, dedup_df['Intensity'].values)
        exp_rt_map = dict(zip(dedup_df['m/z'].values, dedup_df['RT'].values))
        exp_area_map = dict(zip(dedup_df['m/z'].values, dedup_df['Area'].values))

    print(f"[OK] Deduplicated experimental peaks: {len(combined_spectrum)} raw -> {len(dedup_df)} unique (20 ppm clusters)")
    
    print(f"[OK] Loaded {len(exp_spectrum[0])} unique peaks from MS/MS")
    
    # Step 3: Generate theoretical fragments using updated method
    print(f"\n[Step 3] Generating theoretical fragments...")
    try:
        # Create analysis object (internally creates Glycopeptide)
        analysis = GlycoPeptideAnalysis(
            peptide_sequence=peptide_sequence,
            glycan_code=glycan_code,
            glycosylation_site=glycosylation_site,
            glycan_type=glycan_type,
            use_cam=True,
            max_structures=100,
            mod_string=mod_string,
            peptide_fragment_types=peptide_fragment_types,
            glycan_fragment_types=glycan_fragment_types
        )
        
        # All fragments are auto-generated during GlycoPeptideAnalysis initialization
        # Access raw theoretical fragments
        theoretical_table = analysis.theoretical_df
        total_fragments = len(theoretical_table)
        print(f"[OK] Generated {total_fragments} theoretical fragments")
        print(f"   Glycosylation site: position {glycosylation_site} ({peptide_sequence[glycosylation_site-1]})")
        
        # Calculate masses for QC summary (extracted from analysis object)
        peptide_mass = analysis.peptide_mass
        glycan_mass = analysis.glycan_mass
        glycopeptide_mass = analysis.mass

        # Print summary of calculated masses and modifications
        print("\n===== GLYCOPEPTIDE MASS SUMMARY =====")
        print(f"Peptide Sequence: {peptide_sequence}")
        print(f"Glycan Code: {glycan_code}")
        print(f"Glycan Type: {glycan_type}")
        print(f"Glycosylation Site: {glycosylation_site}")
        print(f"Modifications: {analysis.mod_string or ''}")
        print(f"Peptide Mass (Da): {peptide_mass}")
        print(f"Glycan Mass (Da): {glycan_mass}")
        print(f"Glycopeptide Mass (Da): {glycopeptide_mass}")
        print("=====================================")
        
    except Exception as e:
        raise ValueError(f"Failed to generate fragments: {e}")
    
    
    # Step 4: Match fragments  (all tables auto-generated)
    print(f"\n[Step 4] Matching fragments (PPM tolerance: ±{ppm_tolerance}, Charge states: 1-3)...")
    
    matched = analysis.match_fragments(
        experimental_spectrum=exp_spectrum,
        exp_rt_map=exp_rt_map,
        exp_area_map=exp_area_map,
        ppm_tolerance=ppm_tolerance,
        charge_states=[1, 2, 3]
    )
    print(f"[OK] Matched {len(matched)} fragments")
    
    # Step 5: Export results to 3-sheet Excel
    print(f"\n[Step 5] Exporting results...")
    
    if output_file is None:
        output_file = f"Analysis_{peptide_sequence}_{glycan_code}_{precursor_mz:.2f}.xlsx"
    
    # Use analysis.export_to_excel() - much simpler!
    analysis.export_to_excel(output_file, matched)
    
    if best_structure_report and len(matched) > 0:
        # Visualize predicted structures using the updated tree-based visualization
        analysis.visualize_predicted_structures(
            matched_df=matched,
            which="best",
            output_dir=structure_output_dir,
            output_excel=output_file,
            figsize=figsize
        )
   
    print("=" * 80)
    print(f"Analysis complete! Results saved to: {output_file}")
    print("=" * 80)
    
    return matched

# ============================================================================
# Interactive Main
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("GLYCOPEPTIDE MS/MS ANALYSIS - SERUM PRM SAMPLE")
    print("=" * 80)
    
    # ========== HARDCODED PARAMETERS ==========
    # Modify these values to analyze different samples
    
    raw_file = "Serum_PRM_1_6mz_1.raw"  # Path to Thermo RAW file
    
    # Precursor targeting
    precursor_mz = 1500.9452 # Target m/z
    charge = 4               # Charge state
    target_rt = 65.6          # Target retention time (minutes)
    rt_tolerance = 1.0       # RT search window (minutes)
    
    # Glycopeptide information
    peptide_sequence = "SPYYNVSDEISFHCYDGYTLR"  # Peptide sequence
    glycan_code = "HexNAc(4)Hex(5)Fuc(2)"          # Glycan composition code (HexNAc=2, Hex=1, NeuAc=1, Fuc=0) - valid O-glycan
    glycosylation_site = 5         # Position of N in N-X-S/T motif (1-indexed)
    glycan_type = 'N'              # N-glycan or O-glycan
    mod_string = None              # e.g., "N-term:Ac"
    peptide_fragment_types=["by"]
    glycan_fragment_types=["BY"]
    # Fragment matching tolerance
    ppm_tolerance = 20.0     # Fragment m/z tolerance (ppm)
    
    # Optional: Custom output file name
    output_file = f"{glycan_code}_{peptide_sequence}.xlsx"

    # Optional: Structure visualization and reporting
    best_structure_report = True  # Predict and visualize structures
    structure_output_dir = "output/"  # Directory to save structure images (None = display interactively)
    figsize = (5, 4)  # Figure size in inches (width, height) - default is (5, 4)
    
    print(f"\nAnalysis Parameters:")
    print(f"  RAW File: {raw_file}")
    print(f"  Precursor: m/z={precursor_mz:.4f}, z=+{charge}")
    print(f"  RT Range: {target_rt-rt_tolerance:.2f}-{target_rt+rt_tolerance:.2f} min")
    print(f"  Peptide: {peptide_sequence}")
    print(f"  Glycan: {glycan_code} ({glycan_type}-glycan)")
    print(f"  Glycosylation Site: Position {glycosylation_site} ({peptide_sequence[glycosylation_site-1]})")
    print(f"  Fragment Tolerance: ±{ppm_tolerance} ppm")
    
    # Run analysis
    try:
        results = analyze_glycopeptide_from_raw(
            raw_file,
            precursor_mz,
            charge,
            target_rt,
            rt_tolerance,
            peptide_sequence,
            glycan_code,
            glycosylation_site,
            glycan_type=glycan_type,
            ppm_tolerance=ppm_tolerance,
            output_file=output_file,
            mod_string=mod_string,
            best_structure_report=best_structure_report,
            structure_output_dir=structure_output_dir,
            peptide_fragment_types=peptide_fragment_types,
            glycan_fragment_types=glycan_fragment_types,
            figsize=figsize
        )
        # Normalize results: handle DataFrame or list-like return values
        result_list = []
        if results is None:
            result_list = []
        elif hasattr(results, 'to_dict') and callable(getattr(results, 'to_dict')):
            try:
                result_list = results.to_dict('records')
            except Exception:
                try:
                    result_list = list(results)
                except Exception:
                    result_list = []
        else:
            try:
                result_list = list(results)
            except Exception:
                result_list = []

        # Print summary of matched fragments
        if result_list:
            print(f"\n{'=' * 80}")
            print(f"SUMMARY OF MATCHED FRAGMENTS (Top 15)")
            print(f"{'=' * 80}")
            
            # Show top 15 matched fragments (safe key access)
            for i, frag in enumerate(result_list[:15], 1):
                frag_name = frag.get('Fragment Name') or frag.get('Fragment') or frag.get('name') or 'Unknown'
                charge = frag.get('Charge', 1)
                composition = frag.get('Composition', frag.get('composition', 'N/A'))

                theo = frag.get('Theoretical m/z')
                exp = frag.get('Experimental m/z')
                ppm = frag.get('PPM Error')
                intensity = frag.get('Intensity (%)')
                raw_int = frag.get('Raw Intensity')
                rt = frag.get('RT (min)', frag.get('RT', 0.0))

                print(f"\n{i:2d}. {frag_name} (z={charge})")
                print(f"     Composition: {composition}")

                if isinstance(theo, (int, float)):
                    print(f"     Theoretical m/z: {theo:.6f}")
                else:
                    print(f"     Theoretical m/z: {theo}")

                if isinstance(exp, (int, float)):
                    print(f"     Experimental m/z: {exp:.6f}")
                else:
                    print(f"     Experimental m/z: {exp}")

                if isinstance(ppm, (int, float)):
                    print(f"     PPM Error: {ppm:.2f} ppm")
                else:
                    print(f"     PPM Error: {ppm}")

                if isinstance(intensity, (int, float)):
                    raw_str = f" (raw: {int(raw_int)})" if raw_int is not None else ""
                    print(f"     Intensity: {intensity:.1f}%{raw_str}")
                else:
                    print(f"     Intensity: {intensity}")

                if isinstance(rt, (int, float)):
                    print(f"     RT: {rt:.3f} min")
                else:
                    print(f"     RT: {rt}")
            
            print(f"\n{'=' * 80}")
            print(f"Total matched fragments: {len(result_list)}")
            print(f"Results saved to: {output_file}")
        else:
            print("\n[ERROR] No fragments matched!")
    
    except FileNotFoundError:
        print(f"\n[ERROR] RAW file not found: {raw_file}", flush=True)
        print("   Please ensure the file is in the current directory or provide full path.", flush=True)
    
    except Exception as e:
        print(f"\n[ERROR] Error: {e}", flush=True)
        import traceback
        traceback.print_exc()


