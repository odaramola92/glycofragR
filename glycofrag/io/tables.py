"""
Fragment table generation and export utilities.

Provides DataFrame helpers for converting fragment dictionaries into
analysis-ready tables with m/z values and metadata.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.modifications import get_modification_type


MONO_KEYS = ("HexNAc", "Hex", "Fuc", "NeuAc", "NeuGc")


def _normalize_mono_separators(label: str) -> str:
    """Insert hyphens between adjacent monosaccharide tokens in a label."""
    if not label:
        return label

    import re

    token_pattern = r"(HexNAc|NeuAc|NeuGc|Hex|Fuc)(\d*)"
    matches = list(re.finditer(token_pattern, label))
    if len(matches) < 2:
        return label

    parts: List[str] = []
    last_end = 0
    for idx, match in enumerate(matches):
        parts.append(label[last_end:match.start()])
        parts.append(match.group(0))

        next_match = matches[idx + 1] if idx + 1 < len(matches) else None
        if next_match is not None and match.end() == next_match.start():
            parts.append("-")

        last_end = match.end()

    parts.append(label[last_end:])
    return "".join(parts)


def format_fragment_string(composition: Dict[str, int], ion_type: str, modification_type: int = 3) -> str:
    """Create a consistent fragment string with proper naming convention.
    
    Naming convention:
        - Pure reducing end (Y, YY, YYY, Z, ZZ, ZZZ): TYPE-composition-suffix
        - Internal fragments (BY, CZ, BYY, CZZ): Show both ends, no peptide
        - Non-reducing (B, C): composition-TYPE
    
    For pure reducing end ions, appends modification-specific suffix:
        - '-Redend': Reduced end (PERMETHYLATED_REDUCED=3 or REDUCED=1)
        - '-FreeEnd': Free end (PERMETHYLATED_FREE=2 or FREE=0)
        - '-2AB': 2-AB labeled (TWO_AB=4 or TWO_AB_PERMETHYLATED=5)
        - '-PEP': Glycopeptide (GLYCOPEPTIDE=6)
    
    Args:
        composition: Dict of monosaccharide counts
        ion_type: Type of ion ('B', 'Y', 'BY', 'YY', 'YYY', 'Z', 'ZZ', 'ZZZ', 'CZ', etc.)
        modification_type: Modification type (0-6, default 3 for PERMETHYLATED_REDUCED)
    """
    parts: List[str] = []
    for mono in MONO_KEYS:
        count = composition.get(mono, 0)
        if count > 0:
            parts.append(mono if count == 1 else f"{mono}{count}")
    comp_str = "-".join(parts)
    
    ion_type_upper = ion_type.upper()
    
    # Pure reducing end ion types (Y and Z series): use PREFIX format with -PEP/suffix
    if ion_type_upper in ('Y', 'YY', 'YYY', 'YYYY', 'YYYYY', 'Z', 'ZZ', 'ZZZ', 'ZZZZ', 'ZZZZZ'):
        suffix = _get_y_ion_suffix(modification_type)
        if comp_str:
            return f"{ion_type}-{comp_str}{suffix}"
        else:
            return f"{ion_type}{suffix}"
    
    # Internal combination ions (BY and CZ): show both ends, no peptide
    elif ion_type_upper == 'BY':
        # BY: Y cleavage and B cleavage (internal fragment, no reducing end peptide)
        if comp_str:
            return f"Y-{comp_str}-B"
        else:
            return f"Y--B"
    
    elif ion_type_upper == 'CZ':
        # CZ: Z cleavage and C cleavage (internal fragment, no reducing end peptide)
        if comp_str:
            return f"Z-{comp_str}-C"
        else:
            return f"Z--C"
    
    # Non-reducing ion types (B and C series): use SUFFIX format
    else:
        if comp_str:
            return f"{comp_str}-{ion_type}"
        else:
            return f"-{ion_type}"

def _get_y_ion_suffix(modification_type: int) -> str:
    """
    Get the appropriate Y-ion suffix based on modification type.
    
    Args:
        modification_type: Modification type (0-7)
            0: FREE
            1: REDUCED
            2: PERMETHYLATED_FREE
            3: PERMETHYLATED_REDUCED
            4: TWO_AB
            5: TWO_AB_PERMETHYLATED
            6: GLYCOPEPTIDE
            7: CUSTOM
    
    Returns:
        String suffix like '-RedEnd', '-FreeEnd', '-2AB', '-PEP', or '-Custom'
    """
    if modification_type in (1, 3):  # REDUCED or PERMETHYLATED_REDUCED
        return "-RedEnd"
    elif modification_type in (0, 2):  # FREE or PERMETHYLATED_FREE
        return "-FreeEnd"
    elif modification_type in (4, 5):  # TWO_AB or TWO_AB_PERMETHYLATED
        return "-2AB"
    elif modification_type == 6:  # GLYCOPEPTIDE
        return "-PEP"
    elif modification_type == 7:  # CUSTOM
        return "-Custom"
    else:
        # Default to RedEnd for unknown types
        return "-RedEnd"

def extract_fragment_composition(fragment_string: str) -> Dict[str, int]:
    """Extract composition from a fragment label into a dictionary."""
    if not fragment_string:
        return {}

    import re

    # Match monosaccharide tokens anywhere in the fragment string.
    # Order matters: longer names must come before shorter ones (e.g., HexNAc before Hex).
    token_pattern = r"(HexNAc|NeuAc|NeuGc|Hex|Fuc)(\d*)"
    tokens = re.findall(token_pattern, fragment_string)
    if not tokens:
        return {}

    composition: Dict[str, int] = {}
    for mono, count_str in tokens:
        count = int(count_str) if count_str else 1
        composition[mono] = composition.get(mono, 0) + count

    return composition

def format_fragment_name(
    fragment: Dict[str, Any],
    frag_type: str,
    peptide_position: Optional[int] = None,
    modification_type: int = 3,
    is_peptide: bool = False
) -> str:
    """
    Format a fragment name according to GlyPRM conventions.
    
    For glycan fragments:
        B-type: {composition}-B
        Y-type: Y-{composition}{suffix}
        BY-type: {composition}-BY
        YY-type: YY-{composition}{suffix}
        BYY-type: {composition}-BYY
        YYY-type: YYY-{composition}{suffix}
        BYYY-type: {composition}-BYYY
        C-type: {composition}-C
        Z-type: Z-{composition}{suffix}
        CZ-type: {composition}-CZ
        ZZ-type: ZZ-{composition}{suffix}
        CZZ-type: {composition}-CZZ
        ZZZ-type: ZZZ-{composition}{suffix}
        A-type : {composition}-A
        Custom: {_custom_label}
    
    Naming convention:
        - Reducing end ions (Y, YY, YYY, Z, ZZ, ZZZ): PREFIX format TYPE-composition-PEP
        - Non-reducing ions (B, BY, BYY, C, CZ, CZZ): SUFFIX format composition-TYPE
    
    For peptide fragments:
        b-type: b{position}
        y-type: y{position}
        c-type: c{position}
        z-type: z{position}
    
    Args:
        fragment: Fragment dictionary with monosaccharide counts and metadata
        frag_type: Fragment type code ('b_ions', 'y_ions', 'by_ions', etc.) or
                   type category ('B', 'Y', 'BY', 'YY', etc.)
        peptide_position: Position for peptide fragments
        modification_type: Modification type for Y-ion suffix
        is_peptide: Whether this is a peptide fragment
    
    Returns:
        Formatted fragment name string
    """
    # Clean up frag_type: remove _ions suffix and convert to uppercase
    # This allows the function to accept both 'b_ions' and 'B' formats
    frag_type = frag_type.replace('_ions', '').upper()
    
    # Handle custom labels with special processing based on fragment type
    if '_custom_label' in fragment:
        label = fragment['_custom_label']
        
        # For A ions, add -A suffix if not already present
        if frag_type == 'A' and not label.endswith('-A'):
            return f"{label}-A"
        
        # For custom ions, remove Y- or B- prefixes
        if frag_type == 'CUSTOM':
            if label.startswith('Y-'):
                return label[2:]
            elif label.startswith('B-'):
                return label[2:]
        
        return label
    
    # Handle special intact fragment
    if frag_type == 'INTACT':
        return "Intact"
    
    # Handle Y0 and Y1 glycopeptide fragments (check these before generic peptide check)
    if frag_type == 'Y0':
        if peptide_position is not None:
            return f"Y0-y{peptide_position}"
        return "Y0"
    
    if frag_type == 'Y1':
        # Check if fucosylated
        fuc_count = fragment.get('Fuc', 0)
        if fuc_count > 0:
            return "Y1F"
        return "Y1"
    
    if frag_type == 'Y_GLYCAN':  # y_glycan_ions becomes Y_GLYCAN after cleaning
        # Get composition name
        comp_str = _format_composition_string(fragment)
        return f"Y-{comp_str}"
    
    # Handle peptide fragments (only for simple letter types without composition)
    if is_peptide or frag_type in ['b', 'y', 'c', 'z', 'B', 'Y', 'C', 'Z']:
        comp_str = _format_composition_string(fragment)
        # If it has composition, it's likely a glycan fragment, not peptide
        if comp_str:
            # Fall through to glycan handler
            pass
        else:
            # No composition, so treat as simple peptide fragment
            if peptide_position is not None:
                ion_type = frag_type.lower() if frag_type.isupper() else frag_type
                return f"{ion_type}{peptide_position}"
            # Return uppercase version if that's what was passed in, otherwise lowercase
            return frag_type.upper() if frag_type.isupper() else frag_type.lower()
    
    # Handle glycan fragment types - frag_type is already cleaned at the top of the function
    comp_str = _format_composition_string(fragment)
    
    if not comp_str:
        return frag_type  # Return just the type if no composition
    
    # Apply proper format based on type
    if frag_type == 'B':
        return f"{comp_str}-B"
    elif frag_type == 'Y':
        suffix = _get_y_ion_suffix(modification_type)
        return f"Y-{comp_str}{suffix}"
    elif frag_type == 'BY':
        # BY fragments: Y-composition-B (internal fragment, no -PEP)
        return f"Y-{comp_str}-B"
    elif frag_type == 'YY':
        # YY fragments are reducing end, use prefix format
        suffix = _get_y_ion_suffix(modification_type)
        return f"YY-{comp_str}{suffix}"
    elif frag_type == 'BYY':
        # BYY: YY-type internal fragment with B-cleavage
        return f"YY-{comp_str}-B"
    elif frag_type == 'YYY':
        # YYY: pure reducing end (multiple cleavages on arms)
        suffix = _get_y_ion_suffix(modification_type)
        return f"YYY-{comp_str}{suffix}"
    elif frag_type == 'BYYY':
        # BYYY: YYY-type internal fragment with B-cleavage
        return f"YYY-{comp_str}-B"
    elif frag_type == 'C':
        # C-ions formatted like B-ions (non-reducing end)
        return f"{comp_str}-C"
    elif frag_type == 'Z':
        # Z-ions are reducing end like Y-ions
        suffix = _get_y_ion_suffix(modification_type)
        return f"Z-{comp_str}{suffix}"
    elif frag_type == 'CZ':
        # CZ: Z-type internal fragment with C-cleavage
        return f"Z-{comp_str}-C"
    elif frag_type == 'ZZ':
        # ZZ: pure reducing end (multiple cleavages on arms)
        suffix = _get_y_ion_suffix(modification_type)
        return f"ZZ-{comp_str}{suffix}"
    elif frag_type == 'CZZ':
        # CZZ: ZZ-type internal fragment with C-cleavage
        return f"ZZ-{comp_str}-C"
    elif frag_type == 'ZZZ':
        # ZZZ: pure reducing end (multiple cleavages on arms)
        suffix = _get_y_ion_suffix(modification_type)
        return f"ZZZ-{comp_str}{suffix}"
    elif frag_type == 'CZZZ':
        # CZZZ: ZZZ-type internal fragment with C-cleavage
        return f"ZZZ-{comp_str}-C"
    elif frag_type == 'BZZZ':
        # BZZZ: ZZZ-type internal fragment with B-cleavage
        return f"ZZZ-{comp_str}-B"
    elif frag_type == 'A':
        # A-type oxonium ions
        return f"{comp_str}-A"
    else:
        # Default format
        return f"{comp_str}-{frag_type}"

def _format_composition_string(composition: Dict[str, Any]) -> str:
    """
    Format a composition dictionary into a readable string.
    
    Example: {'HexNAc': 2, 'Hex': 3} -> "HexNAc2-Hex3"
    
    Args:
        composition: Dictionary with monosaccharide counts
    
    Returns:
        Formatted composition string
    """
    if not composition:
        return ""
    
    parts: List[str] = []
    for mono in MONO_KEYS:
        count = composition.get(mono, 0)
        if count and count > 0:
            if count == 1:
                parts.append(mono)
            else:
                parts.append(f"{mono}{count}")
    
    return "-".join(parts)

def clean_composition_for_display(composition_str: str, frag_type: str) -> str:
    """
    Clean composition string for display in tables.
    
    Removes prefixes like 'Y-' and 'B-' for custom ions.
    Removes fragment type suffixes like '-B', '-Y', '-A', '-By', etc.
    
    Args:
        composition_str: The composition string to clean
        frag_type: Fragment type (will be normalized to uppercase without _ions)
    
    Returns:
        Cleaned composition string
    """
    # Normalize frag_type
    frag_type_clean = frag_type.replace('_ions', '').upper()
    
    # For custom ions, remove Y- and B- prefixes
    if frag_type_clean == 'CUSTOM':
        if composition_str.startswith('Y-'):
            composition_str = composition_str[2:]
        elif composition_str.startswith('B-'):
            composition_str = composition_str[2:]
    
    # Remove fragment type suffixes for all types
    # Remove -B, -Y, -A, -RedEnd, -FreeEnd, -2AB, -PEP suffixes
    if composition_str.endswith('-B'):
        composition_str = composition_str[:-2]
    elif composition_str.endswith('-Y'):
        composition_str = composition_str[:-2]
    elif composition_str.endswith('-A'):
        composition_str = composition_str[:-2]
    elif composition_str.endswith('-RedEnd'):
        composition_str = composition_str[:-7]
    elif composition_str.endswith('-FreeEnd'):
        composition_str = composition_str[:-8]
    elif composition_str.endswith('-2AB'):
        composition_str = composition_str[:-4]
    elif composition_str.endswith('-PEP'):
        composition_str = composition_str[:-4]
    
    return _normalize_mono_separators(composition_str)

def _format_glycan_composition(fragment: Mapping[str, Any]) -> str:
    glycan_comp = fragment.get("glycan_composition")
    if isinstance(glycan_comp, dict):
        return _format_composition_string(glycan_comp)
    if isinstance(fragment, dict):
        return _format_composition_string(fragment)
    return ""

def _as_mono_dict(fragment: Mapping[str, Any]) -> Dict[str, int]:
    mono_dict: Dict[str, int] = {}
    for mono in MONO_KEYS:
        value = fragment.get(mono)
        if isinstance(value, int):
            mono_dict[mono] = value
    return mono_dict

def _format_peptide_sequence_with_glycan(sequence: str, glycan_comp: str) -> str:
    if not sequence:
        return glycan_comp
    if not glycan_comp:
        return sequence
    return f"{sequence}-{glycan_comp}"

def _format_fragment_fields(
    fragment: Mapping[str, Any],
    fragment_type_clean: str,
    modification_type_int: int = 6,
) -> Tuple[str, str, str, str]:
    """
    Return (fragment_name, composition, type_label, fragment_type_series).
    
    Args:
        fragment: Fragment dictionary
        fragment_type_clean: Fragment type (B, Y, BY, etc.)
        modification_type_int: Modification type as integer (0-6). Only glycopeptide mode (6) gets -PEP suffix
    """
    type_label = "Glycan"
    fragment_type_series = fragment_type_clean
    fragment_name = str(fragment.get("name", ""))
    composition = ""

    if fragment_type_clean.startswith("PEPTIDE_"):
        type_label = "Peptide"
        fragment_name = str(fragment.get("fragment_name", fragment_name))
        composition = str(fragment.get("fragment_sequence", fragment.get("sequence", "")))
        series_value = fragment.get("fragment_type", fragment_type_clean)
        fragment_type_series = str(series_value).lower()
        # Strip neutral loss suffix from Fragment Type (e.g. 'b-co' -> 'b')
        if '-' in fragment_type_series:
            fragment_type_series = fragment_type_series.split('-')[0]
        # Composition already includes loss label from fragment_sequence (e.g. 'LCPD(-H2O)')
        return fragment_name, composition, type_label, fragment_type_series

    if fragment_type_clean in ("Y0", "Y1", "Y1F", "Y_GLYCAN"):
        type_label = "Peptide" if fragment_type_clean == "Y0" else "Glycopeptide"
        fragment_type_series = "y"
        if fragment_type_clean in ("Y1", "Y1F"):
            # Y1 is the intact peptide + core HexNAc — always "Y1"
            fragment_name = "Y1F" if fragment_type_clean == "Y1F" else "Y1"
        else:
            position = fragment.get("position")
            if isinstance(position, int):
                fragment_name = f"y{position}"
            else:
                fragment_name = str(fragment.get("peptide_name", fragment_name))
        sequence = str(fragment.get("sequence", ""))
        glycan_comp = _format_glycan_composition(fragment)
        if fragment.get("glycan_attached"):
            composition = _format_peptide_sequence_with_glycan(sequence, glycan_comp)
        else:
            composition = sequence or "Peptide"
        return fragment_name, composition, type_label, fragment_type_series

    if fragment_type_clean == "INTACT":
        type_label = "Glycopeptide"
        fragment_type_series = "intact"
        fragment_name = "Intact"
        sequence = str(fragment.get("sequence", ""))
        glycan_comp = _format_glycan_composition(fragment)
        composition = _format_peptide_sequence_with_glycan(sequence, glycan_comp)
        return fragment_name, composition, type_label, fragment_type_series

    # A-type oxonium ions and custom diagnostic ions
    if fragment_type_clean in ("A", "CUSTOM"):
        type_label = "Diagnostic"
        fragment_type_series = "custom"
        custom_label = fragment.get("_custom_label", "")
        if custom_label:
            if fragment_type_clean == "A":
                fragment_name = f"Y-{custom_label}-B"
            else:
                fragment_name = custom_label
            fragment_name = _normalize_mono_separators(fragment_name)
            composition = _normalize_mono_separators(
                clean_composition_for_display(fragment_name, "CUSTOM")
            )
        else:
            glycan_comp = _format_glycan_composition(fragment)
            if fragment_type_clean == "A":
                fragment_name = f"Y-{glycan_comp}-B" if glycan_comp else "Y--B"
            else:
                fragment_name = f"{glycan_comp}-{fragment_type_clean}" if glycan_comp else fragment_type_clean
            fragment_name = _normalize_mono_separators(fragment_name)
            composition = _normalize_mono_separators(glycan_comp or "")
        return fragment_name, composition, type_label, fragment_type_series

    # Glycan fragments (B/Y/BY/YY/etc.)
    # Pure reducing end (Y, YY, YYY, Z, ZZ, ZZZ): Type as PREFIX with -PEP in glycopeptide mode
    # Internal fragments (BY, CZ, BYY, CZZ, etc.): Show cleavage types, no peptide
    # Non-reducing (B, C): Type as SUFFIX
    pure_reducing_end = {'Y', 'YY', 'YYY', 'Z', 'ZZ', 'ZZZ'}
    glycan_comp = _format_glycan_composition(fragment)
    if glycan_comp:
        if fragment_type_clean in pure_reducing_end:
            # Pure Y/Z ions: Type as PREFIX with modification-specific suffix
            suffix = _get_y_ion_suffix(modification_type_int)
            fragment_name = f"{fragment_type_clean}-{glycan_comp}{suffix}"
            composition = f"{glycan_comp}{suffix}"  # Include suffix in composition
        elif fragment_type_clean == 'BY':
            # BY fragments: Y-{comp}-B (internal fragment, no peptide)
            fragment_name = f"Y-{glycan_comp}-B"
            composition = glycan_comp
        elif fragment_type_clean == 'CZ':
            # CZ fragments: Z-{comp}-C (internal fragment, no peptide)
            fragment_name = f"Z-{glycan_comp}-C"
            composition = glycan_comp
        elif fragment_type_clean == 'BYY':
            # BYY fragments: YY-{comp}-B (internal fragment derived from YY)
            fragment_name = f"YY-{glycan_comp}-B"
            composition = glycan_comp
        elif fragment_type_clean == 'CZZ':
            # CZZ fragments: ZZ-{comp}-C (internal fragment derived from ZZ)
            fragment_name = f"ZZ-{glycan_comp}-C"
            composition = glycan_comp
        elif fragment_type_clean == 'BYYY':
            # BYYY fragments: YYY-{comp}-B (internal fragment derived from YYY)
            fragment_name = f"YYY-{glycan_comp}-B"
            composition = glycan_comp
        elif fragment_type_clean == 'CZZZ':
            # CZZZ fragments: ZZZ-{comp}-C (internal fragment derived from ZZZ)
            fragment_name = f"ZZZ-{glycan_comp}-C"
            composition = glycan_comp
        elif fragment_type_clean == 'BZZZ':
            # BZZZ fragments: ZZZ-{comp}-B (internal fragment derived from ZZZ)
            fragment_name = f"ZZZ-{glycan_comp}-B"
            composition = glycan_comp
        else:
            # B, C, etc.: Type as SUFFIX
            fragment_name = f"{glycan_comp}-{fragment_type_clean}"
            composition = glycan_comp
    else:
        fragment_name = fragment_type_clean
        composition = ""

    fragment_name = _normalize_mono_separators(fragment_name)
    composition = _normalize_mono_separators(composition)

    return fragment_name, composition, type_label, fragment_type_series

def _deduplicate_fragments(
    data: List[Dict[str, object]]
) -> List[Dict[str, object]]:
    """
    Deduplicate fragments by m/z value, keeping first occurrence but combining structure info.
    
    When the same fragment appears from multiple structures:
    - Keep the first occurrence
    - Combine structure sources (e.g., "1,2,3")
    
    Args:
        data: List of fragment rows with 'Structure' field
        
    Returns:
        Deduplicated list with combined structure info
    """
    # Group by m/z and fragment name
    seen = {}  # key: (m/z, fragment_name) -> index
    combined_structures = {}  # key: (m/z, fragment_name) -> list of structures
    
    for row in data:
        mz_val = row.get("m/z(z=1)")
        frag_name = row.get("Fragment")
        structure = row.get("Structure", "unknown")
        
        key = (mz_val, frag_name)
        
        if key not in seen:
            seen[key] = len(seen)
            combined_structures[key] = [structure]
        else:
            # Track additional structures
            if structure not in combined_structures[key]:
                combined_structures[key].append(structure)
    
    # Collect all numeric structure indices used
    all_struct_nums = set()
    for structs in combined_structures.values():
        for s in structs:
            if s != 'all' and s != 'unknown':
                try:
                    all_struct_nums.add(int(s))
                except (ValueError, TypeError):
                    pass
    all_struct_label = ",".join(str(n) for n in sorted(all_struct_nums)) if all_struct_nums else "all"
    
    # Build deduplicated result
    result = []
    processed = set()
    
    for row in data:
        mz_val = row.get("m/z(z=1)")
        frag_name = row.get("Fragment")
        key = (mz_val, frag_name)
        
        if key not in processed:
            # Add all unique structures for this fragment
            structures = combined_structures[key]
            # Sort if all are integers
            try:
                structures = [str(int(s)) if s != 'all' else s for s in structures]
                structures = sorted(structures, key=lambda x: (x == 'all', int(x) if x != 'all' else 0))
            except (ValueError, TypeError):
                structures = sorted(structures)
            
            row_copy = row.copy()
            struct_str = ",".join(structures)
            # Replace 'all' with actual structure numbers for consistency
            if 'all' in structures and all_struct_nums:
                struct_str = all_struct_label
            row_copy["Structure"] = struct_str
            result.append(row_copy)
            processed.add(key)
    
    return result

def _iter_fragment_rows(
    fragments: Mapping[str, Sequence[Mapping[str, Any]]],
    glycan_code: str,
    calculator: GlycanMassCalculator,
    peptide: Optional[str],
    charge_states: Iterable[int],
    modification_type_int: int = 6,
    structure_label: Optional[str] = None,
) -> List[Dict[str, object]]:
    data: List[Dict[str, object]] = []
    charges = list(charge_states)

    for frag_type, frag_list in fragments.items():
        fragment_type_clean = frag_type.replace("_ions", "").upper()
        for fragment in frag_list:
            if not isinstance(fragment, Mapping):
                continue

            fragment_mass = fragment.get("mass")
            if not isinstance(fragment_mass, (int, float)):
                fragment_mass = fragment.get("fragment_mass")

            if not isinstance(fragment_mass, (int, float)):
                glycan_comp = fragment.get("glycan_composition")
                if isinstance(glycan_comp, dict):
                    fragment_mass = calculator.calculate_fragment_mass(
                        glycan_comp,
                        frag_type,
                        peptide=peptide,
                    )
                else:
                    fragment_mass = calculator.calculate_fragment_mass(
                        _as_mono_dict(fragment),
                        frag_type,
                        peptide=peptide,
                    )
            if fragment_mass is None:
                continue

            fragment_name, composition, type_label, fragment_type_series = _format_fragment_fields(
                fragment,
                fragment_type_clean,
                modification_type_int=modification_type_int,
            )
            if not fragment_name:
                continue

            row: Dict[str, object] = {
                "Fragment": fragment_name,
                "Composition": composition,
                "Type": type_label,
                "Fragment Type": fragment_type_series,
                "Mass(Da)": round(float(fragment_mass), 4),
                "Structure": str(fragment.get("structure", fragment.get("_structure_number"))) if (fragment.get("structure") is not None or fragment.get("_structure_number") is not None) else (structure_label or "unknown"),
            }

            for charge in charges:
                mz = calculator.calculate_mz(fragment_mass, charge)
                row[f"m/z(z={charge})"] = round(float(mz), 4)

            data.append(row)

    return data

def generate_fragment_table(
    fragments: Mapping[str, Sequence[Mapping[str, Any]]],
    glycan_code: str,
    modification_type: Union[str, int],
    peptide: Optional[str] = None,
    use_cam: bool = False,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    mod_string: Optional[str] = None,
    glycan_type: Optional[str] = None,
    charge_states: Iterable[int] = (1, 2, 3),
    deduplicate: bool = True,
    structure_label: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    """
    Generate a DataFrame of fragment ions with m/z values and metadata.
    
    Args:
        fragments: Dictionary mapping fragment type to list of fragment dicts
        glycan_code: Glycan composition code (e.g., '4501')
        modification_type: **REQUIRED**. Reducing end modification type that controls Y-ion labeling.
            Can be:
            - String (case-insensitive):
                - 'free': Free reducing end
                - 'reduced': Reduced (alditol) end
                - 'permethylated_free' or 'pm_free': Permethylated free end
                - 'permethylated_reduced' or 'pm_reduced': Permethylated reduced end
                - 'two_ab' or '2ab': 2-AB labeled
                - 'two_ab_permethylated' or '2ab_permethylated': 2-AB + permethylated
                - 'glycopeptide' or 'pep': Peptide-attached reducing end (use for Glycopeptide)
                - 'custom': Custom modification (use for standalone Glycan with custom masses)
            - Integer 0-7 for backward compatibility
            
            **MUST MATCH your workflow:**
            - For Glycopeptide: pass 'glycopeptide' (or 6) → Y-ions labeled with -PEP suffix
            - For standalone Glycan: pass the type used at Glycan initialization
              (0/'free', 3/'permethylated_reduced', 7/'custom', etc.)
              → Y-ions labeled with -FreeEnd, -RedEnd, -2AB, or -Custom suffix
            
        peptide: Peptide sequence for mass calculations
        use_cam: Whether to apply carbamidomethylation
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Peptide modification string (e.g., 'M:Ox; K:Ac')
        glycan_type: 'N' or 'O' for glycan type
        charge_states: Tuple of charge states to calculate m/z (default (1, 2, 3))
        deduplicate: Whether to deduplicate fragments by m/z (default True)
        structure_label: Optional structure identifier to populate the Structure column
    
    Returns:
        DataFrame with columns: Fragment, Composition, Type, Fragment Type, Mass(Da), Structure, m/z(z=1), m/z(z=2), m/z(z=3)
    """
    # Convert string modification_type to integer
    mod_type_int = get_modification_type(modification_type)
    
    calculator = GlycanMassCalculator(
        modification_type=mod_type_int,
        use_cam=use_cam,
        fixed_mods=fixed_mods or [],
        variable_mods=variable_mods or [],
        mod_string=mod_string,
        peptide=peptide,
    )

    data = _iter_fragment_rows(
        fragments=fragments,
        glycan_code=glycan_code,
        calculator=calculator,
        peptide=peptide,
        charge_states=charge_states,
        modification_type_int=mod_type_int,
        structure_label=str(structure_label) if structure_label is not None else None,
    )
    
    # Deduplicate fragments from multiple structures
    if deduplicate and data:
        data = _deduplicate_fragments(data)

    df = pd.DataFrame(data)
    return df

def generate_all_fragments_table(
    results: Dict[str, dict],
    modification_type: Union[str, int],
    peptide: Optional[str] = None,
    use_cam: bool = False,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    mod_string: Optional[str] = None,
    glycan_type: Optional[str] = None,
    charge_states: Iterable[int] = (1, 2, 3),
    structure_label: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    """
    Generate a DataFrame containing ALL fragments from each structure with identifiers.
    
    Args:
        results: Dict mapping glycan_code -> {"fragments": fragments_dict, ...}
        structure_label: Optional structure identifier for the Structure column
    """
    all_tables: List[pd.DataFrame] = []

    for glycan_code, result in results.items():
        fragments = result.get("fragments")
        if isinstance(fragments, tuple):
            fragments = fragments[0]

        if not isinstance(fragments, dict):
            continue

        df = generate_fragment_table(
            fragments=fragments,
            glycan_code=glycan_code,
            modification_type=modification_type,
            peptide=peptide,
            use_cam=use_cam,
            fixed_mods=fixed_mods,
            variable_mods=variable_mods,
            mod_string=mod_string,
            glycan_type=glycan_type,
            charge_states=charge_states,
            structure_label=(
                result.get("structure_label")
                or result.get("structure_number")
                or result.get("structure")
                or structure_label
            ),
        )

        if df.empty:
            continue

        all_tables.append(df)

    if not all_tables:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(
            columns=[
                "Fragment",
                "Composition",
                "Type",
                "Fragment Type",
                "Mass(Da)",
                "m/z(z=1)",
                "m/z(z=2)",
                "m/z(z=3)",
            ]
        )

    return pd.concat(all_tables, ignore_index=True)

def deduplicate_fragments_by_mass(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate fragments by Fragment_mz per Glycopeptide with basic prioritization.
    """
    if df.empty:
        return df

    df = df.copy()
    if "Fragment_mz" not in df.columns:
        return df

    df["Fragment_mz"] = df["Fragment_mz"].round(4)

    priority_map = {
        "B": 1,
        "Y": 2,
        "BY": 3,
        "YY": 4,
        "BYY": 5,
        "YYY": 6,
        "BYYY": 7,
        "C": 8,
        "Z": 9,
        "CZ": 10,
        "ZZ": 11,
        "CZZ": 12,
        "BZZ": 13,
        "ZZZ": 14,
        "CZZZ": 15,
        "BZZZ": 16,
        "B-PEP": 50,
        "Y-PEP": 51,
        "C-PEP": 52,
        "Z-PEP": 53,
    }

    if "FragmentType" in df.columns:
        df["fragment_priority"] = df["FragmentType"].map(priority_map).fillna(99)
    else:
        df["fragment_priority"] = 99

    if "Ions" in df.columns:
        df["charge_priority"] = (
            df["Ions"].astype(str).str.replace("H+", "", regex=False).astype(float)
        )
    else:
        df["charge_priority"] = 99

    sort_cols = ["fragment_priority", "charge_priority"]
    if "Glycopeptide" in df.columns:
        sort_cols = ["Glycopeptide", "Fragment_mz"] + sort_cols
    else:
        sort_cols = ["Fragment_mz"] + sort_cols

    # Sort ascending (lower priority number = keep first)
    df = df.sort_values(sort_cols, ascending=True)

    subset_cols = ["Fragment_mz"]
    if "Glycopeptide" in df.columns:
        subset_cols = ["Glycopeptide", "Fragment_mz"]

    deduped = df.drop_duplicates(subset=subset_cols, keep="first")
    return deduped.drop(columns=["fragment_priority", "charge_priority"], errors="ignore")

def export_fragment_table_to_excel(
    df: pd.DataFrame,
    df_all: Optional[pd.DataFrame] = None,
    filename: str = "glycan_fragments.xlsx",
) -> str:
    """
    Export fragment tables to Excel with optional all-fragments sheet.
    """
    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Fragments", index=False)
        if df_all is not None:
            df_all.to_excel(writer, sheet_name="All_Fragments", index=False)

    return f"Fragment tables exported to {filename}"


def create_theoretical_summary_sheet(
    peptide_sequence: str,
    glycan_code: str,
    peptide_mass: float,
    glycan_mass: float,
    glycopeptide_mass: float,
    modifications_applied: Optional[Dict[str, List[str]]] = None,
    use_cam: bool = True,
    mod_string: Optional[str] = None,
    glycosylation_site: Optional[int] = None,
    glycan_type: Optional[str] = None,
) -> pd.DataFrame:
    """
    Create a comprehensive summary sheet for theoretical glycopeptide analysis.
    
    This sheet provides users with clear information about:
    - Peptide sequence and calculated mass
    - Glycan code and calculated mass
    - Glycopeptide mass (peptide + glycan)
    - Modifications applied to the peptide and where
    - Troubleshooting info for debugging mass discrepancies
    
    Args:
        peptide_sequence: The peptide amino acid sequence
        glycan_code: Glycan composition code (e.g., '4501')
        peptide_mass: Calculated neutral mass of peptide
        glycan_mass: Calculated neutral mass of glycan
        glycopeptide_mass: Calculated neutral mass of combined glycopeptide
        modifications_applied: Dict mapping modification names to list of positions
        use_cam: Whether CAM was applied
        mod_string: Modification string used
        glycosylation_site: Position of glycan attachment (1-indexed)
        glycan_type: 'N' or 'O' for glycan type
    
    Returns:
        DataFrame with two sections: Summary Info and Modifications Detail
    """
    summary_rows = []
    
    # SECTION 1: BASIC INFORMATION
    summary_rows.append({
        "Category": "--- PEPTIDE INFO ---",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "Peptide",
        "Parameter": "Sequence",
        "Value": peptide_sequence
    })
    summary_rows.append({
        "Category": "Peptide",
        "Parameter": "Length",
        "Value": len(peptide_sequence)
    })
    summary_rows.append({
        "Category": "Peptide",
        "Parameter": "Calculated Mass (Da)",
        "Value": f"{peptide_mass:.6f}"
    })
    
    # SECTION 2: GLYCAN INFORMATION
    summary_rows.append({
        "Category": "",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "--- GLYCAN INFO ---",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "Glycan",
        "Parameter": "Composition Code",
        "Value": glycan_code
    })
    if glycan_type:
        summary_rows.append({
            "Category": "Glycan",
            "Parameter": "Type",
            "Value": f"{glycan_type}-glycan"
        })
    if glycosylation_site:
        site_aa = peptide_sequence[glycosylation_site - 1] if glycosylation_site <= len(peptide_sequence) else "?"
        summary_rows.append({
            "Category": "Glycan",
            "Parameter": "Attachment Site",
            "Value": f"Position {glycosylation_site} ({site_aa})"
        })
    summary_rows.append({
        "Category": "Glycan",
        "Parameter": "Calculated Mass (Da)",
        "Value": f"{glycan_mass:.6f}"
    })
    
    # SECTION 3: GLYCOPEPTIDE MASS
    summary_rows.append({
        "Category": "",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "--- CALCULATED MASSES ---",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "Combined",
        "Parameter": "Peptide Mass (Da)",
        "Value": f"{peptide_mass:.6f}"
    })
    summary_rows.append({
        "Category": "Combined",
        "Parameter": "Glycan Mass (Da)",
        "Value": f"{glycan_mass:.6f}"
    })
    summary_rows.append({
        "Category": "Combined",
        "Parameter": "Total Glycopeptide Mass (Da)",
        "Value": f"{glycopeptide_mass:.6f}"
    })
    
    # SECTION 4: MODIFICATIONS
    summary_rows.append({
        "Category": "",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "--- MODIFICATIONS APPLIED ---",
        "Parameter": "",
        "Value": ""
    })
    
    if use_cam:
        if 'C' in peptide_sequence.upper():
            c_positions = [i + 1 for i, aa in enumerate(peptide_sequence.upper()) if aa == 'C']
            summary_rows.append({
                "Category": "Modifications",
                "Parameter": "CAM (Carbamidomethyl)",
                "Value": f"Applied to C at positions {', '.join(map(str, c_positions))}"
            })
        else:
            summary_rows.append({
                "Category": "Modifications",
                "Parameter": "CAM Status",
                "Value": "NOT APPLIED (no cysteines in sequence)"
            })
    
    # Add other modifications from mod_string or modifications_applied
    if mod_string:
        summary_rows.append({
            "Category": "Modifications",
            "Parameter": "Additional Modifications",
            "Value": mod_string
        })
    
    if modifications_applied:
        for mod_name, positions in sorted(modifications_applied.items()):
            pos_str = ', '.join(positions)
            summary_rows.append({
                "Category": "Modifications",
                "Parameter": f"  {mod_name}",
                "Value": f"Positions: {pos_str}"
            })
    
    # SECTION 5: DEBUG NOTES
    summary_rows.append({
        "Category": "",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "--- DEBUG NOTES ---",
        "Parameter": "",
        "Value": ""
    })
    summary_rows.append({
        "Category": "Debug",
        "Parameter": "Expected Glycopeptide m/z (z=2)",
        "Value": f"{(glycopeptide_mass + 2 * 1.00727646688) / 2:.6f}"
    })
    summary_rows.append({
        "Category": "Debug",
        "Parameter": "Expected Glycopeptide m/z (z=3)",
        "Value": f"{(glycopeptide_mass + 3 * 1.00727646688) / 3:.6f}"
    })
    summary_rows.append({
        "Category": "Debug",
        "Parameter": "To troubleshoot mass issues:",
        "Value": "Check that peptide mass matches, glycan code is correct, and modifications are expected"
    })
    
    return pd.DataFrame(summary_rows)


def build_clean_theoretical(theoretical_df: pd.DataFrame, mass_round: int = 6) -> pd.DataFrame:
    """
    Merge theoretical fragments that have the same rounded mass AND identical Composition.

    Rules:
    - Only merge rows when the Composition values are identical (string-equal).
    - Keep the first occurrence's `Fragment` value.
    - Combine unique `Fragment Type` values into a comma-separated sorted string.
    - Do not change column names; returned DataFrame preserves columns but one row per merged group.

    Args:
        theoretical_df: Original theoretical DataFrame (must contain 'Mass(Da)' and 'Composition')
        mass_round: Number of decimal places to round mass for grouping

    Returns:
        Cleaned DataFrame
    """
    if theoretical_df is None or theoretical_df.empty:
        return theoretical_df.copy()

    if 'Mass(Da)' not in theoretical_df.columns or 'Composition' not in theoretical_df.columns:
        # Nothing to do
        return theoretical_df.copy()

    df = theoretical_df.copy()
    # Normalize Composition to string for comparison
    df['Composition_str'] = df['Composition'].astype(str)
    # Create mass key for grouping
    df['mass_key'] = df['Mass(Da)'].round(mass_round)

    grouped_rows = []

    # Group by mass_key, then by exact Composition string
    for mass_key, mass_grp in df.groupby('mass_key'):
        for comp_val, comp_grp in mass_grp.groupby('Composition_str'):
            # If composition is literal 'nan' or empty, treat as non-mergeable (keep rows individually)
            if comp_val in ("nan", "None", "NoneType", "") and len(comp_grp) == 1:
                grouped_rows.append(comp_grp.iloc[0].to_dict())
                continue

            if len(comp_grp) == 1:
                grouped_rows.append(comp_grp.iloc[0].to_dict())
                continue

            # Multiple rows with same mass_key and identical Composition -> merge
            # Keep first occurrence as base
            base = comp_grp.iloc[0].to_dict()

            # Combine Fragment Type values (unique, sorted)
            if 'Fragment Type' in comp_grp.columns:
                types = comp_grp['Fragment Type'].dropna().astype(str).unique().tolist()
            else:
                types = []

            # Remove empty placeholders
            types = [t for t in types if t and t.lower() != 'nan']
            types_sorted = sorted(types)
            combined_type = ','.join(types_sorted) if types_sorted else base.get('Fragment Type', '')
            base['Fragment Type'] = combined_type

            # For merged rows only: replace the 'Fragment' text with the Composition value
            try:
                comp_text = str(base.get('Composition', ''))
                base['Fragment'] = comp_text
            except Exception:
                pass

            # Other columns: leave as in base row
            grouped_rows.append(base)

    clean_df = pd.DataFrame(grouped_rows)

    # Drop helper cols if present
    if 'Composition_str' in clean_df.columns:
        clean_df = clean_df.drop(columns=['Composition_str'], errors=True)
    if 'mass_key' in clean_df.columns:
        clean_df = clean_df.drop(columns=['mass_key'], errors=True)

    # Preserve original column order where possible; ensure 'Composition' and 'Fragment' remain
    cols = [c for c in theoretical_df.columns if c in clean_df.columns]
    # Add any extra cols that might have appeared
    extra = [c for c in clean_df.columns if c not in cols]
    final_cols = cols + extra

    return clean_df[final_cols]


def _is_structure_id(value: object) -> bool:
    """Check if value is a single structure ID (digit string)."""
    if value is None:
        return False
    text = str(value).strip()
    return text.isdigit()


def _parse_structure_ids(value: object) -> set:
    """Parse structure ID(s) from a value that may be a single ID or comma-separated list.
    
    Examples:
        "1" -> {1}
        "1,2,3" -> {1, 2, 3}
        "1,2,3,4,5" -> {1, 2, 3, 4, 5}
    
    Returns:
        Set of integer structure IDs
    """
    if value is None:
        return set()
    
    text = str(value).strip()
    if not text or text.lower() == 'nan':
        return set()
    
    # Split by comma and parse each ID
    structure_ids = set()
    for part in text.split(','):
        part = part.strip()
        if part.isdigit():
            structure_ids.add(int(part))
    
    return structure_ids


def _fragment_key(row: Mapping[str, object]) -> Tuple[str, str, str]:
    fragment = str(row.get('Fragment', '')).strip()
    composition = str(row.get('Composition', '')).strip()
    frag_type = str(row.get('Fragment Type', '')).strip()
    return fragment, composition, frag_type


def score_structures_by_uniqueness(
    clean_theoretical_df: pd.DataFrame,
    matched_df: pd.DataFrame,
) -> Dict[str, Dict[str, float]]:
    """
    Score structures using weighted fragment exclusivity.

    Algorithm:
    1. Count fragments per structure (all fragments, not just exclusive)
    2. Weight matched fragments by exclusivity:
       - Exclusive to 1 structure: weight = 1.0
       - Shared by 2 structures: weight = 0.5
       - Shared by N structures: weight = 1/N
    3. Calculate weighted score for ranking
    4. Handle edge cases where no exclusive fragments exist

    Returns:
        Dict mapping structure_id to scores:
        - unique_total: Total theoretical fragments for this structure
        - matched_unique: Total matched fragments for this structure (unweighted count)
        - match_percent: (weighted_matched / unique_total) * 100
        - match_share: (weighted_matched / sum_all_weighted) * 100
        - weighted_score: Weighted sum of matched fragments (for ranking)
    """
    if clean_theoretical_df.empty or matched_df.empty:
        return {}

    # Collect all unique structure IDs from the Structure column
    # Structure column may contain comma-separated values like "1,2,3"
    all_structure_ids = set()
    for _, row in clean_theoretical_df.iterrows():
        structure_ids_in_row = _parse_structure_ids(row.get('Structure'))
        all_structure_ids.update(structure_ids_in_row)
    
    structure_ids = sorted([str(sid) for sid in all_structure_ids])

    if not structure_ids:
        return {}

    # Count all fragments per structure (theoretical)
    unique_by_structure: Dict[str, set] = {sid: set() for sid in structure_ids}
    for _, row in clean_theoretical_df.iterrows():
        structure_ids_in_row = _parse_structure_ids(row.get('Structure'))
        frag_key = _fragment_key(row)
        # Add this fragment to all structures it belongs to
        for sid in structure_ids_in_row:
            unique_by_structure[str(sid)].add(frag_key)

    # Count matched fragments per structure and track fragment exclusivity
    matched_unique_by_structure: Dict[str, set] = {sid: set() for sid in structure_ids}
    fragment_to_structures: Dict[Tuple[str, str, str], set] = {}  # fragment_key -> set of structure_ids
    
    for _, row in matched_df.iterrows():
        structure_ids_in_row = _parse_structure_ids(row.get('Structure'))
        if not structure_ids_in_row:
            continue
        frag_key = _fragment_key(row)
        
        # Add this fragment to all structures it belongs to
        for sid in structure_ids_in_row:
            matched_unique_by_structure[str(sid)].add(frag_key)
        
        # Track which structures share this fragment
        if frag_key not in fragment_to_structures:
            fragment_to_structures[frag_key] = set()
        for sid in structure_ids_in_row:
            fragment_to_structures[frag_key].add(str(sid))

    # Calculate weighted scores based on fragment exclusivity
    weighted_matched_by_structure: Dict[str, float] = {}
    
    for sid in structure_ids:
        weighted_sum = 0.0
        for frag_key in matched_unique_by_structure[sid]:
            # Weight = 1 / (number of structures sharing this fragment)
            num_structures_sharing = len(fragment_to_structures.get(frag_key, {sid}))
            weight = 1.0 / num_structures_sharing
            weighted_sum += weight
        
        weighted_matched_by_structure[sid] = weighted_sum

    # Total weighted matched across all structures
    total_weighted_matched = sum(weighted_matched_by_structure.values())
    if total_weighted_matched <= 0:
        total_weighted_matched = 1.0  # Avoid division by zero

    # Calculate final scores
    scores: Dict[str, Dict[str, float]] = {}
    for sid in structure_ids:
        unique_total = len(unique_by_structure[sid])
        matched_unique = len(matched_unique_by_structure[sid])
        weighted_matched = weighted_matched_by_structure[sid]
        
        # Match % based on weighted score
        match_percent = (weighted_matched / unique_total * 100.0) if unique_total else 0.0
        
        # Match Share % based on contribution to total weighted matches
        match_share = (weighted_matched / total_weighted_matched * 100.0) if total_weighted_matched > 0 else 0.0
        
        scores[sid] = {
            'unique_total': float(unique_total),
            'matched_unique': float(matched_unique),
            'match_percent': round(match_percent, 1),
            'match_share': round(match_share, 1),
            'weighted_score': round(weighted_matched, 2),  # For debugging/ranking
        }

    return scores


def create_structure_prediction_sheet(
    clean_theoretical_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    structure_scores: Dict[str, Dict[str, float]],
) -> pd.DataFrame:
    """
    Build structure prediction report based on weighted fragment exclusivity scoring.
    
    Scoring Algorithm:
    - Fragments exclusive to one structure get weight = 1.0
    - Fragments shared by N structures get weight = 1/N
    - Structures with more exclusive fragments rank higher
    
    Confidence Level: (Weighted Match % / Sum of All Match %) * 10
    Range: 0-10, where higher = more confident prediction
    
    Example:
        Structure 1: 5 exclusive fragments + 3 shared (1/2) = 5.0 + 1.5 = 6.5
        Structure 2: 2 exclusive fragments + 3 shared (1/2) + 4 shared (1/4) = 2.0 + 1.5 + 1.0 = 4.5
        → Structure 1 ranks higher (better exclusivity)
    """
    if not structure_scores:
        return pd.DataFrame()

    rows: List[Dict[str, object]] = []
    
    # Sort by weighted score (descending) for ranking
    ordered = sorted(structure_scores.keys(), 
                    key=lambda s: (-structure_scores[s].get('weighted_score', 0), 
                                  -structure_scores[s]['match_share']))
    
    # Calculate sum of all weighted match percentages for relative confidence
    total_match_pct = sum(score['match_percent'] for score in structure_scores.values())
    
    for rank, sid in enumerate(ordered, start=1):
        score = structure_scores[sid]
        match_pct = score['match_percent']
        weighted_score = score.get('weighted_score', 0.0)
        
        # Confidence level relative to all structures: (Match % / Sum) * 10
        if total_match_pct > 0:
            confidence_level = (match_pct / total_match_pct) * 10
            confidence_str = f"{confidence_level:.1f}"
        else:
            confidence_str = "0.0"
        
        rows.append({
            'Structure': sid,
            'Unique Theoretical': int(score['unique_total']),
            'Matched Unique': int(score['matched_unique']),
            'Weighted Score': round(weighted_score, 2),
            'Match % (Weighted/Total)': match_pct,
            'Match Share (%)': score['match_share'],
            'Rank': rank,
            'Confidence Level (0-10)': confidence_str,
        })

    return pd.DataFrame(rows)


def get_best_structure(
    prediction_df: Optional[pd.DataFrame] = None,
) -> Optional[str]:
    if prediction_df is None or prediction_df.empty:
        return None
    top = prediction_df.sort_values('Rank').head(1)
    if top.empty:
        return None
    return str(top.iloc[0]['Structure']).strip()


def match_fragments_from_table(
    theoretical_table: pd.DataFrame,
    experimental_spectrum: Tuple[np.ndarray, np.ndarray],
    exp_rt_map: Optional[Dict[float, float]] = None,
    exp_area_map: Optional[Dict[float, float]] = None,
    ppm_tolerance: float = 20.0,
    charge_states: List[int] = None,
    intensity_threshold: float = 0.01
) -> pd.DataFrame:
    """
    Match theoretical fragments (from DataFrame) with experimental peaks across multiple charge states.
    
    Args:
        theoretical_table: DataFrame with Fragment, Mass(Da), and other columns from generate_fragment_table
        experimental_spectrum: (mz_array, intensity_array) from MS/MS
        exp_rt_map: Dict mapping experimental m/z to retention time
        exp_area_map: Dict mapping experimental m/z to peak area
        ppm_tolerance: PPM tolerance for matching
        charge_states: List of charge states to try (default [1, 2, 3])
        intensity_threshold: Intensity threshold (as fraction of max)
    
    Returns:
        DataFrame with matched fragments
    """
    if charge_states is None:
        charge_states = [1, 2, 3]

    calc = GlycanMassCalculator()
    
    mz_array, intensity_array = experimental_spectrum
    
    # Normalize intensities to 0-100%
    max_intensity = np.max(intensity_array) if len(intensity_array) > 0 else 1
    normalized_intensity = (intensity_array / max_intensity) * 100
    
    matched_list = []
    
    # Iterate through theoretical fragments
    for _, row in theoretical_table.iterrows():
        neutral_mass = row['Mass(Da)']
        frag_name = row['Fragment']
        frag_type = row['Fragment Type']
        composition = row['Composition']
        
        # Try each charge state
        for z in charge_states:
            theoretical_mz = calc.calculate_mz(neutral_mass, z)
            
            # Find matching experimental peaks
            for i, exp_mz in enumerate(mz_array):
                ppm_error = abs(theoretical_mz - exp_mz) / theoretical_mz * 1e6
                
                if abs(ppm_error) <= ppm_tolerance:
                    intensity = normalized_intensity[i]
                    
                    if intensity >= intensity_threshold:
                        # Get RT for this experimental m/z
                        rt = exp_rt_map.get(exp_mz, 0.0) if exp_rt_map else 0.0
                        area = exp_area_map.get(exp_mz, None) if exp_area_map else None

                        matched_list.append({
                            'Fragment': frag_name,
                            'Fragment Type': frag_type,
                            'Composition': composition,
                            'Structure': row.get('Structure', 'unknown'),
                            'Charge': z,
                            'Experimental m/z': round(exp_mz, 4),
                            'Theoretical m/z': round(theoretical_mz, 4),
                            'PPM Error': round(ppm_error, 2),
                            'Intensity (%)': round(intensity, 2),
                            'Raw Intensity': intensity_array[i],
                            'RT (min)': rt,
                            'Area': float(area) if area is not None else None,
                        })
    
    # Convert to DataFrame and sort by intensity (descending)
    matched_df = pd.DataFrame(matched_list)
    if len(matched_df) > 0:
        matched_df = matched_df.sort_values('Intensity (%)', ascending=False).reset_index(drop=True)
    
    return matched_df
