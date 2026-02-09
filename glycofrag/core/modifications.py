"""
Peptide modification parsing and handling.

This module provides functions for parsing modification strings
and applying modifications to peptide sequences.

Modification String Formats:
    - "ModName:Position" - e.g., "Ox:M6" (oxidation at methionine position 6)
    - "ModName:AA" - e.g., "CAM:C" (carbamidomethylation on all cysteines)
    - "AA:ModName" - e.g., "M:Ox" (oxidation on all methionines)
    - "AA1,5:ModName" - e.g., "M4,24:Ox" (targeted positions)
    - "Custom:AA, mass" - e.g., "Custom:A, 45.8978" (custom mass on all A)
    - "(+mass):Position" - e.g., "(+15.99):M6" (custom mass at position)
    - Multiple mods separated by semicolon: "Ox:M6;Phos:S3"
"""

import re
from typing import List, Tuple, Optional, Union
from enum import IntEnum


def parse_modification_string(mod_string: str) -> List[Tuple[str, Union[str, int]]]:
    """
    Parse a modification string into a list of (modification, position) tuples.
    
    Args:
        mod_string: Modification string in various formats:
            - "ModName:Position" (e.g., "Ox:M6")
            - "ModName:AA" (e.g., "CAM:C")
            - "AA:ModName" (e.g., "M:Ox")
            - "AA1,5:ModName" (e.g., "M4,24:Ox")
            - "Custom:AA, mass" (e.g., "Custom:A, 45.8978")
            - "(+mass):Position" (e.g., "(+15.99):M6")
            - Multiple separated by semicolons
    
    Returns:
        List of tuples (modification_name, position) where position can be:
            - int: 1-indexed position number
            - str: amino acid letter, "N-term", or "C-term"
    
    Examples:
        >>> parse_modification_string("Ox:M6")
        [('Ox', 'M6')]
        >>> parse_modification_string("M:Ox")
        [('Ox', 'M')]
        >>> parse_modification_string("M4,24:Ox")
        [('Ox', 'M4'), ('Ox', 'M24')]
        >>> parse_modification_string("Custom:A, 45.8978")
        [('(+45.8978)', 'A')]
        >>> parse_modification_string("CAM:C")
        [('CAM', 'C')]
        >>> parse_modification_string("Ox:M6;Phos:S3")
        [('Ox', 'M6'), ('Phos', 'S3')]
        >>> parse_modification_string("(+15.99):M6")
        [('(+15.99)', 'M6')]
    """
    if not mod_string or (hasattr(mod_string, '__class__') and 
                          mod_string.__class__.__name__ == 'NAType'):
        return []
    
    # Handle pandas NA
    try:
        import pandas as pd
        if pd.isna(mod_string):
            return []
    except ImportError:
        pass
    
    mod_string = str(mod_string).strip()
    if not mod_string:
        return []
    
    result = []
    
    def _looks_like_target(token: str) -> bool:
        token = token.strip()
        if not token:
            return False
        token_lower = token.lower()
        if token_lower in ['n-term', 'nterm', 'n', 'c-term', 'cterm', 'c']:
            return True
        if token.isdigit():
            return True
        if re.fullmatch(r"[A-Za-z](\d+(,\d+)*)?", token):
            return True
        return False

    def _expand_target(target: str) -> List[str]:
        target = target.strip()
        if not target:
            return []
        target_lower = target.lower()
        if target_lower in ['n-term', 'nterm', 'n']:
            return ['N-term']
        if target_lower in ['c-term', 'cterm', 'c']:
            return ['C-term']
        if target.isdigit():
            return [target]
        if re.fullmatch(r"[A-Za-z](\d+(,\d+)*)", target):
            aa = target[0].upper()
            numbers = [num.strip() for num in target[1:].split(',') if num.strip()]
            return [f"{aa}{num}" for num in numbers]
        if len(target) == 1 and target.isalpha():
            return [target.upper()]
        return [target]

    def _parse_custom_target_mass(value: str) -> Tuple[Optional[str], Optional[float]]:
        if not value:
            return None, None
        if ',' in value:
            parts = [part.strip() for part in value.split(',') if part.strip()]
        else:
            parts = [part for part in value.split() if part.strip()]
        if len(parts) < 2:
            return None, None
        target = parts[0]
        try:
            mass = float(parts[1])
        except ValueError:
            return target, None
        return target, mass

    # Split multiple modifications (separated by semicolons)
    for mod_part in mod_string.split(';'):
        mod_part = mod_part.strip()
        if not mod_part:
            continue
        
        if ':' in mod_part:
            parts = mod_part.split(':', 1)
            left = parts[0].strip()
            right = parts[1].strip()

            if left.lower() == 'custom':
                target, mass = _parse_custom_target_mass(right)
                if target and mass is not None:
                    mod_name = f"({mass:+g})"
                    for position in _expand_target(target):
                        result.append((mod_name, position))
                continue

            if _looks_like_target(left) and not _looks_like_target(right):
                target = left
                mod_name = right
            else:
                mod_name = left
                target = right

            for position in _expand_target(target):
                result.append((mod_name, position))
    
    return result

def parse_custom_mass_mod(mod_string: str) -> Optional[float]:
    """
    Parse a custom mass modification in the format "(+15.01)" or "(-2.1)".
    
    Args:
        mod_string: String that may contain a custom mass modification
            in parentheses with a + or - sign.
    
    Returns:
        Float value of the modification mass, or None if not a custom mass.
    
    Examples:
        >>> parse_custom_mass_mod("(+15.99)")
        15.99
        >>> parse_custom_mass_mod("(-17.03)")
        -17.03
        >>> parse_custom_mass_mod("Ox")
        None
    """
    if not mod_string:
        return None
    
    mod_string = str(mod_string).strip()
    
    # Check if the string is enclosed in parentheses and contains a sign
    if (mod_string.startswith('(') and mod_string.endswith(')') and 
        ('+' in mod_string or '-' in mod_string)):
        try:
            # Remove parentheses and convert to float
            mass_str = mod_string[1:-1].strip()
            return float(mass_str)
        except ValueError:
            return None
    
    return None

def apply_modifications_to_positions(
    peptide: str,
    mod_string: Optional[str] = None,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    use_cam: bool = True,
) -> dict:
    """
    Determine which positions in a peptide have modifications applied.
    
    IMPORTANT: This function ONLY applies modifications when the target amino acid
    is present in the peptide. For example, CAM is ONLY applied to positions with 'C',
    and Ox is ONLY applied to positions with 'M'.
    
    Args:
        peptide: Peptide sequence (single letter amino acid codes)
        mod_string: Modification string (e.g., "Ox:M6;Phos:S3")
        fixed_mods: List of fixed modifications (e.g., ["CAM:C"])
        variable_mods: List of variable modifications (e.g., ["Ox:M"])
        use_cam: Whether to automatically apply CAM to cysteines (ONLY if 'C' present)
    
    Returns:
        Dictionary mapping positions to modification names:
        {
            position_index: modification_name,
            'N-term': modification_name,  # if N-terminal mod
            'C-term': modification_name,  # if C-terminal mod
        }
    
    Example:
        >>> apply_modifications_to_positions("LCPDCPLLAPLNDSR", 
        ...                                   fixed_mods=["CAM:C"])
        {1: 'CAM', 4: 'CAM'}
        >>> apply_modifications_to_positions("LPDPLLAPLNDSR",  # No C
        ...                                   fixed_mods=["CAM:C"])
        {}  # No CAM applied - no cysteines
    """
    from glycofrag.core.constants import MODIFICATION_TARGETS
    
    modifications = {}
    fixed_mods = fixed_mods or []
    variable_mods = variable_mods or []
    
    # Priority 1: Parse mod_string (Excel/external modifications)
    if mod_string:
        parsed_mods = parse_modification_string(mod_string)
        for mod_name, position in parsed_mods:
            custom_mass = parse_custom_mass_mod(mod_name)
            
            if isinstance(position, str):
                if position in ['N-term', 'C-term']:
                    if custom_mass is not None:
                        modifications[position] = ('CUSTOM', custom_mass)
                    else:
                        modifications[position] = mod_name
                elif len(position) >= 2 and position[0].isalpha() and position[1:].isdigit():
                    # Format like "M6" - amino acid at specific position
                    amino_acid = position[0].upper()
                    pos = int(position[1:]) - 1  # Convert to 0-indexed
                    
                    # Verify the amino acid matches the sequence
                    if 0 <= pos < len(peptide) and peptide[pos].upper() == amino_acid:
                        if custom_mass is not None:
                            modifications[pos] = ('CUSTOM', custom_mass)
                        else:
                            modifications[pos] = mod_name
                elif len(position) == 1 and position.isalpha():
                    # Single amino acid - apply to all occurrences ONLY if present
                    target_aa = position.upper()
                    # CRITICAL: Only apply modification if the target AA exists in the sequence
                    if target_aa in peptide.upper():
                        for i, aa in enumerate(peptide.upper()):
                            if aa == target_aa and i not in modifications:
                                if custom_mass is not None:
                                    modifications[i] = ('CUSTOM', custom_mass)
                                else:
                                    modifications[i] = mod_name
    
    # Priority 2: Fixed modifications
    # Handle CAM specially - ONLY apply if use_cam is True AND 'C' is in the sequence
    if use_cam and 'C' in peptide.upper():
        # CRITICAL: Only apply CAM if cysteine is actually present
        for i, aa in enumerate(peptide.upper()):
            if aa == 'C' and i not in modifications:
                modifications[i] = 'CAM'
    
    # Process other fixed modifications
    for mod in fixed_mods:
        if mod in ['CAM:C', 'CAM']:
            # Skip - already handled above
            continue
        
        if ':' in mod:
            mod_type, target = mod.split(':', 1)
            if len(target) == 1 and target.isalpha():
                # Apply to all occurrences - ONLY if the target AA is in the sequence
                target_upper = target.upper()
                if target_upper in peptide.upper():  # CRITICAL check
                    for i, aa in enumerate(peptide.upper()):
                        if aa == target_upper and i not in modifications:
                            modifications[i] = mod_type
            elif target in ['N-term', 'C-term']:
                if target not in modifications:
                    modifications[target] = mod_type
        else:
            # Modification without explicit target - use default targets
            if mod in MODIFICATION_TARGETS:
                for target in MODIFICATION_TARGETS[mod]:
                    if len(target) == 1:
                        # CRITICAL: Only apply if target AA is in sequence
                        if target in peptide.upper():
                            for i, aa in enumerate(peptide.upper()):
                                if aa == target and i not in modifications:
                                    modifications[i] = mod
    
    # Priority 3: Variable modifications
    for mod in variable_mods:
        if ':' in mod:
            mod_type, target = mod.split(':', 1)
            if len(target) == 1 and target.isalpha():
                # CRITICAL: Only apply if target AA is in sequence
                target_upper = target.upper()
                if target_upper in peptide.upper():
                    for i, aa in enumerate(peptide.upper()):
                        if aa == target_upper and i not in modifications:
                            modifications[i] = mod_type
            elif target in ['N-term', 'C-term']:
                if target not in modifications:
                    modifications[target] = mod_type
        else:
            if mod in MODIFICATION_TARGETS:
                for target in MODIFICATION_TARGETS[mod]:
                    if len(target) == 1:
                        # CRITICAL: Only apply if target AA is in sequence
                        if target in peptide.upper():
                            for i, aa in enumerate(peptide.upper()):
                                if aa == target and i not in modifications:
                                    modifications[i] = mod
    
    return modifications

def format_modification_summary(modifications: dict, peptide: str) -> str:
    """
    Format a human-readable summary of modifications applied to a peptide.
    
    Args:
        modifications: Dictionary of position -> modification name
        peptide: Peptide sequence
    
    Returns:
        Formatted string describing modifications
    
    Example:
        >>> mods = {1: 'CAM', 4: 'CAM', 5: 'Ox'}
        >>> format_modification_summary(mods, "LCPDCMLLAPLNDSR")
        'CAM(C2,C5), Ox(M6)'
    """
    if not modifications:
        return "No modifications"
    
    # Group by modification type
    mod_groups = {}
    for pos, mod_name in modifications.items():
        if isinstance(mod_name, tuple) and mod_name[0] == 'CUSTOM':
            mod_key = f"Custom({mod_name[1]:+.2f})"
        else:
            mod_key = mod_name
        
        if mod_key not in mod_groups:
            mod_groups[mod_key] = []
        
        if isinstance(pos, int):
            aa = peptide[pos] if pos < len(peptide) else '?'
            mod_groups[mod_key].append(f"{aa}{pos+1}")
        else:
            mod_groups[mod_key].append(pos)
    
    # Format output
    parts = []
    for mod_name, positions in sorted(mod_groups.items()):
        positions_str = ','.join(positions)
        parts.append(f"{mod_name}({positions_str})")
    
    return ', '.join(parts)


def create_modification_detail_report(
    peptide: str,
    modifications: dict,
    use_cam: bool = False,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    mod_string: Optional[str] = None,
) -> List[dict]:
    """
    Create a detailed report of all modifications applied to a peptide.
    
    This function generates a list of dictionaries detailing each modification,
    including:
    - Modification type
    - Target amino acid(s)
    - Positions where applied
    - Whether the position actually has the target amino acid
    - Mass delta
    
    Args:
        peptide: Peptide sequence
        modifications: Dictionary of position -> modification name
        use_cam: Whether CAM was enabled
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Modification string
    
    Returns:
        List of dictionaries with modification details
    """
    from glycofrag.core.constants import MODIFICATION_MASSES, MODIFICATION_TARGETS
    
    report = []
    
    # Track which modifications have been reported
    reported_mods = set()
    
    # Process actual applied modifications
    for pos, mod_name in modifications.items():
        if isinstance(pos, int):
            aa = peptide[pos].upper()
            mass_delta = 0.0
            
            if isinstance(mod_name, tuple) and mod_name[0] == 'CUSTOM':
                mass_delta = mod_name[1]
                mod_display = f"Custom ({mod_name[1]:+.6f} Da)"
            else:
                mass_delta = MODIFICATION_MASSES.get(mod_name, 0.0)
                mod_display = mod_name
            
            report.append({
                "Modification": mod_display,
                "Type": "Applied",
                "Position": f"{pos+1}",
                "Amino Acid": aa,
                "Mass Delta (Da)": f"{mass_delta:.6f}",
            })
            reported_mods.add(mod_name)
    
    # Also report terminal modifications
    for term in ['N-term', 'C-term']:
        if term in modifications:
            mod_name = modifications[term]
            mass_delta = 0.0
            
            if isinstance(mod_name, tuple) and mod_name[0] == 'CUSTOM':
                mass_delta = mod_name[1]
                mod_display = f"Custom ({mod_name[1]:+.6f} Da)"
            else:
                mass_delta = MODIFICATION_MASSES.get(mod_name, 0.0)
                mod_display = mod_name
            
            report.append({
                "Modification": mod_display,
                "Type": "Applied",
                "Position": term,
                "Amino Acid": "N/A",
                "Mass Delta (Da)": f"{mass_delta:.6f}",
            })
            reported_mods.add(mod_name)
    
    # Report CAM status
    if use_cam:
        if 'C' in peptide.upper():
            c_count = peptide.upper().count('C')
            cam_mass = MODIFICATION_MASSES.get('CAM', 0.0)
            report.append({
                "Modification": "CAM (Carbamidomethyl)",
                "Type": "Fixed (Default)",
                "Position": "All C",
                "Amino Acid": "C",
                "Mass Delta (Da)": f"{cam_mass:.6f} per C (x{c_count})",
            })
        else:
            report.append({
                "Modification": "CAM (Carbamidomethyl)",
                "Type": "Fixed (Default) - NOT APPLIED",
                "Position": "N/A",
                "Amino Acid": "C (absent)",
                "Mass Delta (Da)": "0.0",
            })
    
    # Report unmatched fixed modifications
    for mod in (fixed_mods or []):
        if mod not in reported_mods and mod not in ['CAM:C', 'CAM']:
            report.append({
                "Modification": mod,
                "Type": "Fixed (Not Applied)",
                "Position": "N/A",
                "Amino Acid": "Not found",
                "Mass Delta (Da)": "0.0",
            })
    
    return report


# =============================================================================
# Glycan Reducing End Modification Types
# =============================================================================

class ReducingEndType(IntEnum):
    """
    Enumeration of glycan reducing end modification types.
    
    These modifications affect the mass calculation of the reducing end HexNAc.
    
    Usage - you can use either the enum directly or the string API:
        >>> modification_type = get_modification_type("permethylated_reduced")
        >>> # or
        >>> from glycofrag.core.modifications import ReducingEndType
        >>> modification_type = ReducingEndType.PERMETHYLATED_REDUCED

    Public API:
        - FREE
        - REDUCED
        - PERMETHYLATED_FREE
        - PERMETHYLATED_REDUCED
        - TWO_AB
        - TWO_AB_PERMETHYLATED
        - GLYCOPEPTIDE
    """
    FREE = 0                      # Free reducing end
    REDUCED = 1                   # Reduced (alditol)
    PERMETHYLATED_FREE = 2        # Permethylated free reducing end
    PERMETHYLATED_REDUCED = 3     # Permethylated reduced
    TWO_AB = 4                    # 2-AB labeled
    TWO_AB_PERMETHYLATED = 5      # 2-AB labeled + permethylated
    GLYCOPEPTIDE = 6              # Attached to peptide (internal — use Glycopeptide class)
    CUSTOM = 7                    # User-defined custom mass modifications

# Mapping for string-based modification lookup (user-friendly API)
_MODIFICATION_TYPE_MAP = {
    'free': ReducingEndType.FREE,
    'reduced': ReducingEndType.REDUCED,
    'permethylated_free': ReducingEndType.PERMETHYLATED_FREE,
    'permethylated_reduced': ReducingEndType.PERMETHYLATED_REDUCED,
    'two_ab': ReducingEndType.TWO_AB,
    '2ab': ReducingEndType.TWO_AB,
    'two_ab_permethylated': ReducingEndType.TWO_AB_PERMETHYLATED,
    '2ab_permethylated': ReducingEndType.TWO_AB_PERMETHYLATED,
    'glycopeptide': ReducingEndType.GLYCOPEPTIDE,
    'pep': ReducingEndType.GLYCOPEPTIDE,
    'custom': ReducingEndType.CUSTOM,
    # Aliases for common patterns
    'pm_reduced': ReducingEndType.PERMETHYLATED_REDUCED,
    'pm_free': ReducingEndType.PERMETHYLATED_FREE,
}

def get_modification_type(modification: Union[str, int]) -> int:
    """
    Convert a modification type to its integer enum value.
    
    Provides a user-friendly API for specifying glycan modifications.
    
    Args:
        modification: Either a string name or an integer enum value
            String options (case-insensitive):
            - 'free': Free reducing end
            - 'reduced': Reduced (alditol) end
            - 'permethylated_free' or 'pm_free': Permethylated free end
            - 'permethylated_reduced' or 'pm_reduced': Permethylated reduced end
            - 'two_ab' or '2ab': 2-AB labeled
            - 'two_ab_permethylated' or '2ab_permethylated': 2-AB + permethylated
            - 'glycopeptide' or 'pep': Attached to peptide
            
            Or pass an integer 0-6 directly.
    
    Returns:
        Integer modification type (0-6)
    
    Raises:
        ValueError: If modification string is not recognized
    
    Examples:
        >>> get_modification_type('permethylated_reduced')
        3
        >>> get_modification_type('pm_reduced')
        3
        >>> get_modification_type(3)  # Pass through integer
        3
        >>> get_modification_type('FREE')  # Case-insensitive
        0
    """
    if isinstance(modification, int):
        if 0 <= modification <= 7:
            return modification
        raise ValueError(f"Modification type must be 0-7, got {modification}")
    
    if isinstance(modification, str):
        mod_lower = modification.lower().strip()
        if mod_lower in _MODIFICATION_TYPE_MAP:
            return int(_MODIFICATION_TYPE_MAP[mod_lower])
        raise ValueError(
            f"Unknown modification type: '{modification}'\n"
            f"Valid options: {', '.join(sorted(_MODIFICATION_TYPE_MAP.keys()))}"
        )
    
    raise TypeError(f"Modification must be str or int, got {type(modification)}")

def get_modification_name(modification: Union[str, int]) -> str:
    """
    Get human-readable name for a modification type.
    
    Args:
        modification: Either a string name or an integer enum value
    
    Returns:
        Human-readable modification name
    
    Examples:
        >>> get_modification_name(0)
        'Free Reducing End'
        >>> get_modification_name('pm_reduced')
        'Permethylated Reduced'
        >>> get_modification_name(3)
        'Permethylated Reduced'
    """
    # Convert to int if string
    if isinstance(modification, str):
        mod_type = get_modification_type(modification)
    else:
        mod_type = modification
    
    # Map to human-readable names
    names = {
        0: 'Free Reducing End',
        1: 'Reduced (Alditol)',
        2: 'Permethylated Free',
        3: 'Permethylated Reduced',
        4: '2-AB Labeled',
        5: '2-AB + Permethylated',
        6: 'Glycopeptide',
        7: 'Custom',
    }
    
    return names.get(mod_type, f'Unknown ({mod_type})')
