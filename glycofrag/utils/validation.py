"""
Input validation utilities.

This module provides functions for validating inputs to the glycofrag package,
including glycan codes, peptide sequences, and modifications.
"""

import re
from typing import List, Optional, Tuple

from glycofrag.core.constants import AMINO_ACID_MASSES, MODIFICATION_MASSES


class ValidationError(Exception):
    """
    Exception raised for validation errors.

    Public API:
        - ValidationError
    """
    pass


def validate_glycan_code(code: str) -> Tuple[bool, str]:
    """
    Validate a glycan composition code.
    
    Args:
        code: Glycan code to validate
    
    Returns:
        Tuple of (is_valid, message)
        - is_valid: True if code is valid
        - message: Error message if invalid, empty string if valid
    
    Examples:
        >>> validate_glycan_code("4501")
        (True, '')
        >>> validate_glycan_code("HexNAc(4)Hex(5)")
        (True, '')
        >>> validate_glycan_code("invalid")
        (False, 'Invalid glycan code format: invalid')
    """
    if not code:
        return False, "Glycan code cannot be empty"
    
    code = str(code).strip().replace(" ", "")
    
    # Check numeric format (e.g., "4501", "45010")
    if code.isdigit():
        if len(code) < 4:
            return False, f"Numeric glycan code must have at least 4 digits, got {len(code)}"
        if len(code) > 5:
            return False, f"Numeric glycan code must have at most 5 digits, got {len(code)}"
        return True, ""
    
    # Check named format (e.g., "HexNAc(4)Hex(5)Fuc(1)NeuAc(2)")
    if code[0].isalpha():
        # Pattern to match monosaccharide names with counts
        pattern = r'([A-Za-z]+)\((\d+)\)'
        matches = re.findall(pattern, code)
        
        if not matches:
            return False, f"Invalid glycan code format: {code}"
        
        # Validate monosaccharide names
        valid_names = {
            'HEXNAC', 'HEXN', 'GLCNAC', 'GALNAC',
            'HEX', 'HEXOSE', 'GLC', 'GAL', 'MAN',
            'FUC', 'FUCOSE', 'DEOXYHEX',
            'NEUAC', 'SIA', 'SIALIC',
            'NEUGC',
        }
        
        for mono_name, count in matches:
            if mono_name.upper() not in valid_names:
                return False, f"Unknown monosaccharide: {mono_name}"
            if int(count) < 0:
                return False, f"Monosaccharide count cannot be negative: {mono_name}({count})"
        
        return True, ""
    
    return False, f"Invalid glycan code format: {code}"


def validate_peptide_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate a peptide sequence.
    
    Args:
        sequence: Peptide sequence (single letter codes)
    
    Returns:
        Tuple of (is_valid, message)
        - is_valid: True if sequence is valid
        - message: Error message if invalid, empty string if valid
    
    Examples:
        >>> validate_peptide_sequence("PEPTIDE")
        (True, '')
        >>> validate_peptide_sequence("PEPT1DE")
        (False, 'Invalid character in peptide sequence: 1')
    """
    if not sequence:
        return False, "Peptide sequence cannot be empty"
    
    sequence = str(sequence).strip().upper()
    
    # Check for valid amino acid letters only
    valid_aas = set(AMINO_ACID_MASSES.keys())
    
    for i, char in enumerate(sequence):
        if char not in valid_aas:
            return False, f"Invalid character in peptide sequence at position {i+1}: {char}"
    
    return True, ""


def validate_modification(mod_string: str) -> Tuple[bool, str]:
    """
    Validate a modification string.
    
    Args:
        mod_string: Modification string (e.g., "Ox:M", "CAM:C", "(+15.99):M6")
    
    Returns:
        Tuple of (is_valid, message)
        - is_valid: True if modification is valid
        - message: Error message if invalid, empty string if valid
    
    Examples:
        >>> validate_modification("Ox:M")
        (True, '')
        >>> validate_modification("CAM:C")
        (True, '')
        >>> validate_modification("(+15.99):M6")
        (True, '')
        >>> validate_modification("INVALID")
        (False, 'Modification must be in format "Name:Target": INVALID')
    """
    if not mod_string:
        return False, "Modification string cannot be empty"
    
    mod_string = str(mod_string).strip()
    
    # Check for colon separator
    if ':' not in mod_string:
        return False, f'Modification must be in format "Name:Target": {mod_string}'
    
    parts = mod_string.split(':', 1)
    mod_name = parts[0].strip()
    target = parts[1].strip()
    
    # Validate modification name
    if mod_name.startswith('(') and mod_name.endswith(')'):
        # Custom mass modification
        try:
            mass_str = mod_name[1:-1].strip()
            float(mass_str)
        except ValueError:
            return False, f"Invalid custom mass format: {mod_name}"
    elif mod_name not in MODIFICATION_MASSES:
        # Unknown modification - warn but allow (might be user-defined)
        pass  # Just allow it for flexibility
    
    # Validate target
    target_upper = target.upper()
    if target_upper in ['N-TERM', 'NTERM', 'N', 'C-TERM', 'CTERM', 'C']:
        return True, ""
    
    # Check if target is a single amino acid
    if len(target) == 1 and target.upper() in AMINO_ACID_MASSES:
        return True, ""
    
    # Check if target is amino acid + position (e.g., "M6")
    if len(target) >= 2 and target[0].upper() in AMINO_ACID_MASSES:
        position_part = target[1:]
        if position_part.isdigit():
            return True, ""
    
    # Just a position number
    if target.isdigit():
        return True, ""
    
    return False, f"Invalid modification target: {target}"


def validate_charge_state(charge: int) -> Tuple[bool, str]:
    """
    Validate a charge state.
    
    Args:
        charge: Charge state to validate
    
    Returns:
        Tuple of (is_valid, message)
    """
    if not isinstance(charge, int):
        return False, f"Charge must be an integer, got {type(charge).__name__}"
    
    if charge <= 0:
        return False, f"Charge must be positive, got {charge}"
    
    if charge > 10:
        return False, f"Charge state {charge} is unusually high (max recommended: 10)"
    
    return True, ""


def validate_modification_list(modifications: List[str]) -> Tuple[bool, List[str]]:
    """
    Validate a list of modifications.
    
    Args:
        modifications: List of modification strings
    
    Returns:
        Tuple of (all_valid, list of error messages)
    """
    errors = []
    
    for mod in modifications:
        is_valid, message = validate_modification(mod)
        if not is_valid:
            errors.append(message)
    
    return len(errors) == 0, errors


def sanitize_peptide_sequence(sequence: str) -> str:
    """
    Sanitize a peptide sequence by removing invalid characters.
    
    Args:
        sequence: Raw peptide sequence
    
    Returns:
        Sanitized peptide sequence (uppercase, only valid AAs)
    """
    if not sequence:
        return ""
    
    valid_aas = set(AMINO_ACID_MASSES.keys())
    sanitized = ''.join(
        char for char in str(sequence).upper()
        if char in valid_aas
    )
    
    return sanitized
