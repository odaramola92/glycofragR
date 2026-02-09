"""
Cache key generation utilities.

This module provides functions for creating consistent cache keys
used throughout the glycofrag package for performance optimization.
"""

from typing import List, Optional, Dict, Any
from functools import lru_cache
import pickle
import hashlib


def create_peptide_cache_key(
    peptide: str,
    use_cam: bool = True,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    mod_string: Optional[str] = None,
) -> str:
    """
    Create a unique cache key for a peptide with its modifications.
    
    This function generates a consistent string key that uniquely identifies
    a peptide sequence with all its modification parameters. This is used
    for caching peptide mass calculations.
    
    Args:
        peptide: Peptide sequence
        use_cam: Whether carbamidomethylation is applied to cysteines
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Modification string from external source
    
    Returns:
        A unique string key for this peptide + modifications combination
    
    Example:
        >>> key = create_peptide_cache_key("PEPTIDE", use_cam=True, fixed_mods=["CAM:C"])
        >>> print(key)
        'PEPTIDE||cam:True||fixed:CAM:C'
    """
    fixed_mods = fixed_mods or []
    variable_mods = variable_mods or []
    
    # Ensure peptide is a string
    peptide = str(peptide or "")
    
    # Sort modifications for consistent keys
    fixed_sorted = sorted([str(mod) for mod in fixed_mods])
    var_sorted = sorted([str(mod) for mod in variable_mods])
    
    # Build key components
    components = [peptide]
    
    # Always include CAM status
    components.append(f"cam:{use_cam}")
    
    if fixed_sorted:
        components.append("fixed:" + ",".join(fixed_sorted))
    
    if var_sorted:
        components.append("var:" + ",".join(var_sorted))
    
    if mod_string:
        components.append("mod:" + str(mod_string))
    
    return "||".join(components)


def create_fragment_cache_key(
    glycan_code: str,
    peptide: Optional[str] = None,
    modification_type: int = 6,
    use_cam: bool = False,
    fixed_mods: Optional[List[str]] = None,
    variable_mods: Optional[List[str]] = None,
    mod_string: Optional[str] = None,
    glycan_type: Optional[str] = None,
) -> str:
    """
    Create a unique cache key for fragment table generation.
    
    This function generates a consistent string key that uniquely identifies
    all parameters affecting fragment generation for a glycan/glycopeptide.
    
    Args:
        glycan_code: Glycan composition code
        peptide: Peptide sequence (if glycopeptide)
        modification_type: Type of reducing end modification (0-6)
        use_cam: Whether carbamidomethylation is enabled
        fixed_mods: List of fixed modifications
        variable_mods: List of variable modifications
        mod_string: Modification string from external source
        glycan_type: Type of glycan ('N' or 'O')
    
    Returns:
        A unique string key for this fragment generation configuration
    
    Example:
        >>> key = create_fragment_cache_key("4501", "PEPTIDE", glycan_type="N")
        >>> print(key)
        '4501||PEPTIDE||mod_type:6||glycan_type:N||cam:False'
    """
    fixed_mods = fixed_mods or []
    variable_mods = variable_mods or []
    
    # Sort modifications for consistent keys
    fixed_sorted = sorted([str(mod) for mod in fixed_mods])
    var_sorted = sorted([str(mod) for mod in variable_mods])
    
    # Build key components
    components = [str(glycan_code), peptide or ""]
    
    # Add modification type
    components.append(f"mod_type:{modification_type}")
    
    # Add glycan type if available
    if glycan_type:
        components.append(f"glycan_type:{glycan_type}")
    
    # Always include CAM status
    components.append(f"cam:{use_cam}")
    
    if fixed_sorted:
        components.append("fixed:" + ",".join(fixed_sorted))
    
    if var_sorted:
        components.append("var:" + ",".join(var_sorted))
    
    if mod_string:
        components.append("mod:" + str(mod_string))
    
    return "||".join(components)


def create_structure_cache_key(
    glycan_code: str,
    glycan_type: str = "N",
    max_structures: int = 100,
    isomer_sensitive: bool = False
) -> str:
    """
    Create a unique cache key for glycan structure prediction.
    
    Args:
        glycan_code: Glycan composition code (e.g., '4501')
        glycan_type: Type of glycan ('N' or 'O')
        max_structures: Maximum number of structures to predict
        isomer_sensitive: If True, treats mirror images as distinct structures
    
    Returns:
        A unique string key for this structure prediction configuration
    """
    return f"{glycan_code}||{glycan_type}||{max_structures}||isomer_sens:{isomer_sensitive}"


# Global structure cache for glycan structure prediction
_STRUCTURE_CACHE: Dict[str, tuple] = {}


def get_cached_structures(cache_key: str) -> Optional[tuple]:
    """Retrieve cached glycan structures if available."""
    return _STRUCTURE_CACHE.get(cache_key)


def cache_structures(cache_key: str, structures: List[Any], fingerprints: set) -> None:
    """Cache predicted glycan structures for future use."""
    if len(_STRUCTURE_CACHE) >= 128:
        first_key = next(iter(_STRUCTURE_CACHE))
        del _STRUCTURE_CACHE[first_key]
    
    import copy
    structures_copy = [copy.deepcopy(s) for s in structures]
    fingerprints_copy = fingerprints.copy()
    
    _STRUCTURE_CACHE[cache_key] = (structures_copy, fingerprints_copy)


def clear_structure_cache() -> None:
    """Clear all cached structures."""
    _STRUCTURE_CACHE.clear()


def get_structure_cache_stats() -> Dict[str, int]:
    """Get statistics about the structure cache."""
    total_structures = sum(len(structures) for structures, _ in _STRUCTURE_CACHE.values())
    return {
        'cache_entries': len(_STRUCTURE_CACHE),
        'total_structures': total_structures,
        'max_cache_size': 128
    }
