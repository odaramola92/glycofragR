# glycofrag/composition/parser.py
"""Glycan composition parsing and validation."""

from typing import Tuple
from glycofrag.core.mass_calculator import GlycanMassCalculator


class GlycanComposition:
    """
    Represents a glycan composition (monosaccharide counts).
    
    This class is responsible ONLY for parsing and storing composition data.
    It does NOT handle:
    - Structure prediction
    - Core calculations
    - Fragmentation
    - Mass calculations (beyond parsing)
    
    Attributes:
        hexnac: Number of HexNAc residues
        hex: Number of Hex residues
        fuc: Number of Fuc residues
        neuac: Number of NeuAc residues
        neugc: Number of NeuGc residues

    Public API:
        - from_code(...)
        - to_code()
        - total_residues()
        - is_valid()
    
    Example:
        >>> # Parse from code
        >>> comp = GlycanComposition.from_code('4501')
        >>> print(comp)
        GlycanComposition(HexNAc=4, Hex=5, Fuc=0, NeuAc=1)
        
        >>> # Direct instantiation
        >>> comp = GlycanComposition(hexnac=4, hex=5, fuc=0, neuac=1)
        >>> print(comp.to_code())
        4501
    """
    
    def __init__(self, hexnac: int, hex: int, fuc: int, neuac: int, neugc: int = 0):
        """
        Initialize a glycan composition.
        
        Args:
            hexnac: Number of HexNAc residues
            hex: Number of Hex residues
            fuc: Number of Fuc residues
            neuac: Number of NeuAc residues
            neugc: Number of NeuGc residues (default: 0)
        """
        self.hexnac = hexnac
        self.hex = hex
        self.fuc = fuc
        self.neuac = neuac
        self.neugc = neugc
    
    @classmethod
    def from_code(cls, code: str) -> 'GlycanComposition':
        """
        Parse glycan composition from 4-digit code.
        
        Args:
            code: 4-digit composition code (e.g., '4501' for HexNAc4Hex5Fuc0NeuAc1)
        
        Returns:
            GlycanComposition object
        
        Example:
            >>> comp = GlycanComposition.from_code('4501')
            >>> comp.hexnac
            4
            >>> comp.hex
            5
        """
        calculator = GlycanMassCalculator()
        hexnac, hex, fuc, neuac, neugc = calculator.parse_glycan_code(code)
        return cls(hexnac, hex, fuc, neuac, neugc)
    
    def __repr__(self):
        """String representation of composition."""
        parts = [f"HexNAc={self.hexnac}", f"Hex={self.hex}"]
        if self.fuc > 0:
            parts.append(f"Fuc={self.fuc}")
        if self.neuac > 0:
            parts.append(f"NeuAc={self.neuac}")
        if self.neugc > 0:
            parts.append(f"NeuGc={self.neugc}")
        return f"GlycanComposition({', '.join(parts)})"
    
    def to_code(self) -> str:
        """
        Convert composition back to 4-digit code.
        
        Returns:
            4-digit code string
        
        Example:
            >>> comp = GlycanComposition(hexnac=4, hex=5, fuc=0, neuac=1)
            >>> comp.to_code()
            '4501'
        """
        return f"{self.hexnac}{self.hex}{self.fuc}{self.neuac}"
    
    def total_residues(self) -> int:
        """Return total number of monosaccharides."""
        return self.hexnac + self.hex + self.fuc + self.neuac + self.neugc
    
    def is_valid(self) -> bool:
        """
        Check if composition is valid (all counts non-negative).
        
        Returns:
            True if valid, False otherwise
        """
        return all(count >= 0 for count in [self.hexnac, self.hex, self.fuc, self.neuac, self.neugc])
