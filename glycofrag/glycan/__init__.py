"""
Glycan module for glycan structure prediction and fragmentation analysis.

This module provides classes and functions for:
- Building N-glycan and O-glycan core structures
- Predicting possible glycan structures from composition
- Generating fragment ions (B, Y, BY, YY ions)
- Classifying glycan structures

NOTE: This module now re-exports from the new modular architecture.
The monolithic glycan.py has been removed. Use: from glycofrag import Glycan
"""

from glycofrag.glycan.facade import Glycan
from glycofrag.core.modifications import ReducingEndType

__all__ = ['Glycan', 'ReducingEndType']
