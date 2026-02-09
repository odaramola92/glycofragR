"""
Glycopeptide module for combined peptide and glycan fragmentation analysis.

This module provides the Glycopeptide class for analyzing glycosylated peptides,
including Y0 (peptide-only), Y1/Y2/... (peptide+glycan), and intact ions.
"""

from .glycopeptide import Glycopeptide

__all__ = ['Glycopeptide']
