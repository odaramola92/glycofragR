# Changelog

All notable changes to glycofrag will be documented in this file.

## [0.1.0] - 2025-07-15

### Added

- **Glycan fragmentation**: B/Y, C/Z, and A-type ion series with automatic oxonium/diagnostic ion generation
- **Glycopeptide fragmentation**: Combined peptide backbone (b/y) and glycan fragmentation with stub ions
- **Neutral losses**: H₂O, NH₃, CO₂, CH₄OS, CO neutral loss generation matching GlypPRM conventions
- **Structure prediction**: Automatic glycan structure enumeration from monosaccharide codes (e.g., "5411")
- **Structure classification**: Categorizes structures as high-mannose, complex, or hybrid
- **Mass calculation**: Monoisotopic mass computation for glycans, peptides, and glycopeptides
- **Fragment tables**: Publication-ready DataFrame output via `generate_fragment_table()`
- **Batch processing**: Multi-glycopeptide analysis with `BatchProcessor`
- **Excel export**: Optional Excel output with `openpyxl` (`pip install glycofrag[excel]`)
- **Visualization**: SNFG-compliant glycan structure drawing with `GlycanVisualizer`
- **Unified API**: Single import path — `from glycofrag import Glycan, Glycopeptide, GlycanVisualizer`
- **Modification support**: Fucosylation, sialylation (NeuAc/NeuGc), bisecting GlcNAc
- **Caching**: LRU-based fragment caching for repeated analyses
- **354 unit tests** covering core, fragmentation, I/O, integration, and edge cases
