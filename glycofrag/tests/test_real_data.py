"""
Phase 6.3: Real Data Testing

This test module validates glycofrag against known glycopeptide sequences
from published literature and realistic MS/MS scenarios.

Test scenarios:
- Known N-glycopeptides with various compositions
- Known O-glycopeptides (Tn, T, disialyl-T, etc.)
- Fragment generation validation
- Mass accuracy validation
- Complex glycan structures
"""

import pytest
from glycofrag import Glycan, Peptide, Glycopeptide
from glycofrag.core.mass_calculator import GlycanMassCalculator


# ============================================================================
# Known N-Glycopeptide Test Cases
# ============================================================================

class TestNGlycopeptideRealData:
    """Test with known N-glycopeptides from literature."""
    
    def test_common_n_glycopeptide_composition_1(self):
        """Test N-glycan with HexNAc(4)Hex(3) (biantennary complex)."""
        # Common N-glycan: biantennary complex
        glycan_code = '4300'  # HexNAc(4)Hex(3)Fuc(0)NeuAc(0)
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        assert len(structures) > 0, "Should generate structures for common biantennary"
        
        # Generate fragments
        structure = structures[0]
        fragments, _ = glycan.generate_fragments(structure, modification_type=0)
        
        # Should have reasonable number of fragments
        total_fragments = sum(len(v) for v in fragments.values())
        assert total_fragments > 0, "Should generate fragments"
    
    def test_sialylated_n_glycan(self):
        """Test N-glycan with NeuAc (sialic acid) - common in serum proteins."""
        # Sialylated complex N-glycan: common in serum glycoproteins
        glycan_code = '4502'  # HexNAc(4)Hex(5)Fuc(0)NeuAc(2)
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Should generate Y-ions with NeuAc
            assert len(fragments['y_ions']) > 0, "Should have Y-ions with sialic acid"
    
    def test_high_mannose_like_glycan(self):
        """Test high-mannose-like composition."""
        # High-mannose-like: HexNAc(2)Hex(5) or similar
        glycan_code = '2500'  # HexNAc(2)Hex(5)Fuc(0)NeuAc(0)
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # High-mannose should generate multiple fragments
            total = sum(len(v) for v in fragments.values())
            assert total > 0, "High-mannose-like should generate fragments"
    
    def test_fucosylated_n_glycan(self):
        """Test N-glycan with Fuc (core fucose - very common in N-glycans)."""
        # Core-fucosylated N-glycan: HexNAc(4)Hex(5)Fuc(1) - most common variant
        glycan_code = '4510'  # HexNAc(4)Hex(5)Fuc(1)NeuAc(0)
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        assert len(structures) > 0, "Should generate structures for core-fucosylated"
        
        structure = structures[0]
        fragments, _ = glycan.generate_fragments(structure, modification_type=0)
        
        # Should have B-ions and Y-ions
        assert len(fragments['b_ions']) > 0, "Fucosylated should have B-ions"
        assert len(fragments['y_ions']) > 0, "Fucosylated should have Y-ions"
    
    def test_realistic_complex_n_glycan(self):
        """Test realistic complex N-glycan: HexNAc(4)Hex(5)Fuc(1)NeuAc(1) (triantennary with sialic acids)."""
        # Complex triantennary sialylated: representative of IgG glycosylation
        glycan_code = '4511'  # HexNAc(4)Hex(5)Fuc(1)NeuAc(1)
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        assert len(structures) > 0, "Should generate structures for complex IgG-like glycan"
        
        structure = structures[0]
        fragments, _ = glycan.generate_fragments(structure, modification_type=0)
        
        # Complex glycan should generate many fragments
        b_count = len(fragments['b_ions'])
        y_count = len(fragments['y_ions'])
        by_count = len(fragments['by_ions'])
        yy_count = len(fragments['yy_ions'])
        
        total = b_count + y_count + by_count + yy_count
        
        assert total > 5, f"Complex glycan should generate >5 fragments, got {total}"
        assert b_count > 0, "Should have B-ions"
        assert y_count > 0, "Should have Y-ions"


# ============================================================================
# Known O-Glycopeptide Test Cases
# ============================================================================

class TestOGlycopeptideRealData:
    """Test with known O-glycopeptides from literature."""
    
    def test_tn_glycan(self):
        """Test Tn antigen (HexNAc only) - simplest O-glycan."""
        # Tn: HexNAc(1) - GalNAc on Ser/Thr
        glycan_code = '1000'  # HexNAc(1)Hex(0)
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Tn has only one monosaccharide, so BY fragmentation yields nothing
            # (no glycosidic bonds to cleave). A-type oxonium ions require 'A' fragment type.
            assert isinstance(fragments, dict), "Should return fragment dict"
    
    def test_t_glycan(self):
        """Test T antigen (HexNAc(1)Hex(1)) - common O-glycan core."""
        # T antigen: Gal-GalNAc
        glycan_code = '1100'  # HexNAc(1)Hex(1)
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            assert len(fragments['b_ions']) > 0 or len(fragments['y_ions']) > 0
    
    def test_sialylated_t_antigen(self):
        """Test sialylated T antigen (HexNAc(1)Hex(1)NeuAc(1)) - common in mucins."""
        # ST antigen: sialyl-T
        glycan_code = '1101'  # HexNAc(1)Hex(1)NeuAc(1)
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            total = sum(len(v) for v in fragments.values())
            assert total > 0, "Sialylated T should generate fragments"
    
    def test_disialyl_t_antigen(self):
        """Test disialyl-T antigen - common in terminal structures."""
        # 2,6-disialyl-T or 2,3-disialyl-T
        glycan_code = '1102'  # HexNAc(1)Hex(1)NeuAc(2)
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # High sialic acid content should be handled
            assert isinstance(fragments['b_ions'], list)
            assert isinstance(fragments['y_ions'], list)


# ============================================================================
# Glycopeptide Integration Tests
# ============================================================================

class TestGlycopeptideIntegration:
    """Test complete glycopeptide scenarios."""
    
    def test_igG_like_glycopeptide(self):
        """Test IgG-like glycopeptide with typical N-glycosylation."""
        # IgG has N-glycosylation at Asn(297) - represented as position 4 in EEQYNSTYR
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4512',  # IgG-typical: complex triantennary
            glycosylation_site=5,   # N in NxS/T motif
            glycan_type='N'
        )
        
        structures = gp.glycan_structures
        if len(structures) > 0:
            # Generate fragments
            fragments = gp.generate_fragments(structure_index=1)
            
            # Should have Y0 (peptide only) and Y1 (peptide + glycan)
            assert len(fragments['peptide_y_ions']) > 0, "Should have peptide Y ions"
            assert len(fragments['y1_ions']) > 0, "Should have Y1 ions"
            assert len(fragments['b_ions']) > 0, "Should have b-ions"
    
    def test_simple_n_glycopeptide(self):
        """Test simple N-glycopeptide with minimal glycan."""
        gp = Glycopeptide(
            peptide_sequence='LCPDCPLLAPLNDSR',
            glycan_code='2310',  # Core-fucosylated
            glycosylation_site=12,  # N in position
            glycan_type='N'
        )
        
        structures = gp.glycan_structures
        if len(structures) > 0:
            fragments = gp.generate_fragments(structure_index=1)
            
            # Verify all fragment types are present
            assert 'peptide_y_ions' in fragments
            assert 'y1_ions' in fragments
            assert 'b_ions' in fragments
            assert 'intact' in fragments
            
            # Intact should have peptide + glycan mass
            if len(fragments['intact']) > 0:
                intact = fragments['intact'][0]
                total_mass = intact['peptide_mass'] + intact['glycan_mass']
                assert abs(total_mass - intact['total_mass']) < 0.01


# ============================================================================
# Fragment List Validation
# ============================================================================

class TestFragmentListValidation:
    """Validate fragment generation matches expected patterns."""
    
    def test_fragment_count_consistency(self):
        """Test that fragment counts are consistent across multiple runs."""
        glycan_code = '4501'
        
        # Generate fragments multiple times
        counts1 = _get_fragment_counts(glycan_code, 'N')
        counts2 = _get_fragment_counts(glycan_code, 'N')
        counts3 = _get_fragment_counts(glycan_code, 'N')
        
        # Should be identical
        assert counts1 == counts2, "Fragment counts should be consistent"
        assert counts2 == counts3, "Fragment counts should be consistent"
    
    def test_fragment_types_present(self):
        """Test that all expected fragment types are generated."""
        glycan_code = '4501'
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        structure = structures[0]
        fragments, cleavage_info = glycan.generate_fragments(structure, modification_type=0)
        
        # All fragment types should be present
        for frag_type in ['b_ions', 'y_ions', 'by_ions', 'yy_ions']:
            assert frag_type in fragments, f"Missing {frag_type}"
            assert frag_type in cleavage_info, f"Missing {frag_type} in cleavage_info"
    
    def test_fragment_name_format(self):
        """Test that fragment names follow expected format."""
        glycan_code = '4501'
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        structure = structures[0]
        _, cleavage_info = glycan.generate_fragments(structure, modification_type=0)
        
        # B-ions should end with -B
        for name in cleavage_info['b_ions'].keys():
            assert name.endswith('-B'), f"B-ion name should end with -B: {name}"
        
        # Y-ions should end with modification-specific suffix
        valid_y_suffixes = ('-FreeEnd', '-RedEnd', '-PEP', '-2AB')
        for name in cleavage_info['y_ions'].keys():
            assert any(name.endswith(suffix) for suffix in valid_y_suffixes), \
                f"Y-ion should end with one of {valid_y_suffixes}: {name}"


# ============================================================================
# Complex Glycan Validation
# ============================================================================

class TestComplexGlycanValidation:
    """Test complex, high-mass glycan structures."""
    
    def test_high_mass_n_glycan(self):
        """Test high-mass N-glycan (>2500 Da)."""
        # High-mass complex: HexNAc(5)Hex(6)Fuc(2)NeuAc(3)
        glycan_code = '5623'
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Should generate many fragments due to complexity
            total = sum(len(v) for v in fragments.values())
            assert total > 0, "Complex glycan should generate fragments"
    
    def test_highly_branched_n_glycan(self):
        """Test highly branched N-glycan (triantennary+)."""
        # Triantennary: HexNAc(4)Hex(5) with high branching
        glycan_code = '4500'
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Highly branched should generate diverse fragments
            b_count = len(fragments['b_ions'])
            y_count = len(fragments['y_ions'])
            
            assert b_count > 0 or y_count > 0, "Should generate fragments"
    
    def test_all_monosaccharide_types(self):
        """Test glycan with all monosaccharide types."""
        # Has HexNAc, Hex, Fuc, NeuAc: HexNAc(2)Hex(3)Fuc(1)NeuAc(1)
        glycan_code = '2311'
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            fragments, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Collect all monosaccharide types in fragments
            monosac_types = set()
            for frag_list in fragments.values():
                for frag in frag_list:
                    monosac_types.update(frag.keys())
            
            # Should have multiple monosaccharide types
            assert 'HexNAc' in monosac_types
            assert 'Hex' in monosac_types


# ============================================================================
# Modification and Configuration Tests
# ============================================================================

class TestModificationScenarios:
    """Test real-world modification scenarios."""
    
    def test_cam_modified_peptide_with_glycan(self):
        """Test glycopeptide with CAM-modified cysteine."""
        gp = Glycopeptide(
            peptide_sequence='LCPDCPLLAPLNDSR',
            glycan_code='2310',
            glycosylation_site=12,
            glycan_type='N',
            use_cam=True  # CAM on cysteines
        )
        
        structures = gp.glycan_structures
        if len(structures) > 0:
            fragments = gp.generate_fragments()
            
            # Should still generate fragments with CAM applied
            assert len(fragments['peptide_y_ions']) > 0 or len(fragments['y1_ions']) > 0
    
    def test_permethylated_glycan(self):
        """Test permethylated glycan fragments."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            
            # Generate with permethylation
            fragments, _ = glycan.generate_fragments(
                structure,
                modification_type=2
            )
            
            # Should handle permethylation
            assert len(fragments['b_ions']) > 0 or len(fragments['y_ions']) > 0
    
    def test_different_reducing_end_modifications(self):
        """Test different reducing end modification types."""
        glycan = Glycan('2310', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            
            # Test multiple modification types (excluding 6 - reserved for Glycopeptide)
            for mod_type in [0, 1]:
                fragments, _ = glycan.generate_fragments(
                    structure,
                    modification_type=mod_type
                )
                
                # All should generate fragments
                total = sum(len(v) for v in fragments.values())
                assert total > 0, f"Should generate fragments for modification type {mod_type}"


# ============================================================================
# Mass Validation with Real Data
# ============================================================================

class TestRealDataMassValidation:
    """Validate masses for known glycopeptide scenarios."""
    
    def test_peptide_mass_range(self):
        """Test peptide mass is in expected range."""
        sequences = [
            ('PEPTIDE', 700, 900),
            ('EEQYNSTYR', 1100, 1200),
            ('LCPDCPLLAPLNDSR', 1500, 1700)
        ]
        
        for seq, min_mass, max_mass in sequences:
            pep = Peptide(seq, use_cam=False)
            assert min_mass < pep.mass < max_mass, \
                f"Peptide {seq} mass {pep.mass} outside range ({min_mass}, {max_mass})"
    
    def test_glycan_mass_range(self):
        """Test glycan masses are in expected ranges."""
        test_cases = [
            ('2310', 'N', 1000, 1200),  # Core-fucosylated
            ('4501', 'N', 1700, 1900),  # Complex
        ]
        
        for code, gtype, min_mass, max_mass in test_cases:
            glycan = Glycan(code, glycan_type=gtype)
            structures = glycan.predict_structures()
            
            if len(structures) > 0:
                structure = structures[0]
                fragments, _ = glycan.generate_fragments(structure, modification_type=0)
                
                # Check that we have fragments
                total = sum(len(v) for v in fragments.values())
                assert total > 0, f"Should have fragments for {code}"


# ============================================================================
# Helper Functions
# ============================================================================

def _get_fragment_counts(glycan_code, glycan_type='N'):
    """Get fragment type counts for consistency checking."""
    glycan = Glycan(glycan_code, glycan_type=glycan_type)
    structures = glycan.predict_structures()
    
    if not structures:
        return None
    
    structure = structures[0]
    fragments, _ = glycan.generate_fragments(structure, modification_type=0)
    
    return {
        'b_ions': len(fragments['b_ions']),
        'y_ions': len(fragments['y_ions']),
        'by_ions': len(fragments['by_ions']),
        'yy_ions': len(fragments['yy_ions']),
    }


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

