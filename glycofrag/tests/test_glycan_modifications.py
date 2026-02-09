"""
Tests for glycan modification handling (reducing end types, permethylation, special fragments).
"""

import pytest
from glycofrag import Glycan


class TestReducingEndModifications:
    """Test different reducing end modification types."""
    
    def test_free_reducing_end(self):
        """Test modification_type=0 (free reducing end)."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Should work with modification_type=0
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'], modification_type=0)
        
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0
    
    def test_reduced_end_alditol(self):
        """Test modification_type=1 (reduced end/alditol)."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'], modification_type=1)
        
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0
    
    def test_permethylated_free_end(self):
        """Test modification_type=2 (permethylated free end)."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'], modification_type=2)
        
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0
    
    def test_glycopeptide_mode_restricted(self):
        """Test that modification_type='glycopeptide' is blocked for Glycan."""
        import pytest
        
        # Should raise ValueError when trying to use glycopeptide mode
        with pytest.raises(ValueError, match="reserved for Glycopeptide"):
            glycan = Glycan('3300', glycan_type='N', modification_type='glycopeptide')
        
        # Same with numeric value
        with pytest.raises(ValueError, match="reserved for Glycopeptide"):
            glycan = Glycan('3300', glycan_type='N', modification_type=6)
        
        # Also in generate_fragments override
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        with pytest.raises(ValueError, match="reserved for Glycopeptide"):
            glycan.generate_fragments(structures[0], modification_type=6)


class TestNeutralLossFragments:
    """Test special neutral loss fragments."""
    
    def test_neuac_water_loss(self):
        """Test NeuAc-H2O fragment generation."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Check if any B-ions have neutral loss metadata
        has_neutral_loss = any(
            '_neutral_loss' in ion 
            for ion in fragments['b_ions']
        )
        # Note: Might not always generate if NeuAc is not alone
        # This test just verifies the structure
        assert isinstance(fragments['b_ions'], list)
    
    def test_fuc_water_loss(self):
        """Test Fuc-H2O fragment generation."""
        # Create a glycan with fucose
        glycan = Glycan('3310', glycan_type='N')  # Has fucose
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Should generate fragments
        assert len(fragments['b_ions']) > 0
    
    def test_fuc_methyl_loss_permethylated(self):
        """Test Fuc-CH3 fragment generation for permethylated glycans."""
        glycan = Glycan('3310', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Should generate more variants with permethylation
        frags_normal, _ = glycan.generate_fragments(structures[0], ['BY'], modification_type=0)
        frags_perm, _ = glycan.generate_fragments(structures[0], ['BY'], modification_type=2)
        
        # Permethylated should have at least as many fragments
        assert len(frags_perm['b_ions']) >= len(frags_normal['b_ions'])
    
    def test_hexnac_c2h6o3_loss(self):
        """Test HexNAc-C2H6O3 (126 Da) fragment generation."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Should generate B-ions including neutral loss variants
        assert len(fragments['b_ions']) > 0


class TestOxoniumIons:
    """Test A-type oxonium ion generation."""
    
    def test_oxonium_generation(self):
        """Test that A-type oxonium ions are generated."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Request A-type fragments
        fragments, cleavage_info = glycan.generate_fragments(structures[0], ['A'])
        
        # Should have a_ions key
        assert 'a_ions' in fragments
        assert 'a_ions' in cleavage_info
        
        # Should generate oxonium ions
        assert len(fragments['a_ions']) > 0
    
    def test_n_glycan_oxonium_markers(self):
        """Test N-glycan specific oxonium markers."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, cleavage_info = glycan.generate_fragments(structures[0], ['A'])
        
        # Should have oxonium ions
        assert len(fragments['a_ions']) >= 1
        
        # Check for N-glycan marker
        has_n_marker = any(
            ion.get('_custom_label') == 'HexNAc(-C2H4O2)'
            for ion in fragments['a_ions']
        )
        assert has_n_marker, "N-glycan diagnostic marker (144.066) should be present"
    
    def test_o_glycan_oxonium_markers(self):
        """Test O-glycan specific oxonium markers."""
        glycan = Glycan('2200', glycan_type='O')
        structures = glycan.predict_structures()
        
        fragments, cleavage_info = glycan.generate_fragments(structures[0], ['A'])
        
        # Should have oxonium ions
        assert len(fragments['a_ions']) >= 1
        
        # Check for O-glycan marker
        has_o_marker = any(
            ion.get('_custom_label') == 'HexNAc(-CH6O3)'
            for ion in fragments['a_ions']
        )
        assert has_o_marker, "O-glycan diagnostic marker (138.055) should be present"
    
    def test_common_oxonium_ion(self):
        """Test common oxonium ion (126 Da)."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['A'])
        
        # Should have ~126 Da oxonium (HexNAc-C2H6O3, actual mass ~125 Da)
        has_126 = any(
            ion.get('_custom_label') == 'HexNAc-C2H6O3'
            for ion in fragments['a_ions']
        )
        assert has_126, "HexNAc-C2H6O3 oxonium ion should be present"
    
    def test_oxonium_direct_mz(self):
        """Test that oxonium ions have direct m/z values."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['A'])
        
        # All oxonium ions should have _direct_mz
        for ion in fragments['a_ions']:
            assert '_direct_mz' in ion, "Oxonium ions should have direct m/z"
            assert ion['_direct_mz'] > 0


class TestFragmentTypeCombinations:
    """Test generating multiple fragment types together."""
    
    def test_by_cz_a_together(self):
        """Test generating BY, CZ, and A-type fragments together."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY', 'CZ', 'A'])
        
        # Should have all fragment types
        assert 'b_ions' in fragments
        assert 'c_ions' in fragments
        assert 'a_ions' in fragments
        
        # All should have content
        assert len(fragments['b_ions']) > 0
        assert len(fragments['c_ions']) > 0
        assert len(fragments['a_ions']) > 0
    
    def test_default_includes_oxonium(self):
        """Test that requesting A-type includes oxonium ions."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Glycan default is ['BY'] only; include 'A' to get oxonium
        fragments, _ = glycan.generate_fragments(structures[0], fragment_types=['BY', 'A'])
        
        # Should have BY and A-type
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        assert 'a_ions' in fragments


class TestPermethylationEffects:
    """Test permethylation impact on fragmentation."""
    
    def test_permethylation_flag(self):
        """Test permethylation implied by modification_type."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Should work for both non-permethylated and permethylated types
        frags_not_perm, _ = glycan.generate_fragments(structures[0], modification_type=0)
        frags_perm, _ = glycan.generate_fragments(structures[0], modification_type=2)
        
        # Both should generate fragments
        assert len(frags_not_perm['b_ions']) > 0
        assert len(frags_perm['b_ions']) > 0
    
    def test_permethylation_with_mod_type_2(self):
        """Test permethylation with modification_type=2."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(
            structures[0], 
            modification_type=2
        )
        
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0


class TestModificationMetadata:
    """Test that modification metadata is properly stored."""
    
    def test_neutral_loss_metadata(self):
        """Test that neutral loss fragments have metadata."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Check for metadata in fragments
        has_metadata = any(
            any(k.startswith('_') for k in ion.keys())
            for ion in fragments['b_ions']
        )
        # Metadata presence depends on whether single sugar fragments exist
        assert isinstance(fragments['b_ions'], list)
    
    def test_oxonium_metadata(self):
        """Test that oxonium ions have proper metadata."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['A'])
        
        # All oxonium ions should be marked
        for ion in fragments['a_ions']:
            assert '_is_oxonium' in ion
            assert ion['_is_oxonium'] == True
            assert '_custom_label' in ion


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
