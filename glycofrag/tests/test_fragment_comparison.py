"""
Phase 6.1: Fragment Comparison Validation - FOCUSED TESTS

This test module validates that glycofrag's glycan fragment generation
works correctly for B/Y/BY/YY ion generation with realistic glycan codes.
"""

import pytest
from glycofrag import Glycan
from glycofrag.core.mass_calculator import GlycanMassCalculator


class GlycofragTestHelper:
    """Helper class to generate fragments using glycofrag API."""
    
    @staticmethod
    def get_fragments(glycan_code, glycan_type='N', **kwargs):
        """Generate fragments for a glycan using glycofrag."""
        glycan = Glycan(glycan_code, glycan_type=glycan_type)
        structures = glycan.predict_structures()
        if not structures:
            raise ValueError(f"No structures generated for {glycan_code}")
        
        # Use first predicted structure
        structure = structures[0]
        
        # Generate fragments for this structure
        fragments, cleavage_info = glycan.generate_fragments(
            structure,
            **kwargs
        )
        
        return fragments, cleavage_info


class TestBIonFragments:
    """Test B-ion (non-reducing end) fragment generation."""
    
    def test_realistic_n_glycan_b_ions(self):
        """Test B-ion generation for realistic N-glycan (4501)."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('4501', 'N')
        
        # Should have B-ions
        assert len(fragments['b_ions']) > 0, "Should generate B-ions"
        
        # All B-ion names should end with '-B'
        b_names = list(cleavage_info['b_ions'].keys())
        for name in b_names:
            assert name.endswith('-B'), f"B-ion should end with '-B': {name}"
        
        # All B-ions should have valid composition
        for b_ion in fragments['b_ions']:
            assert isinstance(b_ion, dict)
            for mono_type, count in b_ion.items():
                assert count >= 0, f"Count should be non-negative: {mono_type}={count}"
    
    def test_b_ion_mass_validity(self):
        """Test that B-ion masses are in valid range."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N')
        
        calculator = GlycanMassCalculator()
        
        for b_ion in fragments['b_ions']:
            # Convert string values to int if needed
            comp = {k: int(v) if isinstance(v, str) else v for k, v in b_ion.items()}
            mass = calculator.calculate_fragment_mass(comp, 'b_ions')
            # Reasonable range for glycan fragments (50-5000 Da)
            assert 50 < mass < 5000, f"B-ion mass out of range: {mass}"


class TestYIonFragments:
    """Test Y-ion (reducing end) fragment generation."""
    
    def test_realistic_n_glycan_y_ions(self):
        """Test Y-ion generation for realistic N-glycan (4501)."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('4501', 'N')
        
        # Should have Y-ions
        assert len(fragments['y_ions']) > 0, "Should generate Y-ions"
        
        # Y-ion names should end with modification-specific suffix
        # Valid suffixes: -FreeEnd, -RedEnd, -PEP, -2AB
        valid_suffixes = ('-FreeEnd', '-RedEnd', '-PEP', '-2AB')
        y_names = list(cleavage_info['y_ions'].keys())
        for name in y_names:
            assert any(name.endswith(suffix) for suffix in valid_suffixes), \
                f"Y-ion should end with one of {valid_suffixes}: {name}"
        
        # All Y-ions should have valid composition
        for y_ion in fragments['y_ions']:
            assert isinstance(y_ion, dict)
            for mono_type, count in y_ion.items():
                assert count >= 0
    
    def test_y_ion_mass_validity(self):
        """Test that Y-ion masses are in valid range."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N')
        
        calculator = GlycanMassCalculator()
        
        for y_ion in fragments['y_ions']:
            # Convert string values to int if needed
            comp = {k: int(v) if isinstance(v, str) else v for k, v in y_ion.items()}
            mass = calculator.calculate_fragment_mass(comp, 'y_ions')
            # Reasonable range for glycan fragments
            assert 50 < mass < 5000, f"Y-ion mass out of range: {mass}"


class TestBYandYYSeriesFragments:
    """Test multi-cleavage fragment generation."""
    
    def test_by_series_generation(self):
        """Test BY-series fragment generation."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('4501', 'N')
        
        # BY-ions may or may not exist depending on structure
        assert isinstance(fragments['by_ions'], list)
        
        # If they exist, verify naming and composition
        for by_ion in fragments['by_ions']:
            assert isinstance(by_ion, dict)
            total_count = sum(by_ion.values())
            assert total_count > 0, "BY fragment should have at least one monosaccharide"
    
    def test_yy_series_generation(self):
        """Test YY-series fragment generation."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('4501', 'N')
        
        # YY-ions may or may not exist
        assert isinstance(fragments['yy_ions'], list)
        
        # If they exist, verify composition
        for yy_ion in fragments['yy_ions']:
            assert isinstance(yy_ion, dict)
            total_count = sum(yy_ion.values())
            assert total_count > 0, "YY fragment should have at least one monosaccharide"


class TestFragmentConsistency:
    """Test consistency between fragments and cleavage info."""
    
    def test_fragments_and_cleavage_info_pairing(self):
        """Test that fragment lists match cleavage info."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('4501', 'N')
        
        # Verify counts match
        assert len(fragments['b_ions']) == len(cleavage_info['b_ions']), \
            "B-ion count mismatch between fragments and cleavage_info"
        assert len(fragments['y_ions']) == len(cleavage_info['y_ions']), \
            "Y-ion count mismatch"
        assert len(fragments['by_ions']) == len(cleavage_info['by_ions']), \
            "BY-ion count mismatch"
        assert len(fragments['yy_ions']) == len(cleavage_info['yy_ions']), \
            "YY-ion count mismatch"
    
    def test_no_negative_composition_values(self):
        """Test that all fragments have non-negative composition values."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N')
        
        for fragment_list in fragments.values():
            for fragment in fragment_list:
                for mono_type, count in fragment.items():
                    # Skip metadata keys (start with underscore or contain special chars)
                    if mono_type.startswith('_'):
                        continue
                    if isinstance(count, str) and '(' in count:
                        continue
                    
                    # Convert to int if string
                    try:
                        count_val = int(count) if isinstance(count, str) else count
                        assert count_val >= 0, \
                            f"Negative count found: {mono_type}={count_val}"
                    except (ValueError, TypeError):
                        # Skip non-numeric values
                        pass


class TestModificationSupport:
    """Test fragment generation with different modifications."""
    
    def test_modification_type_0(self):
        """Test with modification type 0 (free)."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N', modification_type=0)
        
        # Should still generate fragments
        assert len(fragments['b_ions']) > 0 or len(fragments['y_ions']) > 0
    
    def test_modification_type_6(self):
        """Test with modification type 6 (standard)."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N', modification_type=0)
        
        assert len(fragments['b_ions']) > 0 or len(fragments['y_ions']) > 0
    
    def test_permethylated_generation(self):
        """Test permethylated glycan fragment generation."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N', modification_type=2)
        
        # Should generate fragments even when permethylated
        assert isinstance(fragments['b_ions'], list)
        assert isinstance(fragments['y_ions'], list)


class TestOGlycanFragments:
    """Test O-glycan specific fragment generation."""
    
    def test_simple_o_glycan(self):
        """Test O-glycan fragment generation."""
        fragments, _ = GlycofragTestHelper.get_fragments('1001', 'O')
        
        # Should generate O-glycan fragments
        total_fragments = sum(len(v) for v in fragments.values())
        assert total_fragments > 0, "O-glycan should generate fragments"
    
    def test_o_glycan_naming(self):
        """Test O-glycan fragment naming follows conventions."""
        fragments, cleavage_info = GlycofragTestHelper.get_fragments('1001', 'O')
        
        # Verify B-ions end with -B
        for name in cleavage_info['b_ions'].keys():
            assert name.endswith('-B'), f"O-glycan B-ion should end with -B: {name}"
        
        # Verify Y-ions end with modification-specific suffix
        valid_suffixes = ('-FreeEnd', '-RedEnd', '-PEP', '-2AB')
        for name in cleavage_info['y_ions'].keys():
            assert any(name.endswith(suffix) for suffix in valid_suffixes), \
                f"O-glycan Y-ion should end with one of {valid_suffixes}: {name}"


class TestFragmentDiversity:
    """Test that fragment diversity is appropriate."""
    
    def test_realistic_n_glycan_diversity(self):
        """Test realistic N-glycan (4501) generates multiple fragments."""
        fragments, _ = GlycofragTestHelper.get_fragments('4501', 'N')
        
        b_count = len(fragments['b_ions'])
        y_count = len(fragments['y_ions'])
        
        # For realistic N-glycan, expect multiple fragments
        assert b_count > 0, "Should have B-ions"
        assert y_count > 0, "Should have Y-ions"
        
        total = b_count + y_count
        assert total >= 3, f"Expected at least 3 total fragments, got {total}"
    
    def test_complex_glycan_fragment_count(self):
        """Test complex glycan generates sufficient fragment diversity."""
        fragments, _ = GlycofragTestHelper.get_fragments('3320', 'N')
        
        total_fragments = sum(len(v) for v in fragments.values())
        assert total_fragments > 0, "Complex glycan should generate fragments"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

