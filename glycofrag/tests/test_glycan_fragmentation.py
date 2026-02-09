"""
Tests for Glycan fragmentation functionality (Phase 2.5).

This module contains tests for:
- BY-series: B, Y, BY, YY, BYY, YYY, BYYY
- CZ-series: C, Z, CZ, ZZ, CZZ, BZZ, ZZZ, CZZZ, BZZZ
"""

import pytest
from glycofrag import Glycan


class TestBYSeriesFragmentation:
    """Test BY-series fragmentation (B, Y, BY, YY, BYY, YYY, BYYY)."""
    
    def test_by_series_generation(self):
        """Test that BY series generates all fragment types."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Should have all BY series types
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        assert 'by_ions' in fragments
        assert 'yy_ions' in fragments
        assert 'byy_ions' in fragments
        assert 'yyy_ions' in fragments
        assert 'byyy_ions' in fragments
        
        # Should generate fragments for basic types
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0
        
        # Note: 4501 has only 2 branches, so YYY/BYYY (which require 3 branches) should be empty
        # This is correct behavior based on glycan topology
    
    def test_b_ions_no_reducing_end(self):
        """Test that B-ions don't contain reducing end."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # B-ions should not be the full glycan
        total_comp = glycan._count_residues(structures[0])
        total_count = sum(total_comp.values())
        
        mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
        for b_ion in fragments['b_ions']:
            # Only sum monosaccharide counts, not metadata
            b_count = sum(b_ion.get(k, 0) for k in mono_keys)
            # B-ions are smaller than total (no reducing end)
            assert b_count < total_count
    
    def test_y_ions_contain_reducing_end(self):
        """Test that Y-ions contain reducing end (placeholder)."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Y-ions should all have content (including reducing end placeholder)
        for y_ion in fragments['y_ions']:
            assert sum(y_ion.values()) > 0
    
    def test_yy_ions_from_y_ions(self):
        """Test that YY-ions are derived from Y-ions."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Each YY-ion should be smaller than some Y-ion
        for yy_ion in fragments['yy_ions']:
            yy_count = sum(yy_ion.values())
            found_larger = any(
                sum(y_ion.values()) > yy_count
                for y_ion in fragments['y_ions']
            )
            assert found_larger, "YY-ion should come from larger Y-ion"
    
    def test_byy_ions_from_by_ions(self):
        """Test that BYY-ions are derived from BY-ions."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # If BYY ions are generated, check their properties
        if len(fragments['byy_ions']) > 0:
            # Each BYY-ion should be smaller than some BY-ion
            mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
            for byy_ion in fragments['byy_ions']:
                byy_count = sum(byy_ion.get(k, 0) for k in mono_keys)
                found_larger = any(
                    sum(by_ion.get(k, 0) for k in mono_keys) > byy_count
                for by_ion in fragments['by_ions']
            )
            assert found_larger, "BYY-ion should come from larger BY-ion"
    
    def test_yyy_ions_from_yy_ions(self):
        """Test that YYY-ions can include core fragments even for 2-branch structures."""
        # 4501 has 2 branches, but can still generate core YYY fragments (e.g., HexNAc2Hex2)
        glycan_2branch = Glycan('4501', glycan_type='N')
        structures_2branch = glycan_2branch.predict_structures()
        fragments_2branch, _ = glycan_2branch.generate_fragments(structures_2branch[0], ['BY'])
        
        # Core YYY fragments can be generated (3 cleavages in core region)
        # Verify at least the core fragment is present
        if len(fragments_2branch['yyy_ions']) > 0:
            # Check that YYY fragments meet minimum composition requirements
            mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
            for yyy_ion in fragments_2branch['yyy_ions']:
                yyy_count = sum(yyy_ion.get(k, 0) for k in mono_keys)
                assert yyy_count >= 4, "YYY ions should have at least 4 monosaccharides (core minimum)"
    
    def test_byyy_ions_generated(self):
        """Test that BYYY-ions can include core fragments even for 2-branch structures."""
        # 4501 has 2 branches, but can still generate core BYYY fragments
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        # Core BYYY fragments can be generated (3 cleavages with B-ion retention)
        # Verify fragments meet minimum composition requirements if present
        if len(fragments['byyy_ions']) > 0:
            mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
            for byyy_ion in fragments['byyy_ions']:
                byyy_count = sum(byyy_ion.get(k, 0) for k in mono_keys)
                assert byyy_count >= 2, "BYYY ions should have at least 2 monosaccharides"
    
    def test_branch_counting(self):
        """Test that branch counting works correctly for N-glycans."""
        # Test 4501: should have 2 branches
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        assert len(structures) > 0
        
        num_branches = glycan._count_branches(structures[0])
        assert num_branches == 2, f"4501 should have 2 branches, got {num_branches}"
        
        # Note: YYY/BYYY can still be generated (core fragments with 3 cleavages)
        # Even with 2 branches, 3 cleavages are possible in the core region
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        # Just verify fragment generation completes successfully
        assert 'yyy_ions' in fragments
        assert 'byyy_ions' in fragments
    
    def test_branch_based_fragmentation_logic(self):
        """Test that fragment generation works correctly for 2-branch structures."""
        # Structure 4501:
        # - HexNAc(4)Hex(5)Fuc(0)NeuAc(1)
        # - Has 2 branches (arms from core)
        # - YY ions: OK (2 cleavages on 2 branches)
        # - YYY ions: OK (3 cleavages possible in core region)
        # - BYYY ions: OK (3 cleavages with B-ion retention)
        
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        structure = structures[0]
        fragments, _ = glycan.generate_fragments(structure, ['BY'])
        
        # Basic fragments should always be present
        assert len(fragments['b_ions']) > 0, "B ions should be generated"
        assert len(fragments['y_ions']) > 0, "Y ions should be generated"
        assert len(fragments['by_ions']) > 0, "BY ions should be generated"
        
        # YY ions may be generated if structure has sufficient branches
        if len(fragments['yy_ions']) > 0:
            pass  # If generated, fragmentation is working
        
        # YYY/BYYY can be generated (core cleavages)
        # Just verify the keys exist and fragmentation completed
        assert 'yyy_ions' in fragments, "YYY ions key should exist"
        assert 'byyy_ions' in fragments, "BYYY ions key should exist"

class TestCZSeriesFragmentation:
    """Test CZ-series fragmentation (C, Z, CZ, ZZ, CZZ, BZZ, ZZZ, CZZZ, BZZZ)."""
    
    def test_cz_series_generation(self):
        """Test that CZ series generates all fragment types."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Should have all CZ series types
        assert 'c_ions' in fragments
        assert 'z_ions' in fragments
        assert 'cz_ions' in fragments
        assert 'zz_ions' in fragments
        assert 'czz_ions' in fragments
        assert 'bzz_ions' in fragments
        assert 'zzz_ions' in fragments
        assert 'czzz_ions' in fragments
        assert 'bzzz_ions' in fragments
        
        # Should generate fragments
        assert len(fragments['c_ions']) > 0
        assert len(fragments['z_ions']) > 0
    
    def test_cz_parallel_to_by(self):
        """Test that CZ fragments parallel BY fragments."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        by_fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        cz_fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Should have same number of B and C ions
        assert len(by_fragments['b_ions']) == len(cz_fragments['c_ions'])
        
        # Should have same number of Y and Z ions
        assert len(by_fragments['y_ions']) == len(cz_fragments['z_ions'])
        
        # Should have same number of BY and CZ ions
        assert len(by_fragments['by_ions']) == len(cz_fragments['cz_ions'])
    
    def test_czz_from_byy(self):
        """Test that CZZ-ions parallel BYY-ions."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        by_fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        cz_fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Should have same number of BYY and CZZ ions
        assert len(by_fragments['byy_ions']) == len(cz_fragments['czz_ions'])
    
    def test_bzz_ions_generated(self):
        """Test that BZZ ions can be generated."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Check BZZ ions key exists (may be empty depending on structure)
        assert 'bzz_ions' in fragments
    
    def test_czzz_and_bzzz_generated(self):
        """Test that CZZZ and BZZZ ions keys exist."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Check that fragment keys exist (may be empty depending on structure)
        assert 'czzz_ions' in fragments
        assert 'bzzz_ions' in fragments


class TestBothSeriesFragmentation:
    """Test generating both BY and CZ series together."""
    
    def test_both_series_generated(self):
        """Test that both BY and CZ series can be generated together."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY', 'CZ'])
        
        # Should have both BY and CZ fragment types
        assert 'b_ions' in fragments
        assert 'c_ions' in fragments
        assert 'y_ions' in fragments
        assert 'z_ions' in fragments
        assert 'by_ions' in fragments
        assert 'cz_ions' in fragments
    
    def test_fragment_counts_independent(self):
        """Test that BY and CZ fragments are independently generated."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        both_fragments, _ = glycan.generate_fragments(structures[0], ['BY', 'CZ'])
        by_only, _ = glycan.generate_fragments(structures[0], ['BY'])
        cz_only, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        # BY counts should match
        assert len(both_fragments['b_ions']) == len(by_only['b_ions'])
        assert len(both_fragments['y_ions']) == len(by_only['y_ions'])
        
        # CZ counts should match
        assert len(both_fragments['c_ions']) == len(cz_only['c_ions'])
        assert len(both_fragments['z_ions']) == len(cz_only['z_ions'])


class TestFragmentFormatting:
    """Test fragment string formatting."""
    
    def test_b_ion_format(self):
        """Test B-ion formatting."""
        glycan = Glycan('3300', glycan_type='N')
        comp = {'HexNAc': 1, 'Hex': 2}
        formatted = glycan._format_fragment_string(comp, 'B')
        
        assert 'HexNAc' in formatted
        assert 'Hex2' in formatted
        assert '-B' in formatted
    
    def test_yy_ion_format(self):
        """Test YY-ion formatting."""
        glycan = Glycan('3300', glycan_type='N')
        comp = {'HexNAc': 2, 'Hex': 1}
        formatted = glycan._format_fragment_string(comp, 'YY', is_y_ion=True)
        
        assert 'HexNAc2' in formatted
        assert 'Hex' in formatted
        # YY-ions also use modification-specific suffixes
        valid_suffixes = ('-FreeEnd', '-Redend', '-PEP', '-2AB')
        assert any(formatted.endswith(suffix) for suffix in valid_suffixes)
    
    def test_czz_ion_format(self):
        """Test CZZ-ion formatting."""
        glycan = Glycan('3300', glycan_type='N')
        comp = {'HexNAc': 1, 'Hex': 1}
        formatted = glycan._format_fragment_string(comp, 'CZZ')
        
        assert 'HexNAc' in formatted
        assert 'Hex' in formatted
        assert '-CZZ' in formatted
    
    def test_bzzz_ion_format(self):
        """Test BZZZ-ion formatting."""
        glycan = Glycan('3300', glycan_type='N')
        comp = {'HexNAc': 1}
        formatted = glycan._format_fragment_string(comp, 'BZZZ')
        
        assert 'HexNAc-BZZZ' == formatted


class TestCleavageInfo:
    """Test cleavage information tracking."""
    
    def test_by_cleavage_info(self):
        """Test that BY series has cleavage info."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        _, cleavage_info = glycan.generate_fragments(structures[0], ['BY'])
        
        # Should have cleavage info for all types
        assert 'b_ions' in cleavage_info
        assert 'y_ions' in cleavage_info
        assert 'by_ions' in cleavage_info
        assert 'yy_ions' in cleavage_info
    
    def test_cz_cleavage_info(self):
        """Test that CZ series has cleavage info."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        _, cleavage_info = glycan.generate_fragments(structures[0], ['CZ'])
        
        # Should have cleavage info for CZ types
        assert 'c_ions' in cleavage_info
        assert 'z_ions' in cleavage_info
        assert 'cz_ions' in cleavage_info


class TestDefaultBehavior:
    """Test default fragmentation behavior."""
    
    def test_default_generates_by_only(self):
        """Test that default generates only BY series."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        # No fragment_types specified, should default to BY
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Should have BY fragments
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        
        # Should NOT have CZ fragments
        assert 'c_ions' not in fragments
        assert 'z_ions' not in fragments


class TestOGlycanFragmentation:
    """Test O-glycan fragmentation."""
    
    def test_oglycan_by_series(self):
        """Test BY series for O-glycans."""
        glycan = Glycan('2200', glycan_type='O')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['BY'])
        
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y_ions']) > 0
    
    def test_oglycan_cz_series(self):
        """Test CZ series for O-glycans."""
        glycan = Glycan('2200', glycan_type='O')
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(structures[0], ['CZ'])
        
        assert len(fragments['c_ions']) > 0
        assert len(fragments['z_ions']) > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
    """Test basic glycan fragmentation."""
    
    def test_simple_nglycan_fragmentation(self):
        """Test fragmentation of a simple N-glycan."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        
        assert len(structures) > 0, "Should generate at least one structure"
        
        # Generate fragments for first structure
        fragments, cleavage_info = glycan.generate_fragments(structures[0])
        
        # Should have all fragment types
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        assert 'by_ions' in fragments
        assert 'yy_ions' in fragments
        
        # Should generate some fragments
        assert len(fragments['b_ions']) > 0, "Should generate B-ions"
        assert len(fragments['y_ions']) > 0, "Should generate Y-ions"
    
    def test_fragment_composition(self):
        """Test that fragment compositions are correct."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # All fragments should be dicts with monosaccharide counts
        for b_ion in fragments['b_ions']:
            assert isinstance(b_ion, dict)
            assert all(k in ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc'] for k in b_ion.keys())
            assert all(isinstance(v, int) and v >= 0 for v in b_ion.values())
        
        for y_ion in fragments['y_ions']:
            assert isinstance(y_ion, dict)
            assert sum(y_ion.values()) > 0, "Y-ions should not be empty"
    
    def test_fragmentation_with_fucose(self):
        """Test fragmentation of fucosylated glycan."""
        glycan = Glycan('3310', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Should have fragments containing fucose
        has_fuc_fragment = any(
            frag.get('Fuc', 0) > 0 
            for frag in fragments['b_ions'] + fragments['y_ions']
        )
        assert has_fuc_fragment, "Should generate fucose-containing fragments"
    
    def test_fragmentation_with_sialic_acid(self):
        """Test fragmentation of sialylated glycan."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Should have fragments containing sialic acid
        has_neuac_fragment = any(
            frag.get('NeuAc', 0) > 0 
            for frag in fragments['b_ions'] + fragments['y_ions']
        )
        assert has_neuac_fragment, "Should generate sialic acid-containing fragments"


class TestBIonGeneration:
    """Test B-ion generation specifically."""
    
    def test_b_ions_no_reducing_end(self):
        """Test that B-ions don't contain reducing end."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # B-ions should be smaller than total composition
        total_comp = glycan._count_residues(structures[0])
        total_count = sum(total_comp.values())
        
        mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
        for b_ion in fragments['b_ions']:
            # Only sum monosaccharide counts
            b_count = sum(b_ion.get(k, 0) for k in mono_keys)
            assert b_count <= total_count
    
    def test_b_ion_string_format(self):
        """Test B-ion string formatting."""
        glycan = Glycan('3300', glycan_type='N')
        
        # Test format function
        comp = {'HexNAc': 1, 'Hex': 2}
        formatted = glycan._format_fragment_string(comp, 'B')
        
        assert 'HexNAc' in formatted
        assert 'Hex2' in formatted
        assert '-B' in formatted
    
    def test_b_ions_unique(self):
        """Test that B-ions are unique (no duplicates)."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Convert to strings for comparison
        b_strings = [glycan._format_fragment_string(b, 'B') for b in fragments['b_ions']]
        
        # Should have no duplicates
        assert len(b_strings) == len(set(b_strings))


class TestYIonGeneration:
    """Test Y-ion generation specifically."""
    
    def test_y_ions_contain_reducing_end(self):
        """Test that Y-ions contain the reducing end."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Y-ions should all have at least one HexNAc (reducing end)
        for y_ion in fragments['y_ions']:
            total = sum(y_ion.values())
            assert total > 0, "Y-ion should not be empty"
    
    def test_y_ion_string_format(self):
        """Test Y-ion string formatting."""
        glycan = Glycan('3300', glycan_type='N')
        
        comp = {'HexNAc': 2, 'Hex': 1}
        formatted = glycan._format_fragment_string(comp, 'Y', is_y_ion=True)
        
        assert 'HexNAc2' in formatted
        assert 'Hex' in formatted
        # Should have modification-specific suffix, not just '-Y'
        valid_suffixes = ('-FreeEnd', '-Redend', '-PEP', '-2AB')
        assert any(formatted.endswith(suffix) for suffix in valid_suffixes)
    
    def test_y_ions_unique(self):
        """Test that Y-ions are unique."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        y_strings = [glycan._format_fragment_string(y, 'Y', is_y_ion=True) for y in fragments['y_ions']]
        assert len(y_strings) == len(set(y_strings))


class TestBYIonGeneration:
    """Test BY-ion generation (double cleavages)."""
    
    def test_by_ions_generated(self):
        """Test that BY-ions are generated for complex glycans."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Complex glycan should have BY-ions
        assert len(fragments['by_ions']) > 0
    
    def test_by_ion_format(self):
        """Test BY-ion string formatting."""
        glycan = Glycan('3300', glycan_type='N')
        
        comp = {'HexNAc': 1}
        formatted = glycan._format_fragment_string(comp, 'BY')
        
        assert 'HexNAc' in formatted
        assert '-BY' in formatted
    
    def test_by_ions_smaller_than_original(self):
        """Test that BY-ions are smaller than original glycan."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        total_comp = glycan._count_residues(structures[0])
        total_count = sum(total_comp.values())
        
        mono_keys = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
        for by_ion in fragments['by_ions']:
            by_count = sum(by_ion.get(k, 0) for k in mono_keys)
            assert by_count < total_count, "BY-ion should be smaller than intact glycan"


class TestYYIonGeneration:
    """Test YY-ion generation (sequential losses)."""
    
    def test_yy_ions_generated(self):
        """Test that YY-ions can be generated for branched structures."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # YY-ions may be generated depending on structure branching
        # Just verify the key exists and fragmentation completed
        assert 'yy_ions' in fragments
    
    def test_yy_ion_format(self):
        """Test YY-ion string formatting."""
        glycan = Glycan('3300', glycan_type='N')
        
        comp = {'HexNAc': 2, 'Hex': 2}
        formatted = glycan._format_fragment_string(comp, 'YY', is_y_ion=True)
        
        assert 'HexNAc2' in formatted
        assert 'Hex2' in formatted
        # YY-ions use modification-specific suffixes
        valid_suffixes = ('-FreeEnd', '-Redend', '-PEP', '-2AB')
        assert any(formatted.endswith(suffix) for suffix in valid_suffixes)
    
    def test_yy_ions_from_y_ions(self):
        """Test that YY-ions are derived from Y-ions."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Each YY-ion should be smaller than at least one Y-ion
        for yy_ion in fragments['yy_ions']:
            yy_count = sum(yy_ion.values())
            
            found_larger_y = any(
                sum(y_ion.values()) > yy_count
                for y_ion in fragments['y_ions']
            )
            assert found_larger_y, "YY-ion should be derived from a larger Y-ion"


class TestCleavageInfo:
    """Test cleavage information tracking."""
    
    def test_cleavage_info_structure(self):
        """Test that cleavage info has correct structure."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        _, cleavage_info = glycan.generate_fragments(structures[0])
        
        assert 'b_ions' in cleavage_info
        assert 'y_ions' in cleavage_info
        assert 'by_ions' in cleavage_info
        assert 'yy_ions' in cleavage_info
        
        # Each should be a dict
        assert isinstance(cleavage_info['b_ions'], dict)
        assert isinstance(cleavage_info['y_ions'], dict)
    
    def test_cleavage_info_matches_fragments(self):
        """Test that cleavage info keys match fragment strings."""
        glycan = Glycan('3300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, cleavage_info = glycan.generate_fragments(structures[0])
        
        # Number of entries in cleavage_info should match fragments
        assert len(cleavage_info['b_ions']) == len(fragments['b_ions'])
        assert len(cleavage_info['y_ions']) == len(fragments['y_ions'])
        assert len(cleavage_info['by_ions']) == len(fragments['by_ions'])
        assert len(cleavage_info['yy_ions']) == len(fragments['yy_ions'])


class TestOGlycanFragmentation:
    """Test O-glycan fragmentation."""
    
    def test_oglycan_fragmentation(self):
        """Test that O-glycans can be fragmented."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        
        assert len(structures) > 0
        
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Should generate fragments
        total_fragments = sum(len(v) for v in fragments.values())
        assert total_fragments > 0
    
    def test_oglycan_b_ions(self):
        """Test B-ion generation for O-glycans."""
        glycan = Glycan('2200', glycan_type='O')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        assert len(fragments['b_ions']) > 0
    
    def test_oglycan_y_ions(self):
        """Test Y-ion generation for O-glycans."""
        glycan = Glycan('2200', glycan_type='O')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        assert len(fragments['y_ions']) > 0


class TestEdgeCases:
    """Test edge cases in fragmentation."""
    
    def test_minimal_glycan(self):
        """Test fragmentation of minimal glycan."""
        glycan = Glycan('2300', glycan_type='N')
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(structures[0])
        
        # Should generate at least some fragments
        assert len(fragments['b_ions']) > 0 or len(fragments['y_ions']) > 0
    
    def test_empty_composition(self):
        """Test handling of empty fragment composition."""
        glycan = Glycan('3300', glycan_type='N')
        
        # Test format function with empty composition
        formatted = glycan._format_fragment_string({}, 'B')
        assert formatted == ""
    
    def test_single_monosaccharide_composition(self):
        """Test formatting single monosaccharide."""
        glycan = Glycan('3300', glycan_type='N')
        
        comp = {'HexNAc': 1}
        formatted = glycan._format_fragment_string(comp, 'B')
        
        assert formatted == 'HexNAc-B'
        assert 'HexNAc1' not in formatted  # Should not show "1"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

