"""
Tests for the Peptide class and peptide fragmentation functionality.

This module contains tests for:
- Peptide initialization
- b-ion generation
- y-ion generation
- c-ion generation
- z-ion generation
- Modification handling
- Fragment retrieval
"""

import pytest
from glycofrag.peptide import Peptide


class TestPeptideInitialization:
    """Test Peptide class initialization."""
    
    def test_simple_peptide(self):
        """Test initializing a simple peptide."""
        peptide = Peptide('PEPTIDE')
        assert peptide.sequence == 'PEPTIDE'
        assert len(peptide) == 7
        assert peptide.use_cam is True
    
    def test_peptide_with_cam(self):
        """Test peptide with CAM modification."""
        peptide = Peptide('PEPTIDEC', use_cam=True)
        assert peptide.use_cam is True
        # CAM adds 57.021464 to cysteine
        assert peptide.mass > 0
    
    def test_peptide_without_cam(self):
        """Test peptide without CAM modification."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        assert peptide.use_cam is False
    
    def test_lowercase_sequence(self):
        """Test that sequences are converted to uppercase."""
        peptide = Peptide('peptide')
        assert peptide.sequence == 'PEPTIDE'
    
    def test_peptide_mass_calculated(self):
        """Test that peptide mass is calculated during initialization."""
        peptide = Peptide('AAA')
        assert peptide.mass > 0
        # 3 alanines: 3 * 71.037114 = 213.111342
        # Plus reducing end modification
        assert peptide.mass > 200


class TestBIonGeneration:
    """Test b-ion generation."""
    
    def test_simple_b_ions(self):
        """Test generating b-ions for a simple peptide."""
        peptide = Peptide('AAA', use_cam=False)
        b_ions = peptide.generate_b_ions()
        
        assert len(b_ions) == 3
        assert b_ions[0]['fragment_name'] == 'b1'
        assert b_ions[1]['fragment_name'] == 'b2'
        assert b_ions[2]['fragment_name'] == 'b3'
    
    def test_b_ion_sequences(self):
        """Test that b-ion sequences are correct."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        b_ions = peptide.generate_b_ions()
        
        assert b_ions[0]['fragment_sequence'] == 'P'
        assert b_ions[1]['fragment_sequence'] == 'PE'
        assert b_ions[2]['fragment_sequence'] == 'PEP'
        assert b_ions[6]['fragment_sequence'] == 'PEPTIDE'
    
    def test_b_ion_positions(self):
        """Test that b-ion positions are correct."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        b_ions = peptide.generate_b_ions()
        
        for i, ion in enumerate(b_ions):
            assert ion['fragment_position'] == i + 1
            assert ion['position'] == i + 1
    
    def test_b_ion_masses(self):
        """Test that b-ion masses are calculated correctly."""
        peptide = Peptide('A', use_cam=False)
        b_ions = peptide.generate_b_ions()
        
        # b1 for A: 71.037114 (Ala) + 1.007825 (proton)
        expected_mass = 71.037114 + 1.007825
        assert abs(b_ions[0]['fragment_mass'] - expected_mass) < 0.001
    
    def test_b_ion_type(self):
        """Test that b-ions have correct type."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        b_ions = peptide.generate_b_ions()
        
        for ion in b_ions:
            assert ion['fragment_type'] == 'b'


class TestYIonGeneration:
    """Test y-ion generation."""
    
    def test_simple_y_ions(self):
        """Test generating y-ions for a simple peptide."""
        peptide = Peptide('AAA', use_cam=False)
        y_ions = peptide.generate_y_ions(include_full_peptide=False)
        
        # Should have y1 and y2 (not Y0)
        assert len(y_ions) == 2
        # y-ions are generated from C to N, so y1 comes first
        assert y_ions[0]['fragment_name'] == 'y1'
        assert y_ions[1]['fragment_name'] == 'y2'
    
    def test_y_ion_with_full_peptide(self):
        """Test Y0 (full peptide) generation."""
        peptide = Peptide('AAA', use_cam=False)
        y_ions = peptide.generate_y_ions(include_full_peptide=True)
        
        # Should have y1, y2, Y0 (in that order - C to N)
        assert len(y_ions) == 3
        # First ion is y1, last is Y0 (full peptide)
        assert y_ions[-1]['fragment_name'] == 'Y0'
    
    def test_y_ion_sequences(self):
        """Test that y-ion sequences are correct."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        y_ions = peptide.generate_y_ions(include_full_peptide=False)
        
        # y-ions go from C-terminus: y1=E, y2=DE, etc.
        assert y_ions[0]['fragment_sequence'] == 'E'
        assert y_ions[1]['fragment_sequence'] == 'DE'
        assert y_ions[-1]['fragment_sequence'] == 'EPTIDE'
    
    def test_y_ion_positions(self):
        """Test that y-ion positions are from C-terminus."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        y_ions = peptide.generate_y_ions(include_full_peptide=True)
        
        # y1 has position 1 (first from C-terminus)
        assert y_ions[0]['fragment_position'] == 1
        # Y0 (full peptide) has position 7
        assert y_ions[-1]['fragment_position'] == 7
    
    def test_y_ion_includes_water(self):
        """Test that y-ions include H2O mass."""
        peptide = Peptide('A', use_cam=False)
        y_ions = peptide.generate_y_ions(include_full_peptide=True)
        
        # Y0 for A: 71.037114 (Ala) + 18.010565 (H2O) + 1.007825 (proton)
        expected_mass = 71.037114 + 18.010564684 + 1.007825
        assert abs(y_ions[0]['fragment_mass'] - expected_mass) < 0.001
    
    def test_y_ion_type(self):
        """Test that y-ions have correct type."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        y_ions = peptide.generate_y_ions()
        
        for ion in y_ions:
            assert ion['fragment_type'] == 'y'


class TestCIonGeneration:
    """Test c-ion generation."""
    
    def test_simple_c_ions(self):
        """Test generating c-ions for a simple peptide."""
        peptide = Peptide('AAA', use_cam=False)
        c_ions = peptide.generate_c_ions()
        
        assert len(c_ions) == 3
        assert c_ions[0]['fragment_name'] == 'c1'
        assert c_ions[1]['fragment_name'] == 'c2'
    
    def test_c_ion_type(self):
        """Test that c-ions have correct type."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        c_ions = peptide.generate_c_ions()
        
        for ion in c_ions:
            assert ion['fragment_type'] == 'c'
    
    def test_c_ion_includes_nh3(self):
        """Test that c-ions include NH3 mass."""
        peptide = Peptide('A', use_cam=False)
        c_ions = peptide.generate_c_ions()
        
        # c1 for A: 71.037114 (Ala) + 17.026549 (NH3) + 1.007825 (proton)
        expected_mass = 71.037114 + 17.026549 + 1.007825
        assert abs(c_ions[0]['fragment_mass'] - expected_mass) < 0.001


class TestZIonGeneration:
    """Test z-ion generation."""
    
    def test_simple_z_ions(self):
        """Test generating z-ions for a simple peptide."""
        peptide = Peptide('AAA', use_cam=False)
        z_ions = peptide.generate_z_ions()
        
        assert len(z_ions) == 3
        # z-ions generated from C to N: z1, z2, z3
        assert z_ions[0]['fragment_name'] == 'z1'
        assert z_ions[-1]['fragment_name'] == 'z3'
    
    def test_z_ion_type(self):
        """Test that z-ions have correct type."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        z_ions = peptide.generate_z_ions()
        
        for ion in z_ions:
            assert ion['fragment_type'] == 'z'


class TestAllFragments:
    """Test generation of all fragment types."""
    
    def test_generate_all_fragments(self):
        """Test generating all fragment types at once."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        fragments = peptide.generate_all_fragments(fragment_types=['by', 'cz'])
        
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        assert 'c_ions' in fragments
        assert 'z_ions' in fragments
        
        assert len(fragments['b_ions']) == 7
        assert len(fragments['y_ions']) == 7  # includes Y0
        assert len(fragments['c_ions']) == 7
        assert len(fragments['z_ions']) == 7
    
    def test_fragment_counts(self):
        """Test that fragment counts match peptide length."""
        peptide = Peptide('ABCDEF', use_cam=False)
        fragments = peptide.generate_all_fragments(fragment_types=['by', 'cz'])
        
        # All fragment types should have same count as peptide length
        assert len(fragments['b_ions']) == 6
        assert len(fragments['c_ions']) == 6
        assert len(fragments['z_ions']) == 6


class TestFragmentRetrieval:
    """Test retrieval of specific fragments."""
    
    def test_get_fragment_by_name_b(self):
        """Test retrieving a b-ion by name."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        fragment = peptide.get_fragment_by_name('b3')
        
        assert fragment is not None
        assert fragment['fragment_name'] == 'b3'
        assert fragment['fragment_sequence'] == 'PEP'
    
    def test_get_fragment_by_name_y(self):
        """Test retrieving a y-ion by name."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        fragment = peptide.get_fragment_by_name('y2')
        
        assert fragment is not None
        assert fragment['fragment_name'] == 'y2'
    
    def test_get_fragment_by_name_Y0(self):
        """Test retrieving Y0 (full peptide)."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        fragment = peptide.get_fragment_by_name('Y0')
        
        assert fragment is not None
        assert fragment['fragment_name'] == 'Y0'
        assert fragment['fragment_sequence'] == 'PEPTIDE'
    
    def test_get_fragment_invalid_name(self):
        """Test retrieving non-existent fragment."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        fragment = peptide.get_fragment_by_name('x99')
        
        assert fragment is None


class TestModifications:
    """Test peptide modifications."""
    
    def test_cam_modification(self):
        """Test that CAM modifies cysteine."""
        peptide_no_cam = Peptide('PEPTIDEC', use_cam=False)
        peptide_with_cam = Peptide('PEPTIDEC', use_cam=True)
        
        # Peptide with CAM should be heavier
        assert peptide_with_cam.mass > peptide_no_cam.mass
        
        # Difference should be approximately 57.021464 (CAM mass)
        mass_diff = peptide_with_cam.mass - peptide_no_cam.mass
        assert abs(mass_diff - 57.021464) < 0.1
    
    def test_modifications_applied_to_fragments(self):
        """Test that modifications affect fragment masses."""
        peptide = Peptide('PEC', use_cam=True)
        b_ions = peptide.generate_b_ions()
        
        # b3 should include CAM on cysteine
        b3 = b_ions[2]
        # Should be heavier than without CAM
        assert b3['fragment_mass'] > 300


class TestPeptideRepresentation:
    """Test string representation methods."""
    
    def test_repr(self):
        """Test __repr__ method."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        repr_str = repr(peptide)
        
        assert 'PEPTIDE' in repr_str
        assert 'Peptide' in repr_str
        assert 'mass=' in repr_str
    
    def test_len(self):
        """Test __len__ method."""
        peptide = Peptide('PEPTIDE')
        assert len(peptide) == 7
        
        peptide2 = Peptide('AAA')
        assert len(peptide2) == 3


class TestEdgeCases:
    """Test edge cases and special scenarios."""
    
    def test_single_amino_acid(self):
        """Test peptide with single amino acid."""
        peptide = Peptide('A', use_cam=False)
        
        b_ions = peptide.generate_b_ions()
        y_ions = peptide.generate_y_ions(include_full_peptide=True)
        
        assert len(b_ions) == 1
        assert len(y_ions) == 1
        assert b_ions[0]['fragment_name'] == 'b1'
        assert y_ions[0]['fragment_name'] == 'Y0'
    
    def test_two_amino_acids(self):
        """Test peptide with two amino acids."""
        peptide = Peptide('AA', use_cam=False)
        
        b_ions = peptide.generate_b_ions()
        y_ions = peptide.generate_y_ions(include_full_peptide=False)
        
        assert len(b_ions) == 2
        assert len(y_ions) == 1  # y1 only (no Y0)
    
    def test_all_same_amino_acid(self):
        """Test peptide with repeated amino acids."""
        peptide = Peptide('AAAAA', use_cam=False)
        fragments = peptide.generate_all_fragments()
        
        assert len(fragments['b_ions']) == 5
        assert len(fragments['y_ions']) == 5


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
