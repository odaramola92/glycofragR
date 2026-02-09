"""
Tests for core constants module.
"""

import pytest
from glycofrag.core.constants import (
    MONOSACCHARIDE_MASSES,
    AMINO_ACID_MASSES,
    MODIFICATION_MASSES,
    MODIFICATION_TARGETS,
    PROTON_MASS,
    WATER_MASS,
    REDUCING_END_MASSES,
    OXONIUM_IONS,
)


class TestMonosaccharideMasses:
    """Tests for monosaccharide mass constants."""
    
    def test_hexnac_mass(self):
        """Test HexNAc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['HexNAc'] - 203.079373) < 0.0001
    
    def test_hex_mass(self):
        """Test Hex mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['Hex'] - 162.052824) < 0.0001
    
    def test_fuc_mass(self):
        """Test Fuc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['Fuc'] - 146.057909) < 0.0001
    
    def test_neuac_mass(self):
        """Test NeuAc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['NeuAc'] - 291.095417) < 0.0001
    
    def test_neugc_mass(self):
        """Test NeuGc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['NeuGc'] - 307.090331) < 0.0001
    
    def test_all_monosaccharides_present(self):
        """Test all expected monosaccharides are defined."""
        expected = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc']
        for mono in expected:
            assert mono in MONOSACCHARIDE_MASSES


class TestAminoAcidMasses:
    """Tests for amino acid mass constants."""
    
    def test_all_standard_amino_acids_present(self):
        """Test all 20 standard amino acids are defined."""
        expected = 'ACDEFGHIKLMNPQRSTVWY'
        for aa in expected:
            assert aa in AMINO_ACID_MASSES
    
    def test_glycine_mass(self):
        """Test glycine (smallest AA) mass is correct."""
        assert abs(AMINO_ACID_MASSES['G'] - 57.021464) < 0.0001
    
    def test_tryptophan_mass(self):
        """Test tryptophan (largest AA) mass is correct."""
        assert abs(AMINO_ACID_MASSES['W'] - 186.079313) < 0.0001
    
    def test_cysteine_mass(self):
        """Test cysteine mass (important for CAM)."""
        assert abs(AMINO_ACID_MASSES['C'] - 103.009185) < 0.0001


class TestModificationMasses:
    """Tests for modification mass constants."""
    
    def test_cam_mass(self):
        """Test carbamidomethylation mass is correct."""
        assert abs(MODIFICATION_MASSES['CAM'] - 57.021464) < 0.0001
    
    def test_oxidation_mass(self):
        """Test oxidation mass is correct."""
        assert abs(MODIFICATION_MASSES['Ox'] - 15.994915) < 0.0001
    
    def test_phosphorylation_mass(self):
        """Test phosphorylation mass is correct."""
        assert abs(MODIFICATION_MASSES['Phos'] - 79.966331) < 0.0001
    
    def test_deamidation_mass(self):
        """Test deamidation mass is correct."""
        assert abs(MODIFICATION_MASSES['Deam'] - 0.984016) < 0.0001


class TestModificationTargets:
    """Tests for modification target definitions."""
    
    def test_cam_targets_cysteine(self):
        """Test CAM targets cysteine."""
        assert 'C' in MODIFICATION_TARGETS['CAM']
    
    def test_oxidation_targets_methionine(self):
        """Test oxidation targets methionine."""
        assert 'M' in MODIFICATION_TARGETS['Ox']
    
    def test_phosphorylation_targets(self):
        """Test phosphorylation targets STY."""
        targets = MODIFICATION_TARGETS['Phos']
        assert 'S' in targets
        assert 'T' in targets
        assert 'Y' in targets


class TestPhysicalConstants:
    """Tests for physical constants."""
    
    def test_proton_mass(self):
        """Test proton mass is correct (high precision)."""
        assert abs(PROTON_MASS - 1.00727646688) < 1e-9
    
    def test_water_mass(self):
        """Test water mass is correct."""
        assert abs(WATER_MASS - 18.010564684) < 0.0001


class TestReducingEndMasses:
    """Tests for reducing end mass constants."""
    
    def test_free_reducing_end(self):
        """Test free reducing end mass (type 0)."""
        assert abs(REDUCING_END_MASSES[0] - 18.010564684) < 0.0001
    
    def test_glycopeptide_mode(self):
        """Test glycopeptide mode has zero reducing end mass (type 6)."""
        assert REDUCING_END_MASSES[6] == 0.0


class TestOxoniumIons:
    """Tests for oxonium ion masses."""
    
    def test_hexnac_oxonium(self):
        """Test HexNAc oxonium ion mass."""
        assert abs(OXONIUM_IONS['HexNAc'] - 204.086646) < 0.0001
    
    def test_neuac_oxonium(self):
        """Test NeuAc oxonium ion mass."""
        assert abs(OXONIUM_IONS['NeuAc'] - 292.102690) < 0.0001


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
