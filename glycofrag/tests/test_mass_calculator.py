"""
Tests for mass calculator module.
"""

import pytest
from glycofrag.core.mass_calculator import (
    GlycanMassCalculator,
    calculate_glycan_mass,
    calculate_peptide_mass,
    calculate_mz,
)


class TestGlycanCodeParsing:
    """Tests for glycan code parsing."""
    
    def test_parse_numeric_code_4_digit(self):
        """Test parsing 4-digit numeric code."""
        calc = GlycanMassCalculator()
        result = calc.parse_glycan_code("4501")
        assert result == (4, 5, 0, 1, 0)
    
    def test_parse_numeric_code_5_digit(self):
        """Test parsing 5-digit numeric code."""
        calc = GlycanMassCalculator()
        result = calc.parse_glycan_code("45012")
        assert result == (4, 5, 0, 1, 2)
    
    def test_parse_named_code(self):
        """Test parsing named glycan code."""
        calc = GlycanMassCalculator()
        result = calc.parse_glycan_code("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)")
        assert result == (4, 5, 1, 2, 0)
    
    def test_parse_named_code_partial(self):
        """Test parsing partial named glycan code."""
        calc = GlycanMassCalculator()
        result = calc.parse_glycan_code("HexNAc(2)Hex(3)")
        assert result == (2, 3, 0, 0, 0)
    
    def test_parse_code_with_spaces(self):
        """Test parsing code with spaces."""
        calc = GlycanMassCalculator()
        result = calc.parse_glycan_code(" 4501 ")
        assert result == (4, 5, 0, 1, 0)


class TestGlycanMassCalculation:
    """Tests for glycan mass calculation."""
    
    def test_simple_glycan_mass(self):
        """Test mass calculation for a simple glycan."""
        calc = GlycanMassCalculator(modification_type=0)
        mass = calc.calculate_glycan_mass("2200")  # HexNAc2 Hex2
        
        # Expected: 2*203.079 + 2*162.053 + 18.011 (water) = 748.275
        expected = 2 * 203.079373 + 2 * 162.052824 + 18.010564684
        assert abs(mass - expected) < 0.01
    
    def test_complex_glycan_mass(self):
        """Test mass calculation for a complex glycan."""
        calc = GlycanMassCalculator(modification_type=0)
        mass = calc.calculate_glycan_mass("4501")
        
        # Expected: 4*HexNAc + 5*Hex + 0*Fuc + 1*NeuAc + water
        expected = (4 * 203.079373 + 5 * 162.052824 + 
                   0 * 146.057909 + 1 * 291.095417 + 18.010564684)
        assert abs(mass - expected) < 0.01
    
    def test_glycan_with_fucose(self):
        """Test mass calculation for glycan with fucose."""
        calc = GlycanMassCalculator(modification_type=0)
        mass = calc.calculate_glycan_mass("4511")
        
        expected = (4 * 203.079373 + 5 * 162.052824 + 
                   1 * 146.057909 + 1 * 291.095417 + 18.010564684)
        assert abs(mass - expected) < 0.01


class TestPeptideMassCalculation:
    """Tests for peptide mass calculation."""
    
    def test_simple_peptide_mass(self):
        """Test mass calculation for simple peptide."""
        calc = GlycanMassCalculator(use_cam=False)
        mass = calc.calculate_peptide_mass("GGG")  # Three glycines
        
        # Expected: 3 * Gly + H2O = 3 * 57.021464 + 18.010564684
        expected = 3 * 57.021464 + 18.010564684
        assert abs(mass - expected) < 0.01
    
    def test_peptide_with_cam(self):
        """Test peptide mass with CAM on cysteine."""
        calc = GlycanMassCalculator(use_cam=True)
        mass = calc.calculate_peptide_mass("GCG")  # Gly-Cys-Gly
        
        # Expected: 2*Gly + Cys + CAM + H2O
        expected = 2 * 57.021464 + 103.009185 + 57.021464 + 18.010564684
        assert abs(mass - expected) < 0.01
    
    def test_peptide_without_cam(self):
        """Test peptide mass without CAM."""
        calc = GlycanMassCalculator(use_cam=False)
        mass = calc.calculate_peptide_mass("GCG")
        
        # Expected: 2*Gly + Cys + H2O (no CAM)
        expected = 2 * 57.021464 + 103.009185 + 18.010564684
        assert abs(mass - expected) < 0.01
    
    def test_peptide_with_oxidation(self):
        """Test peptide mass with methionine oxidation."""
        calc = GlycanMassCalculator(use_cam=False, variable_mods=["Ox:M"])
        mass = calc.calculate_peptide_mass("GMG")
        
        # Expected: 2*Gly + Met + Ox + H2O
        expected = 2 * 57.021464 + 131.040485 + 15.994915 + 18.010564684
        assert abs(mass - expected) < 0.01
    
    def test_realistic_peptide(self):
        """Test a realistic peptide sequence."""
        calc = GlycanMassCalculator(use_cam=True)
        mass = calc.calculate_peptide_mass("LCPDCPLLAPLNDSR")
        
        # This peptide has 2 cysteines, so CAM should add 2 * 57.021
        # Just verify mass is in reasonable range
        assert mass > 1500  # Should be > 1500 Da
        assert mass < 2000  # Should be < 2000 Da


class TestMzCalculation:
    """Tests for m/z calculation."""
    
    def test_mz_charge_1(self):
        """Test m/z at charge state 1."""
        calc = GlycanMassCalculator()
        mz = calc.calculate_mz(1000.0, 1)
        
        # Expected: (1000 + 1*1.00728) / 1 = 1001.00728
        expected = 1000.0 + 1.00727646688
        assert abs(mz - expected) < 0.0001
    
    def test_mz_charge_2(self):
        """Test m/z at charge state 2."""
        calc = GlycanMassCalculator()
        mz = calc.calculate_mz(1000.0, 2)
        
        # Expected: (1000 + 2*1.00728) / 2 = 501.00728
        expected = (1000.0 + 2 * 1.00727646688) / 2
        assert abs(mz - expected) < 0.0001
    
    def test_mz_charge_3(self):
        """Test m/z at charge state 3."""
        calc = GlycanMassCalculator()
        mz = calc.calculate_mz(1000.0, 3)
        
        expected = (1000.0 + 3 * 1.00727646688) / 3
        assert abs(mz - expected) < 0.0001
    
    def test_mz_invalid_charge(self):
        """Test m/z with invalid charge raises error."""
        calc = GlycanMassCalculator()
        with pytest.raises(ValueError):
            calc.calculate_mz(1000.0, 0)
        with pytest.raises(ValueError):
            calc.calculate_mz(1000.0, -1)


class TestFragmentMassCalculation:
    """Tests for fragment mass calculation."""
    
    def test_b_ion_mass(self):
        """Test B-ion fragment mass calculation."""
        calc = GlycanMassCalculator(modification_type=0)
        composition = {'HexNAc': 1, 'Hex': 1}
        mass = calc.calculate_fragment_mass(composition, 'b_ions')
        
        # B-ion: just the glycan residues
        expected = 203.079373 + 162.052824
        assert abs(mass - expected) < 0.01
    
    def test_y_ion_mass_glycopeptide(self):
        """Test Y-ion fragment mass in glycopeptide mode."""
        calc = GlycanMassCalculator(
            modification_type=6,  # Internal calculator can use mode 6
            use_cam=False,
            peptide="GGG"
        )
        composition = {'HexNAc': 1}
        mass = calc.calculate_fragment_mass(composition, 'y_ions', peptide="GGG")
        
        # Y-ion in glycopeptide mode: glycan + peptide
        expected = 374.14376500000003
        assert abs(mass - expected) < 0.1


class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    def test_calculate_glycan_mass_function(self):
        """Test standalone glycan mass function."""
        mass = calculate_glycan_mass("2200")
        assert mass > 0
    
    def test_calculate_peptide_mass_function(self):
        """Test standalone peptide mass function."""
        mass = calculate_peptide_mass("PEPTIDE")
        assert mass > 0
    
    def test_calculate_mz_function(self):
        """Test standalone m/z function."""
        mz = calculate_mz(1000.0, 2)
        assert abs(mz - 501.00728) < 0.01


class TestCaching:
    """Tests for caching behavior."""
    
    def test_peptide_mass_caching(self):
        """Test that peptide mass calculation is cached."""
        calc = GlycanMassCalculator(use_cam=True)
        
        # First calculation
        mass1 = calc.calculate_peptide_mass("PEPTIDE")
        
        # Second calculation (should use cache)
        mass2 = calc.calculate_peptide_mass("PEPTIDE")
        
        assert mass1 == mass2
    
    def test_cache_different_mods(self):
        """Test that different modifications produce different cache entries."""
        calc = GlycanMassCalculator()
        
        mass_cam = calc.calculate_peptide_mass("GCG", use_cam=True)
        mass_no_cam = calc.calculate_peptide_mass("GCG", use_cam=False)
        
        # Should be different (CAM adds ~57 Da)
        assert abs(mass_cam - mass_no_cam) > 50


class TestGlycopeptideMass:
    """Tests for glycopeptide mass calculation."""
    
    def test_glycopeptide_mass(self):
        """Test glycopeptide mass calculation."""
        calc = GlycanMassCalculator(use_cam=False)
        mass = calc.calculate_glycopeptide_mass("GGG", "2200")
        
        # Should be peptide + glycan
        peptide_mass = calc.calculate_peptide_mass("GGG")
        
        # Glycan residue mass (no reducing end for glycopeptide)
        glycan_residue = 2 * 203.079373 + 2 * 162.052824
        
        expected = peptide_mass + glycan_residue
        assert abs(mass - expected) < 0.1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

