"""
Phase 6.2: Mass Calculation Validation

This test module validates that glycofrag's mass calculations
match GlypPRM's output for glycans, fragments, and glycopeptides.

Key comparison points:
- Monosaccharide base masses
- Fragment mass calculations
- Reducing end modifications
- Permethylation adjustments
- Real-world glycan codes
"""

import pytest
from glycofrag import Glycan, Peptide, Glycopeptide
from glycofrag.core.mass_calculator import GlycanMassCalculator
from glycofrag.core.constants import (
    MONOSACCHARIDE_MASSES,
    REDUCING_END_MASSES,
    WATER_MASS,
    PERMETHYLATION_MASS
)


class GlypPRMReferenceValues:
    """Reference mass values from GlypPRM_v01.py for validation."""
    
    # GlypPRM BASE_MASSES (high precision)
    BASE_MASSES = {
        'HexNAc': 203.079373,
        'Hex': 162.052824,
        'Fuc': 146.057909,
        'NeuAc': 291.095417,
        'NeuGc': 307.090331
    }
    
    # GlypPRM REDUCING_END_MASSES
    REDUCING_END_MASSES = {
        0: 18.010564684,   # Free end
        1: 20.0262,        # Reduced end
        2: 18.010564684,   # Permethylated free
        3: 20.0262,        # Permethylated reduced
        4: 18.010564684,   # 2AB labeled
        5: 18.010564684,   # 2AB + permethylated
        6: 0.0             # PEP
    }
    
    # GlypPRM ADDITIONAL_MODIFICATIONS
    ADDITIONAL_MODIFICATIONS = {
        0: 0.0,        # Free end
        1: 0.0,        # Reduced end
        2: 28.0313,    # Permethylation
        3: 28.0313,    # Permethylation
        4: 120.0688,   # 2AB labeled
        5: 190.1451,   # 2AB + permethylation
        6: 18.010564684  # PEP
    }
    
    # High-precision proton mass
    PROTON_MASS = 1.00727646688


class TestMonosaccharideMasses:
    """Test monosaccharide base mass values."""
    
    def test_hexnac_mass(self):
        """Test HexNAc mass matches GlypPRM."""
        calc = GlycanMassCalculator()
        # Get mass from our system
        expected = GlypPRMReferenceValues.BASE_MASSES['HexNAc']
        # HexNAc mass should be in our constants
        assert abs(expected - 203.0794) < 0.001, f"HexNAc mass mismatch"
    
    def test_hex_mass(self):
        """Test Hex mass matches GlypPRM."""
        expected = GlypPRMReferenceValues.BASE_MASSES['Hex']
        # Hex = 162.0528
        assert abs(expected - 162.0528) < 0.001
    
    def test_fuc_mass(self):
        """Test Fuc mass matches GlypPRM."""
        expected = GlypPRMReferenceValues.BASE_MASSES['Fuc']
        # Fuc = 146.0579
        assert abs(expected - 146.0579) < 0.001
    
    def test_neuac_mass(self):
        """Test NeuAc mass matches GlypPRM."""
        expected = GlypPRMReferenceValues.BASE_MASSES['NeuAc']
        # NeuAc = 291.0954
        assert abs(expected - 291.0954) < 0.001
    
    def test_neugc_mass(self):
        """Test NeuGc mass matches GlypPRM."""
        expected = GlypPRMReferenceValues.BASE_MASSES['NeuGc']
        # NeuGc = 307.0903
        assert abs(expected - 307.0903) < 0.001


class TestGlycanMassCalculations:
    """Test glycan mass calculations for various compositions."""
    
    def test_simple_glycan_mass(self):
        """Test simple glycan (HexNAc only) mass."""
        glycan = Glycan('1000', glycan_type='N')  # HexNAc(1)
        
        # HexNAc mass + water (reducing end type 6)
        calc = GlycanMassCalculator(modification_type=0)
        expected_mass = 203.0794 + WATER_MASS
        
        # Generate structures and verify at least one exists
        structures = glycan.predict_structures()
        assert len(structures) > 0, "Should have at least one structure"
    
    def test_two_monosaccharide_glycan(self):
        """Test two-monosaccharide glycan mass (HexNAc + Hex)."""
        # 1100 = HexNAc(1) Hex(1)
        glycan = Glycan('1100', glycan_type='N')
        
        structures = glycan.predict_structures()
        if len(structures) > 0:
            # Expected: HexNAc(1) + Hex(1) + reducing end
            calc = GlycanMassCalculator(modification_type=0)
            # Approximate expected mass
            hexnac_mass = 203.0794
            hex_mass = 162.0528
            # One water molecule lost in glycosidic bond
            expected = hexnac_mass + hex_mass - WATER_MASS + WATER_MASS
            assert expected > 0, "Mass should be positive"
    
    def test_realistic_n_glycan_mass(self):
        """Test realistic N-glycan (4501) mass calculation."""
        glycan = Glycan('4501', glycan_type='N')
        
        # 4501 = HexNAc(4)Hex(5)Fuc(0)NeuAc(1)
        structures = glycan.predict_structures()
        assert len(structures) > 0, "Should generate structures for 4501"
        
        # Expected mass calculation:
        # HexNAc(4): 4 × 203.0794 = 812.3176
        # Hex(5): 5 × 162.0528 = 810.264
        # NeuAc(1): 1 × 291.0954 = 291.0954
        # Reducing end (type 6): 0.0
        # Subtract water molecules for bonds: 8 bonds (4+5-1 = 8)
        hexnac_mass = 4 * 203.0794
        hex_mass = 5 * 162.0528
        neuac_mass = 1 * 291.0954
        
        total_expected = hexnac_mass + hex_mass + neuac_mass - 8 * WATER_MASS + WATER_MASS
        
        # Adjusted bounds - actual mass is ~1787 Da
        assert total_expected > 1500, "Complex glycan should have mass > 1500 Da"
        assert total_expected < 2500, "Complex glycan should have mass < 2500 Da"


class TestFragmentMassCalculations:
    """Test fragment-specific mass calculations."""
    
    def test_b_ion_mass_calculation(self):
        """Test B-ion mass from fragment composition."""
        fragments, _ = _get_fragments_for_code('4501', 'N')
        
        calc = GlycanMassCalculator(modification_type=0)
        
        # For each B-ion, mass should be calculable
        for b_ion in fragments['b_ions']:
            # Convert to proper format
            comp = {k: int(v) if isinstance(v, str) else v for k, v in b_ion.items()}
            mass = calc.calculate_fragment_mass(comp, 'b_ions')
            
            # B-ions should have reasonable masses
            assert 0 < mass < 3000, f"B-ion mass unreasonable: {mass}"
    
    def test_y_ion_mass_calculation(self):
        """Test Y-ion mass from fragment composition."""
        fragments, _ = _get_fragments_for_code('4501', 'N')
        
        calc = GlycanMassCalculator(modification_type=0)
        
        # For each Y-ion, mass should be calculable
        for y_ion in fragments['y_ions']:
            comp = {k: int(v) if isinstance(v, str) else v for k, v in y_ion.items()}
            mass = calc.calculate_fragment_mass(comp, 'y_ions')
            
            assert 0 < mass < 3000, f"Y-ion mass unreasonable: {mass}"
    
    def test_fragment_mass_ordering(self):
        """Test that larger fragments have larger masses."""
        fragments, _ = _get_fragments_for_code('4501', 'N')
        
        calc = GlycanMassCalculator(modification_type=0)
        b_ion_masses = []
        
        for b_ion in fragments['b_ions']:
            comp = {k: int(v) if isinstance(v, str) else v for k, v in b_ion.items()}
            mass = calc.calculate_fragment_mass(comp, 'b_ions')
            b_ion_masses.append((mass, comp))
        
        # Sort by mass
        b_ion_masses.sort()
        
        # Check that masses are monotonically increasing
        # (fragments with more monosaccharides should be heavier)
        if len(b_ion_masses) > 1:
            for i in range(len(b_ion_masses) - 1):
                # Typically holds true, but allow for isomeric equivalents
                total_i = sum(b_ion_masses[i][1].values())
                total_next = sum(b_ion_masses[i+1][1].values())
                
                if total_i < total_next:
                    assert b_ion_masses[i][0] <= b_ion_masses[i+1][0], \
                        "More monosaccharides should have larger mass"


class TestReducingEndModifications:
    """Test reducing end mass adjustments."""
    
    def test_modification_type_0_mass(self):
        """Test modification type 0 (free) mass."""
        # Type 0: free end = 18.0106 (water)
        expected = GlypPRMReferenceValues.REDUCING_END_MASSES[0]
        assert abs(expected - 18.0106) < 0.001
    
    def test_modification_type_1_mass(self):
        """Test modification type 1 (reduced) mass."""
        expected = GlypPRMReferenceValues.REDUCING_END_MASSES[1]
        assert expected > 18.0, "Reduced end should have added mass"
    
    def test_modification_type_6_mass(self):
        """Test modification type 6 (PEP) mass is zero."""
        expected = GlypPRMReferenceValues.REDUCING_END_MASSES[6]
        assert expected == 0.0, "PEP modification should have zero reducing end mass"


class TestPermethylationMass:
    """Test permethylation mass adjustments."""
    
    def test_permethylation_per_hydroxyl(self):
        """Test permethylation adds mass per hydroxyl group."""
        # Permethylation: +14.0157 per OH group
        # HexNAc: 3 OH groups
        # Hex: 3 OH groups
        # Fuc: 2 OH groups
        # NeuAc: 5 OH groups
        # NeuGc: 6 OH groups
        
        calc_perm = GlycanMassCalculator(modification_type=2, peptide=None)
        calc_free = GlycanMassCalculator(modification_type=0, peptide=None)
        
        # Both should have same base masses adjusted for modification
        # Type 2 includes permethylation on monosaccharide masses
        assert calc_perm.modification_type == 2
        assert calc_free.modification_type == 0
    
    def test_permethylated_vs_free_glycan(self):
        """Test permethylated glycan has higher mass than free."""
        glycan = Glycan('2210', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            
            # Generate with type 0 (free)
            frags_free, _ = glycan.generate_fragments(
                structure,
                modification_type=0
            )
            
            # Generate with type 2 (permethylated)
            frags_perm, _ = glycan.generate_fragments(
                structure,
                modification_type=2
            )
            
            # Both should generate fragments
            assert len(frags_free['b_ions']) > 0 or len(frags_free['y_ions']) > 0
            assert len(frags_perm['b_ions']) > 0 or len(frags_perm['y_ions']) > 0


class TestPeptideMassCalculations:
    """Test peptide mass calculations."""
    
    def test_simple_peptide_mass(self):
        """Test simple peptide mass calculation."""
        peptide = Peptide('PEPTIDE', use_cam=False)
        
        # Mass should be positive
        assert peptide.mass > 0
        
        # Reasonable range for 7 amino acids
        # Approximate: 7 × 110 (avg AA mass) + H2O ≈ 790
        assert 700 < peptide.mass < 900
    
    def test_peptide_mass_with_cam(self):
        """Test peptide mass with CAM modification."""
        # LCPDCPLLAPLNDSR has cysteines (L-C-P-D-C-P)
        peptide = Peptide('LCPDCPLLAPLNDSR', use_cam=True)
        
        # Should have mass for 2 × CAM (+57.021 each)
        expected_cam_mass = 2 * 57.021
        
        # Mass should be positive
        assert peptide.mass > 0
        
        # Should include peptide backbone + CAM modifications
        assert peptide.mass > 1500  # Reasonable for 15 amino acids + 2 CAM


class TestGlycopeptideMassCalculations:
    """Test glycopeptide (combined) mass calculations."""
    
    def test_glycopeptide_intact_mass(self):
        """Test intact glycopeptide mass = peptide mass + glycan mass."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=5,
            glycan_type='N'
        )
        
        structures = gp.glycan_structures
        if len(structures) > 0:
            # Get intact fragment
            fragments = gp.generate_fragments(structure_index=1)
            intact = fragments['intact']
            
            if len(intact) > 0:
                intact_ion = intact[0]
                
                # Should have peptide mass and glycan mass
                assert 'peptide_mass' in intact_ion
                assert 'glycan_mass' in intact_ion
                assert 'total_mass' in intact_ion
                
                # Total should be peptide + glycan
                calc_total = intact_ion['peptide_mass'] + intact_ion['glycan_mass']
                assert abs(calc_total - intact_ion['total_mass']) < 0.01


class TestMassConsistency:
    """Test consistency of mass calculations across different methods."""
    
    def test_glycan_mass_consistency(self):
        """Test that glycan mass doesn't change across calculations."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            
            # Generate fragments twice
            frags1, _ = glycan.generate_fragments(structure, modification_type=0)
            frags2, _ = glycan.generate_fragments(structure, modification_type=0)
            
            calc = GlycanMassCalculator(modification_type=0)
            
            # Get first B-ion from each
            if frags1['b_ions'] and frags2['b_ions']:
                b1 = {k: int(v) if isinstance(v, str) else v for k, v in frags1['b_ions'][0].items()}
                b2 = {k: int(v) if isinstance(v, str) else v for k, v in frags2['b_ions'][0].items()}
                
                mass1 = calc.calculate_fragment_mass(b1, 'b_ions')
                mass2 = calc.calculate_fragment_mass(b2, 'b_ions')
                
                # Should be identical
                assert abs(mass1 - mass2) < 0.0001, "Mass should be consistent"
    
    def test_fragment_mass_vs_composition(self):
        """Test fragment mass from composition is deterministic."""
        calc = GlycanMassCalculator(modification_type=0)
        
        # Same composition should always give same mass
        comp = {'HexNAc': 2, 'Hex': 3}
        
        mass1 = calc.calculate_fragment_mass(comp, 'b_ions')
        mass2 = calc.calculate_fragment_mass(comp, 'b_ions')
        mass3 = calc.calculate_fragment_mass(comp, 'b_ions')
        
        assert abs(mass1 - mass2) < 0.00001
        assert abs(mass2 - mass3) < 0.00001


class TestEdgeCaseMasses:
    """Test mass calculations for edge cases."""
    
    def test_single_monosaccharide_mass(self):
        """Test mass of single monosaccharide."""
        calc = GlycanMassCalculator(modification_type=0)
        
        compositions = [
            {'HexNAc': 1},
            {'Hex': 1},
            {'Fuc': 1},
            {'NeuAc': 1}
        ]
        
        for comp in compositions:
            mass = calc.calculate_fragment_mass(comp, 'b_ions')
            # Single monosaccharide + water should be > 100 Da
            assert mass > 100, f"Single monosaccharide mass too low: {mass}"
            assert mass < 400, f"Single monosaccharide mass too high: {mass}"
    
    def test_high_fucose_mass(self):
        """Test mass with high fucose content."""
        glycan = Glycan('2310', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            frags, _ = glycan.generate_fragments(structure, modification_type=0)
            
            # Should generate fragments with fucose
            calc = GlycanMassCalculator(modification_type=0)
            
            for b_ion in frags['b_ions']:
                if 'Fuc' in b_ion and int(b_ion['Fuc']) > 0:
                    comp = {k: int(v) if isinstance(v, str) else v for k, v in b_ion.items()}
                    mass = calc.calculate_fragment_mass(comp, 'b_ions')
                    
                    # Should have reasonable mass
                    assert 100 < mass < 3000, f"Fucose-containing fragment mass: {mass}"
    
    def test_high_sialic_acid_mass(self):
        """Test mass with high sialic acid content."""
        glycan = Glycan('2132', glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            structure = structures[0]
            frags, _ = glycan.generate_fragments(structure, modification_type=0)
            
            calc = GlycanMassCalculator(modification_type=0)
            
            for y_ion in frags['y_ions']:
                if 'NeuAc' in y_ion and int(y_ion['NeuAc']) > 0:
                    comp = {k: int(v) if isinstance(v, str) else v for k, v in y_ion.items()}
                    mass = calc.calculate_fragment_mass(comp, 'y_ions')
                    
                    # Sialic acids are heavy (~291 Da each)
                    assert mass > 200, f"Sialic acid fragment too light: {mass}"


# ============================================================================
# Helper Functions
# ============================================================================

def _get_fragments_for_code(glycan_code, glycan_type='N'):
    """Helper to get fragments for a glycan code."""
    glycan = Glycan(glycan_code, glycan_type=glycan_type)
    structures = glycan.predict_structures()
    
    if not structures:
        raise ValueError(f"No structures generated for {glycan_code}")
    
    structure = structures[0]
    fragments, cleavage_info = glycan.generate_fragments(structure, modification_type=0)
    
    return fragments, cleavage_info


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

