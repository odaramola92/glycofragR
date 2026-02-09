"""
Tests for Glycopeptide module (Phase 4).

This module tests combined peptide and glycan fragmentation including:
- Peptide fragments (b, y, c, z ions)
- Glycan fragments (B, Y, BY, YY, C, Z ions with -PEP suffix)
- Oxonium/diagnostic ions (A-type and custom ions)
- Intact glycopeptide ions
"""

import pytest
from glycofrag.glycopeptide import Glycopeptide


class TestGlycopeptideInitialization:
    """Test Glycopeptide initialization."""
    
    def test_init_simple_nglycopeptide(self):
        """Test initialization of simple N-glycopeptide."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,  # N in NxS/T motif
            glycan_type='N'
        )
        
        assert gp.peptide_sequence == 'EEQYNSTYR'
        assert gp.glycan_code == '4501'
        assert gp.glycosylation_site == 4
        assert gp.glycan_type == 'N'
    
    def test_init_oglycopeptide(self):
        """Test initialization of O-glycopeptide."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='2200',
            glycosylation_site=2,  # T or S
            glycan_type='O'
        )
        
        assert gp.glycan_type == 'O'
        assert gp.glycosylation_site == 2
    
    def test_invalid_glycosylation_site(self):
        """Test that invalid glycosylation site raises error."""
        with pytest.raises(ValueError):
            Glycopeptide(
                peptide_sequence='PEPTIDE',
                glycan_code='4501',
                glycosylation_site=10,  # Out of range
                glycan_type='N'
            )
    
    def test_negative_glycosylation_site(self):
        """Test that negative glycosylation site raises error."""
        with pytest.raises(ValueError):
            Glycopeptide(
                peptide_sequence='PEPTIDE',
                glycan_code='4501',
                glycosylation_site=-1,
                glycan_type='N'
            )
    
    def test_glycan_structures_predicted(self):
        """Test that glycan structures are predicted on init."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        assert len(gp.glycan_structures) > 0


class TestPeptideBackboneFragmentGeneration:
    """Test peptide backbone fragment generation (peptide_b_ions, peptide_y_ions)."""
    
    def test_backbone_generation(self):
        """Test that peptide backbone ions are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'peptide_y_ions' in fragments
        assert len(fragments['peptide_y_ions']) > 0
    
    def test_peptide_y_ions_are_backbone(self):
        """Test that peptide y-ions are backbone fragments."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['peptide_y_ions']:
            assert 'fragment_mass' in ion
            assert 'fragment_type' in ion
            assert ion['fragment_type'] == 'y'
    
    def test_peptide_y_ion_naming(self):
        """Test peptide y-ion naming convention."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['peptide_y_ions']:
            # y-ions are named y1, y2, ... or Y0
            assert ion['fragment_name'].startswith('y') or ion['fragment_name'] == 'Y0'
    
    def test_peptide_y_ion_mass(self):
        """Test that peptide y-ions have correct mass."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['peptide_y_ions']:
            assert 'fragment_mass' in ion
            assert ion['fragment_mass'] > 0


class TestPeptideBIonGeneration:
    """Test peptide b-ion generation."""
    
    def test_b_ion_generation(self):
        """Test that peptide b-ions are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'peptide_b_ions' in fragments
        assert len(fragments['peptide_b_ions']) > 0
    
    def test_b_ions_are_peptide_backbone(self):
        """Test that peptide b-ions are backbone fragments."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['peptide_b_ions']:
            assert ion['fragment_type'] == 'b'
    
    def test_b_ion_position_info(self):
        """Test that peptide b-ions have position info."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=5,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['peptide_b_ions']:
            assert 'fragment_position' in ion or 'position' in ion


class TestGlycanFragmentGeneration:
    """Test glycan fragment generation (B, Y, BY, YY, C, Z ions with -PEP)."""
    
    def test_glycan_generation(self):
        """Test that glycan ions with peptide are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['Y1'])
        
        assert 'y1_ions' in fragments
        assert len(fragments['y1_ions']) > 0
    
    def test_y1_has_core_glycan(self):
        """Test that Y1 ions have core glycan type."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['Y1'])
        
        for ion in fragments['y1_ions']:
            assert ion['type'] == 'Y1'
            assert ion['glycan_attached'] == True
            assert ion['glycan_type'] == 'core'
            assert 'glycan_composition' in ion
            assert 'glycan_mass' in ion
            assert 'total_mass' in ion
    
    def test_y1_mass_calculation(self):
        """Test Y1 ion mass = peptide + glycan core."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['Y1'])
        
        for ion in fragments['y1_ions']:
            # Y1 total_mass = mass field (already calculated with water loss)
            assert ion['total_mass'] > 0
            assert ion['total_mass'] == ion['mass']
    
    def test_y1_naming(self):
        """Test Y1 ion naming convention."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['Y1'])
        
        for ion in fragments['y1_ions']:
            assert ion['name'] == 'Y1'


class TestGlycanBYIonGeneration:
    """Test glycan BY-ion generation (glycan fragments with peptide)."""
    
    def test_glycan_by_generation(self):
        """Test that glycan BY ions are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'y_ions' in fragments
        assert len(fragments['y_ions']) > 0

    def test_glycan_b_generation(self):
        """Test that glycan B ions are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'b_ions' in fragments
        assert len(fragments['b_ions']) > 0
    
    def test_glycan_by_naming(self):
        """Test glycan BY-ion naming includes composition."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['by_ions']:
            assert 'composition' in ion
            assert len(ion['composition']) > 0


class TestOxoniumIonGeneration:
    """Test A-type oxonium and custom diagnostic ion generation."""
    
    def test_a_ions_generated(self):
        """Test that A-type oxonium ions are generated in glycopeptide."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'a_ions' in fragments
        assert len(fragments['a_ions']) > 0
    
    def test_custom_ions_generated(self):
        """Test that custom diagnostic ions are generated in glycopeptide."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'custom_ions' in fragments
        assert len(fragments['custom_ions']) > 0
    
    def test_a_ions_have_mass(self):
        """Test that A-type oxonium ions have correct mass."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        for ion in fragments['a_ions']:
            assert 'mass' in ion
            assert ion['mass'] > 0
    
    def test_standard_oxonium_ions_present(self):
        """Test that the 3 standard HexNAc oxonium ions are present."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        labels = [ion.get('_custom_label', '') for ion in fragments['a_ions']]
        assert 'HexNAc(-CH6O3)' in labels
        assert 'HexNAc(-C2H4O2)' in labels
        assert 'HexNAc-C2H6O3' in labels


class TestIntactGlycopeptide:
    """Test intact glycopeptide generation."""
    
    def test_intact_generation(self):
        """Test that intact ion is generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['peptide-b', 'peptide-y', 'Y1', 'intact'])
        
        assert 'intact' in fragments
        assert len(fragments['intact']) == 1
    
    def test_intact_properties(self):
        """Test intact ion has correct properties."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['peptide-b', 'peptide-y', 'Y1', 'intact'])
        intact = fragments['intact'][0]
        
        assert intact['type'] == 'intact'
        assert intact['name'] == 'Intact'
        assert 'peptide_mass' in intact
        assert 'glycan_mass' in intact
        assert 'total_mass' in intact
        assert 'sequence' in intact
        assert 'glycan_composition' in intact
    
    def test_intact_mass_calculation(self):
        """Test intact mass = peptide + glycan."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(fragment_types=['peptide-b', 'peptide-y', 'Y1', 'intact'])
        intact = fragments['intact'][0]
        
        expected_total = intact['peptide_mass'] + intact['glycan_mass']
        assert abs(intact['total_mass'] - expected_total) < 0.001


class TestMultipleFragmentTypes:
    """Test generating multiple fragment types together."""
    
    def test_all_fragment_types(self):
        """Test generating all fragment types."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(
            fragment_types=['peptide-b', 'peptide-y', 'Y1', 'intact']
        )
        
        # Should have all fragment types
        assert 'peptide_y_ions' in fragments
        assert 'peptide_b_ions' in fragments
        assert 'b_ions' in fragments
        assert 'y_ions' in fragments
        assert 'y1_ions' in fragments
        assert 'a_ions' in fragments
        assert 'custom_ions' in fragments
        assert 'intact' in fragments
        
        assert len(fragments['peptide_y_ions']) > 0
        assert len(fragments['peptide_b_ions']) > 0
        assert len(fragments['b_ions']) > 0
        assert len(fragments['y1_ions']) > 0
        assert len(fragments['a_ions']) > 0
        assert len(fragments['intact']) == 1
    
    def test_selective_fragment_types(self):
        """Test generating only specific fragment types."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        # Only Y1 (no peptide fragments)
        fragments = gp.generate_fragments(
            fragment_types=['Y1'],
            peptide_fragment_types=[],
            glycan_fragment_types=[]
        )
        
        assert len(fragments['y1_ions']) > 0
        assert len(fragments['peptide_b_ions']) == 0


class TestModificationSupport:
    """Test support for different modifications."""
    
    def test_cam_modification(self):
        """Test CAM modification support."""
        gp = Glycopeptide(
            peptide_sequence='PEPTCIDE',  # Has cysteine
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N',
            use_cam=True
        )
        
        assert gp.use_cam == True
        fragments = gp.generate_fragments()
        assert len(fragments['peptide_y_ions']) > 0
    
    def test_default_glycopeptide_mode(self):
        """Test that glycopeptide always uses glycopeptide modification mode."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='3300',
            glycosylation_site=4,
            glycan_type='N'
        )
        fragments = gp.generate_fragments()
        assert len(fragments['y1_ions']) > 0
    
    def test_modification_mode_is_glycopeptide(self):
        """Test that the modification type is always glycopeptide mode (6)."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        assert gp.modification_type == 6
        frags = gp.generate_fragments()
        assert len(frags['y1_ions']) > 0


class TestGlycanStructureSelection:
    """Test selecting different predicted glycan structures."""
    
    def test_structure_index_1(self):
        """Test using first predicted structure."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments(structure_index=1)
        assert len(fragments['y1_ions']) > 0
    
    def test_invalid_structure_index(self):
        """Test that invalid structure index raises error."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='3300',
            glycosylation_site=2,
            glycan_type='N'
        )
        
        with pytest.raises(ValueError):
            gp.generate_fragments(structure_index=999)


class TestGlycopeptideRepresentation:
    """Test string representation."""
    
    def test_repr(self):
        """Test __repr__ method."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=5,
            glycan_type='N'
        )
        
        repr_str = repr(gp)
        assert 'EEQYNSTYR' in repr_str
        assert '4501' in repr_str
        assert 'N5' in repr_str


class TestOGlycopeptideFragmentation:
    """Test O-glycopeptide specific cases."""
    
    def test_oglycopeptide_fragments(self):
        """Test O-glycopeptide fragmentation."""
        gp = Glycopeptide(
            peptide_sequence='PEPTIDE',
            glycan_code='2200',
            glycosylation_site=2,
            glycan_type='O'
        )
        
        fragments = gp.generate_fragments(
            fragment_types=['peptide-b', 'peptide-y', 'Y1', 'intact']
        )
        
        assert len(fragments['peptide_y_ions']) > 0
        assert len(fragments['y1_ions']) > 0
        assert len(fragments['intact']) >= 1  # One per structure


class TestNeutralLossGeneration:
    """Test peptide neutral loss generation in glycopeptide context."""
    
    def test_neutral_losses_generated(self):
        """Test that peptide neutral losses are generated."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        assert 'peptide_neutral_losses' in fragments
        assert len(fragments['peptide_neutral_losses']) > 0
    
    def test_neutral_loss_types(self):
        """Test that expected neutral loss types are present."""
        gp = Glycopeptide(
            peptide_sequence='EEQYNSTYR',
            glycan_code='4501',
            glycosylation_site=4,
            glycan_type='N'
        )
        
        fragments = gp.generate_fragments()
        
        loss_types = set()
        for nl in fragments['peptide_neutral_losses']:
            if 'type' in nl:
                loss_types.add(nl['type'])
        
        assert len(loss_types) > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
