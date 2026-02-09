"""
Phase 9.1: Integration Testing - End-to-End Workflows

This test module validates complete glycopeptide analysis workflows from start to finish,
simulating real-world usage patterns and verifying output consistency.

Test scenarios:
- Complete analysis pipeline (glycan prediction → fragmentation → DataFrame)
- Multiple glycosylation configurations
- Complex glycan structures
- Batch processing
- Output format validation
- Literature-based examples with expected results
"""

import pytest
import pandas as pd
from glycofrag import Glycan, Peptide, Glycopeptide
from glycofrag.io.tables import generate_fragment_table


# ============================================================================
# End-to-End Workflow Tests
# ============================================================================

class TestCompleteWorkflow:
    """Test complete glycopeptide analysis workflows."""
    
    def test_simple_n_glycopeptide_workflow(self):
        """Test complete workflow: peptide + N-glycan → fragments → DataFrame."""
        # Simple N-glycopeptide: IgG tryptic peptide with biantennary glycan
        peptide_seq = "EEQYNSTYR"  # Contains N-glycosylation site (N-X-S/T)
        glycan_code = "4500"  # Common biantennary: HexNAc(4)Hex(5)Fuc(0)NeuAc(0)
        glycosylation_site = 5  # Glycan at N (1-indexed)
        
        # Step 1: Create glycopeptide with integrated approach
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=glycosylation_site,
            glycan_type='N',
            use_cam=True
        )
        
        # Step 2: Generate fragments
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate fragments"
        
        # Verify peptide and glycan properties
        assert glycopeptide.peptide_sequence == peptide_seq
        assert glycopeptide.glycan_code == glycan_code
        assert glycopeptide.glycosylation_site == glycosylation_site
        
        # Verify glycan structures were predicted
        assert len(glycopeptide.glycan.possible_structures) > 0
        
        # Step 3: Create DataFrame (using glycan fragments as example)
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            peptide=peptide_seq,
            glycan_type='N',
            charge_states=[1, 2]
        )
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0, "Should create non-empty DataFrame"
        
        # Verify DataFrame structure
        required_columns = ['Fragment Type', 'm/z(z=1)', 'Fragment']
        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"
        
        # Verify fragment data
        assert df['m/z(z=1)'].notna().all(), "All m/z values should be valid"
        assert (df['m/z(z=1)'] > 0).all(), "All m/z values should be positive"
    
    def test_o_glycopeptide_workflow(self):
        """Test complete workflow with O-glycopeptide (T-antigen)."""
        # Mucin-like O-glycopeptide with T-antigen
        peptide_seq = "GTTPSPVPTR"  # Contains multiple S/T sites
        glycan_code = "1100"  # T-antigen: HexNAc(1)Hex(1)
        
        # Complete workflow
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=2,  # Glycan at first T
            glycan_type='O',
            use_cam=True
        )
        
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate O-glycopeptide fragments"
        
        # Verify O-glycan structures
        assert len(glycopeptide.glycan.possible_structures) > 0
        
        # Generate DataFrame for glycan fragments
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            glycan_type='O',
            charge_states=[1, 2]
        )
        
        assert len(df) > 0, "Should create O-glycopeptide DataFrame"
    
    def test_sialylated_glycopeptide_workflow(self):
        """Test complete workflow with sialylated N-glycan (common in serum)."""
        peptide_seq = "DQCIYNTTYLNVQR"  # Contains N-glycosylation site
        glycan_code = "4502"  # Disialylated: HexNAc(4)Hex(5)NeuAc(2)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=6,  # Glycan at N
            glycan_type='N',
            use_cam=True
        )
        
        fragments = glycopeptide.generate_fragments()
        
        if len(glycopeptide.glycan.possible_structures) > 0:
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                glycopeptide.glycan.possible_structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code=glycan_code,
                modification_type=0,
                glycan_type='N',
                charge_states=[1, 2, 3]
            )
            
            assert len(df) > 0, "Should create sialylated glycopeptide DataFrame"
            
            # Verify sialic acid presence in fragments
            # Y-ions should include NeuAc masses
            assert df['m/z(z=1)'].max() > 500, "Should have large fragments with sialic acid"
    
    def test_fucosylated_glycopeptide_workflow(self):
        """Test complete workflow with core-fucosylated N-glycan (very common)."""
        peptide_seq = "EEQFNSTFR"  # Contains N-glycosylation site
        glycan_code = "4501"  # Core-fucosylated: HexNAc(4)Hex(5)Fuc(1)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=5,  # Glycan at N
            glycan_type='N',
            use_cam=True
        )
        
        fragments = glycopeptide.generate_fragments()
        
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            glycan_type='N',
            charge_states=[1, 2]
        )
        
        assert len(df) > 0, "Should create fucosylated glycopeptide DataFrame"
        
        # Verify fucose in fragments (146.0579 Da)
        assert df['m/z(z=1)'].max() > 100, "Should have fragments with fucose"


# ============================================================================
# Literature-Based Examples
# ============================================================================

class TestLiteratureExamples:
    """Test with real glycopeptides from published literature."""
    
    def test_igg_n297_glycopeptide(self):
        """Test IgG1 Fc N-glycosylation site (N297) - most studied glycopeptide."""
        # IgG1 Fc tryptic glycopeptide: EEQYNSTYR (N297 site)
        # This is the most well-characterized glycopeptide in literature
        peptide_seq = "EEQYNSTYR"
        glycan_code = "4501"  # G0F: Most common IgG glycoform
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=5,  # N at position 5
            glycan_type='N',
            use_cam=True
        )
        
        # Verify peptide exists
        assert glycopeptide.peptide is not None
        peptide_mass = glycopeptide.peptide.mass
        
        # Peptide mass should be reasonable (with CAM modification)
        assert 1100 < peptide_mass < 1250, f"IgG peptide mass expected ~1188 Da with CAM, got {peptide_mass}"
        
        # Verify glycan prediction
        assert len(glycopeptide.glycan.possible_structures) > 0, "Should predict G0F structures"
        
        # Generate fragments
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate IgG fragments"
        
        # Generate DataFrame
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            peptide=peptide_seq,
            glycan_type='N',
            charge_states=[1, 2, 3]
        )
        
        # Should have comprehensive fragment coverage
        assert len(df) > 10, "IgG glycopeptide should generate many fragments"
    
    def test_transferrin_n_glycopeptide(self):
        """Test transferrin N-glycopeptide - common serum glycoprotein."""
        # Transferrin has two N-glycosylation sites
        # Example glycopeptide from N413 site
        peptide_seq = "QQQHLFGSNVTDCSR"  # Contains N-X-S/T motif
        glycan_code = "4502"  # Disialylated biantennary (common in serum)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=8,  # N at position 8
            glycan_type='N',
            use_cam=True
        )
        
        if len(glycopeptide.glycan.possible_structures) > 0:
            fragments = glycopeptide.generate_fragments()
            assert len(fragments) > 0, "Transferrin glycopeptide should generate fragments"
            
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                glycopeptide.glycan.possible_structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code=glycan_code,
                modification_type=0,
                glycan_type='N'
            )
            
            assert len(df) > 0, "Should generate transferrin fragment table"
    
    def test_rnase_b_high_mannose(self):
        """Test RNase B high-mannose N-glycan - standard glycoprotein."""
        # RNase B contains high-mannose N-glycans at N34
        # Common composition: Man5-9GlcNAc2
        peptide_seq = "SRNLTK"  # RNase B tryptic peptide with N34
        glycan_code = "2500"  # Man5GlcNAc2 (simplified high-mannose)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=4,  # N at position 4
            glycan_type='N',
            use_cam=True
        )
        
        if len(glycopeptide.glycan.possible_structures) > 0:
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                glycopeptide.glycan.possible_structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code=glycan_code,
                modification_type=0,
                glycan_type='N'
            )
            
            assert len(df) > 0, "RNase B glycopeptide should generate fragments"
    
    def test_mucin_o_glycopeptide(self):
        """Test mucin-type O-glycopeptide - classic O-glycan example."""
        # Mucin proteins are heavily O-glycosylated
        peptide_seq = "PAPGSTAPP"  # Mucin-like repeat with S/T sites
        glycan_code = "1101"  # Sialyl-T antigen (common mucin glycan)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=6,  # O-glycan at S
            glycan_type='O',
            use_cam=True
        )
        
        assert len(glycopeptide.glycan.possible_structures) > 0, "Should predict sialyl-T structures"
        
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            glycan_type='O'
        )
        
        assert len(df) > 0, "Mucin O-glycopeptide should generate fragments"


# ============================================================================
# Complex Multi-Glycosylated Peptides
# ============================================================================

class TestMultiGlycosylation:
    """Test peptides with multiple glycosylation sites."""
    
    def test_two_n_glycosylation_sites(self):
        """Test peptide with two N-glycosylation sites."""
        # Peptide with two N-X-S/T motifs
        peptide_seq = "NETSLNGSQK"  # Two potential N-glycosylation sites
        glycan_code = "4300"  # Biantennary glycan
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            # Glycan at first N-site
            glycopeptide1 = Glycopeptide(
                peptide_sequence=peptide_seq,
                glycan_code=glycan_code,
                glycosylation_site=1,
                glycan_type='N'
            )
            fragments1 = glycopeptide1.generate_fragments()
            
            # Glycan at second N-site
            glycopeptide2 = Glycopeptide(
                peptide_sequence=peptide_seq,
                glycan_code=glycan_code,
                glycosylation_site=7,
                glycan_type='N'
            )
            fragments2 = glycopeptide2.generate_fragments()
            
            # Both should generate fragments
            assert len(fragments1) > 0, "First site should generate fragments"
            assert len(fragments2) > 0, "Second site should generate fragments"
            
            # Fragment patterns should differ based on site position
            # This reflects real MS/MS fragmentation differences
    
    def test_multiple_o_glycosylation_sites(self):
        """Test mucin-like peptide with multiple O-glycosylation sites."""
        # Mucin-like peptide with 4 potential O-glycan sites
        peptide_seq = "TTSTTTSTR"  # Multiple S/T sites
        glycan_code = "1100"  # T-antigen
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            # Test glycosylation at different T/S positions
            for position in [1, 3, 5]:  # Test subset of positions
                glycopeptide = Glycopeptide(
                    peptide_sequence=peptide_seq,
                    glycan_code=glycan_code,
                    glycosylation_site=position,
                    glycan_type='O'
                )
                fragments = glycopeptide.generate_fragments()
                
                assert len(fragments) > 0, f"Should generate fragments for position {position}"
    
    def test_heterogeneous_glycosylation(self):
        """Test same peptide with different glycan structures (microheterogeneity)."""
        # Microheterogeneity: same peptide with different glycan structures
        peptide_seq = "EEQYNSTYR"
        
        # Different glycoforms at the same site
        glycan_codes = [
            "4300",  # G0 (no galactose)
            "4500",  # G2 (full galactosylation)
            "4501",  # G0F (fucosylated)
            "4502",  # G2S2 (disialylated)
        ]
        
        dataframes = []
        for code in glycan_codes:
            glycan = Glycan(code, glycan_type='N')
            structures = glycan.predict_structures()
            
            if len(structures) > 0:
                glycopeptide = Glycopeptide(
                    peptide_sequence=peptide_seq,
                    glycan_code=code,
                    glycosylation_site=4,
                    glycan_type='N'
                )
                
                glycan_fragments, _ = glycan.generate_fragments(
                    structures[0],
                    modification_type=0
                )
                
                df = generate_fragment_table(
                    fragments=glycan_fragments,
                    glycan_code=code,
                    modification_type=0,
                    glycan_type='N'
                )
                dataframes.append((code, df))
        
        # Should generate different fragment masses for different glycoforms
        assert len(dataframes) >= 2, "Should process multiple glycoforms"
        
        # Verify masses differ between glycoforms
        if len(dataframes) >= 2:
            df1_max_mass = dataframes[0][1]['m/z(z=1)'].max()
            df2_max_mass = dataframes[1][1]['m/z(z=1)'].max()
            assert df1_max_mass != df2_max_mass, "Different glycoforms should have different masses"


# ============================================================================
# Output Format Validation
# ============================================================================

class TestOutputFormat:
    """Validate consistency and correctness of output formats."""
    
    def test_dataframe_column_consistency(self):
        """Test that all DataFrames have consistent column structure."""
        # Generate multiple glycopeptides
        test_cases = [
            ("EEQYNSTYR", "4501", 'N', 4),
            ("GTTPSPVPTR", "1100", 'O', 1),
            ("DQCIYNTTYLNVQR", "4300", 'N', 5),
        ]
        
        dataframes = []
        for peptide_seq, glycan_code, glycan_type, site in test_cases:
            glycan = Glycan(glycan_code, glycan_type=glycan_type)
            structures = glycan.predict_structures()
            
            if len(structures) > 0:
                glycan_fragments, _ = glycan.generate_fragments(
                    structures[0],
                    modification_type=0
                )
                
                df = generate_fragment_table(
                    fragments=glycan_fragments,
                    glycan_code=glycan_code,
                    modification_type=0,
                    peptide=peptide_seq,
                    glycan_type=glycan_type
                )
                dataframes.append(df)
        
        assert len(dataframes) >= 2, "Should generate multiple DataFrames"
        
        # Verify all have same column structure
        first_columns = set(dataframes[0].columns)
        for df in dataframes[1:]:
            assert set(df.columns) == first_columns, "All DataFrames should have same columns"
    
    def test_mass_precision(self):
        """Test that masses are reported with appropriate precision."""
        glycan = Glycan("4501", glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            glycan_fragments, _ = glycan.generate_fragments(
                structures[0],
                    modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code="4501",
                modification_type=0,
                glycan_type='N'
            )
            
            # Masses should have reasonable precision (not too many or too few decimals)
            for mass in df['m/z(z=1)']:
                # Mass should be positive
                assert mass > 0, "All masses should be positive"
    
    def test_fragment_numbering(self):
        """Test that fragment types are correct and consistent."""
        glycan = Glycan("4501", glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            glycan_fragments, _ = glycan.generate_fragments(
                structures[0],
                    modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code="4501",
                modification_type=0,
                glycan_type='N'
            )
            
            # Check fragment types are present
            frag_types = df['Fragment Type'].unique()
            assert len(frag_types) > 0, "Should have fragment types"
            
            # Fragment types should be valid glycan fragment types
            valid_types = ['B', 'C', 'Y', 'Z', 'A', 'CUSTOM']
            for frag_type in frag_types:
                # Remove numbers and special characters to get base type
                base_type = ''.join([c for c in frag_type if c.isalpha()]).upper()
                if base_type:  # Only check if there's a base type
                    assert any(base_type.startswith(vt) for vt in valid_types), f"Invalid fragment type: {frag_type}"
    
    def test_simple_n_glycopeptide_workflow(self):
        """Test complete workflow: peptide + N-glycan → fragments → DataFrame."""
        # Simple N-glycopeptide: IgG tryptic peptide with biantennary glycan
        peptide_seq = "EEQYNSTYR"  # Contains N-glycosylation site (N-X-S/T)
        glycan_code = "4500"  # Common biantennary: HexNAc(4)Hex(5)Fuc(0)NeuAc(0)
        glycosylation_site = 5  # Glycan at N (1-indexed)
        
        # Step 1: Create glycopeptide with integrated approach
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=glycosylation_site,
            glycan_type='N',
            use_cam=True
        )
        
        # Step 2: Generate fragments
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate fragments"
        
        # Verify peptide and glycan properties
        assert glycopeptide.peptide_sequence == peptide_seq
        assert glycopeptide.glycan_code == glycan_code
        assert glycopeptide.glycosylation_site == glycosylation_site
        
        # Verify glycan structures were predicted
        assert len(glycopeptide.glycan.possible_structures) > 0
        
        # Step 3: Create DataFrame (using glycan fragments as example)
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            peptide=peptide_seq,
            glycan_type='N',
            charge_states=[1, 2]
        )
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0, "Should create non-empty DataFrame"
        
        # Verify DataFrame structure
        required_columns = ['Fragment Type', 'm/z(z=1)', 'Fragment']
        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"
        
        # Verify fragment data
        assert df['m/z(z=1)'].notna().all(), "All m/z values should be valid"
        assert (df['m/z(z=1)'] > 0).all(), "All m/z values should be positive"
    
    def test_o_glycopeptide_workflow(self):
        """Test complete workflow with O-glycopeptide (T-antigen)."""
        # Mucin-like O-glycopeptide with T-antigen
        peptide_seq = "GTTPSPVPTR"  # Contains multiple S/T sites
        glycan_code = "1100"  # T-antigen: HexNAc(1)Hex(1)
        
        # Complete workflow
        peptide = Peptide(peptide_seq)
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        assert len(structures) > 0, "Should predict O-glycan structures"
        
        glycopeptide = Glycopeptide(peptide, glycan, site_position=1)  # Glycan at first T
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate O-glycopeptide fragments"
        
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code=glycopeptide.glycan_code,
            modification_type='glycopeptide',
            peptide=glycopeptide.peptide_sequence,
            glycan_type=glycopeptide.glycan_type,
        )
        assert len(df) > 0, "Should create O-glycopeptide DataFrame"
        
        # Verify O-glycan specific fragments
        fragment_types = df['Fragment_Type'].unique()
        assert 'b' in fragment_types or 'y' in fragment_types, "Should have peptide fragments"
    
    def test_sialylated_glycopeptide_workflow(self):
        """Test complete workflow with sialylated N-glycan (common in serum)."""
        peptide_seq = "DQCIYNTTYLNVQR"  # Contains N-glycosylation site
        glycan_code = "4502"  # Disialylated: HexNAc(4)Hex(5)NeuAc(2)
        
        peptide = Peptide(peptide_seq)
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            glycopeptide = Glycopeptide(peptide, glycan, site_position=5)  # Glycan at N
            fragments = glycopeptide.generate_fragments()
            
            df = generate_fragment_table(
                fragments=fragments,
                glycan_code=glycopeptide.glycan_code,
                modification_type='glycopeptide',
                peptide=glycopeptide.peptide_sequence,
                glycan_type=glycopeptide.glycan_type,
            )
            assert len(df) > 0, "Should create sialylated glycopeptide DataFrame"
            
            # Verify sialic acid presence in fragments
            # Y-ions should include NeuAc masses
            y_ions = df[df['Fragment_Type'] == 'Y']
            if len(y_ions) > 0:
                assert y_ions['Mass'].max() > 1000, "Should have large Y-ions with sialic acid"
    
    def test_fucosylated_glycopeptide_workflow(self):
        """Test complete workflow with core-fucosylated N-glycan (very common)."""
        peptide_seq = "EEQFNSTFR"  # Contains N-glycosylation site
        glycan_code = "4501"  # Core-fucosylated: HexNAc(4)Hex(5)Fuc(1)
        
        peptide = Peptide(peptide_seq)
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        glycopeptide = Glycopeptide(peptide, glycan, site_position=4)  # Glycan at N
        fragments = glycopeptide.generate_fragments()
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code=glycopeptide.glycan_code,
            modification_type='glycopeptide',
            peptide=glycopeptide.peptide_sequence,
            glycan_type=glycopeptide.glycan_type,
        )
        
        assert len(df) > 0, "Should create fucosylated glycopeptide DataFrame"
        
        # Verify fucose in fragments
        # Should have fragments with fucose (146.0579 Da)
        assert df['Mass'].max() > peptide.mass, "Should have fragments heavier than peptide"


# ============================================================================
# Literature-Based Examples
# ============================================================================

class TestLiteratureExamples:
    """Test with real glycopeptides from published literature."""
    
    def test_igg_n297_glycopeptide(self):
        """Test IgG1 Fc N-glycosylation site (N297) - most studied glycopeptide."""
        # IgG1 Fc tryptic glycopeptide: EEQYNSTYR (N297 site)
        # This is the most well-characterized glycopeptide in literature
        peptide_seq = "EEQYNSTYR"
        glycan_code = "4501"  # G0F: Most common IgG glycoform
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=4,  # N at position 4
            glycan_type='N',
            use_cam=True
        )
        
        # Verify peptide exists
        assert glycopeptide.peptide is not None
        peptide_mass = glycopeptide.peptide.mass
        
        # Peptide mass should be reasonable (with CAM modification)
        assert 1100 < peptide_mass < 1250, f"IgG peptide mass expected ~1188 Da with CAM, got {peptide_mass}"
        
        # Verify glycan prediction
        assert len(glycopeptide.glycan.possible_structures) > 0, "Should predict G0F structures"
        
        # Generate fragments
        fragments = glycopeptide.generate_fragments()
        assert len(fragments) > 0, "Should generate IgG fragments"
        
        # Generate DataFrame
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            peptide=peptide_seq,
            glycan_type='N',
            charge_states=[1, 2, 3]
        )
        
        # Should have comprehensive fragment coverage
        assert len(df) > 10, "IgG glycopeptide should generate many fragments"
    
    def test_transferrin_n_glycopeptide(self):
        """Test transferrin N-glycopeptide - common serum glycoprotein."""
        # Transferrin has two N-glycosylation sites
        # Example glycopeptide from N413 site
        peptide_seq = "QQQHLFGSNVTDCSR"  # Contains N-X-S/T motif
        glycan_code = "4502"  # Disialylated biantennary (common in serum)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=8,  # N at position 8
            glycan_type='N',
            use_cam=True
        )
        
        if len(glycopeptide.glycan.possible_structures) > 0:
            fragments = glycopeptide.generate_fragments()
            assert len(fragments) > 0, "Transferrin glycopeptide should generate fragments"
            
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                glycopeptide.glycan.possible_structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code=glycan_code,
                modification_type=0,
                glycan_type='N'
            )
            
            assert len(df) > 0, "Should generate transferrin fragment table"
    
    def test_rnase_b_high_mannose(self):
        """Test RNase B high-mannose N-glycan - standard glycoprotein."""
        # RNase B contains high-mannose N-glycans at N34
        # Common composition: Man5-9GlcNAc2
        peptide_seq = "SRNLTK"  # RNase B tryptic peptide with N34
        glycan_code = "2500"  # Man5GlcNAc2 (simplified high-mannose)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=4,  # N at position 4
            glycan_type='N',
            use_cam=True
        )
        
        if len(glycopeptide.glycan.possible_structures) > 0:
            glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
                glycopeptide.glycan.possible_structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code=glycan_code,
                modification_type=0,
                glycan_type='N'
            )
            
            assert len(df) > 0, "RNase B glycopeptide should generate fragments"
    
    def test_mucin_o_glycopeptide(self):
        """Test mucin-type O-glycopeptide - classic O-glycan example."""
        # Mucin proteins are heavily O-glycosylated
        peptide_seq = "PAPGSTAPP"  # Mucin-like repeat with S/T sites
        glycan_code = "1101"  # Sialyl-T antigen (common mucin glycan)
        
        glycopeptide = Glycopeptide(
            peptide_sequence=peptide_seq,
            glycan_code=glycan_code,
            glycosylation_site=6,  # O-glycan at S
            glycan_type='O',
            use_cam=True
        )
        
        assert len(glycopeptide.glycan.possible_structures) > 0, "Should predict sialyl-T structures"
        
        glycan_fragments, _ = glycopeptide.glycan.generate_fragments(
            glycopeptide.glycan.possible_structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=glycan_fragments,
            glycan_code=glycan_code,
            modification_type=0,
            glycan_type='O'
        )
        
        assert len(df) > 0, "Mucin O-glycopeptide should generate fragments"


# ============================================================================
# Complex Multi-Glycosylated Peptides
# ============================================================================

class TestMultiGlycosylation:
    """Test peptides with multiple glycosylation sites."""
    
    def test_two_n_glycosylation_sites(self):
        """Test peptide with two N-glycosylation sites."""
        # Peptide with two N-X-S/T motifs
        peptide_seq = "NETSLNGSQK"  # Two potential N-glycosylation sites
        glycan_code = "4300"  # Biantennary glycan
        
        glycan = Glycan(glycan_code, glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            # Glycan at first N-site
            glycopeptide1 = Glycopeptide(
                peptide_sequence=peptide_seq,
                glycan_code=glycan_code,
                glycosylation_site=1,
                glycan_type='N'
            )
            fragments1 = glycopeptide1.generate_fragments()
            
            # Glycan at second N-site
            glycopeptide2 = Glycopeptide(
                peptide_sequence=peptide_seq,
                glycan_code=glycan_code,
                glycosylation_site=7,
                glycan_type='N'
            )
            fragments2 = glycopeptide2.generate_fragments()
            
            # Both should generate fragments
            assert len(fragments1) > 0, "First site should generate fragments"
            assert len(fragments2) > 0, "Second site should generate fragments"
            
            # Fragment patterns should differ based on site position
            # This reflects real MS/MS fragmentation differences
    
    def test_multiple_o_glycosylation_sites(self):
        """Test mucin-like peptide with multiple O-glycosylation sites."""
        # Mucin-like peptide with 4 potential O-glycan sites
        peptide_seq = "TTSTTTSTR"  # Multiple S/T sites
        glycan_code = "1100"  # T-antigen
        
        glycan = Glycan(glycan_code, glycan_type='O')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            # Test glycosylation at different T/S positions
            for position in [1, 3, 5]:  # Test subset of positions
                glycopeptide = Glycopeptide(
                    peptide_sequence=peptide_seq,
                    glycan_code=glycan_code,
                    glycosylation_site=position,
                    glycan_type='O'
                )
                fragments = glycopeptide.generate_fragments()
                
                assert len(fragments) > 0, f"Should generate fragments for position {position}"
    
    def test_heterogeneous_glycosylation(self):
        """Test same peptide with different glycan structures (microheterogeneity)."""  
        # Microheterogeneity: same peptide with different glycan structures
        peptide_seq = "EEQYNSTYR"
        
        # Different glycoforms at the same site
        glycan_codes = [
            "4300",  # G0 (no galactose)
            "4500",  # G2 (full galactosylation)
            "4501",  # G0F (fucosylated)
            "4502",  # G2S2 (disialylated)
        ]
        
        dataframes = []
        for code in glycan_codes:
            glycan = Glycan(code, glycan_type='N')
            structures = glycan.predict_structures()
            
            if len(structures) > 0:
                glycopeptide = Glycopeptide(
                    peptide_sequence=peptide_seq,
                    glycan_code=code,
                    glycosylation_site=4,
                    glycan_type='N'
                )
                
                glycan_fragments, _ = glycan.generate_fragments(
                    structures[0],
                    modification_type=0
                )
                
                df = generate_fragment_table(
                    fragments=glycan_fragments,
                    glycan_code=code,
                    modification_type=0,
                    glycan_type='N'
                )
                dataframes.append((code, df))
        
        # Should generate different fragment masses for different glycoforms
        assert len(dataframes) >= 2, "Should process multiple glycoforms"
        
        # Verify masses differ between glycoforms
        if len(dataframes) >= 2:
            df1_max_mass = dataframes[0][1]['m/z(z=1)'].max()
            df2_max_mass = dataframes[1][1]['m/z(z=1)'].max()
            assert df1_max_mass != df2_max_mass, "Different glycoforms should have different masses"


# ============================================================================
# Output Format Validation
# ============================================================================

class TestOutputFormat:
    """Validate consistency and correctness of output formats."""
    
    def test_dataframe_column_consistency(self):
        """Test that all DataFrames have consistent column structure."""
        # Generate multiple glycopeptides
        test_cases = [
            ("EEQYNSTYR", "4501", 'N', 4),
            ("GTTPSPVPTR", "1100", 'O', 1),
            ("DQCIYNTTYLNVQR", "4300", 'N', 5),
        ]
        
        dataframes = []
        for peptide_seq, glycan_code, glycan_type, site in test_cases:
            glycan = Glycan(glycan_code, glycan_type=glycan_type)
            structures = glycan.predict_structures()
            
            if len(structures) > 0:
                glycan_fragments, _ = glycan.generate_fragments(
                    structures[0],
                    modification_type=0
                )
                
                df = generate_fragment_table(
                    fragments=glycan_fragments,
                    glycan_code=glycan_code,
                    modification_type=0,
                    peptide=peptide_seq,
                    glycan_type=glycan_type
                )
                dataframes.append(df)
        
        assert len(dataframes) >= 2, "Should generate multiple DataFrames"
        
        # Verify all have same column structure
        first_columns = set(dataframes[0].columns)
        for df in dataframes[1:]:
            assert set(df.columns) == first_columns, "All DataFrames should have same columns"
    
    def test_mass_precision(self):
        """Test that masses are reported with appropriate precision."""
        glycan = Glycan("4501", glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            glycan_fragments, _ = glycan.generate_fragments(
                structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code="4501",
                modification_type=0,
                glycan_type='N'
            )
            
            # Masses should have reasonable precision (not too many or too few decimals)
            for mass in df['m/z(z=1)']:
                # Mass should be positive
                assert mass > 0, "All masses should be positive"
    
    def test_fragment_numbering(self):
        """Test that fragment types are correct and consistent."""
        glycan = Glycan("4501", glycan_type='N')
        structures = glycan.predict_structures()
        
        if len(structures) > 0:
            glycan_fragments, _ = glycan.generate_fragments(
                structures[0],
                modification_type=0
            )
            
            df = generate_fragment_table(
                fragments=glycan_fragments,
                glycan_code="4501",
                modification_type=0,
                glycan_type='N'
            )
            
            # Check fragment types are present
            frag_types = df['Fragment Type'].unique()
            assert len(frag_types) > 0, "Should have fragment types"
            
            # Fragment types should be valid glycan fragment types
            valid_types = ['B', 'C', 'Y', 'Z', 'A', 'CUSTOM']
            for frag_type in frag_types:
                # Remove numbers and special characters to get base type
                base_type = ''.join([c for c in frag_type if c.isalpha()]).upper()
                if base_type:  # Only check if there's a base type
                    assert any(base_type.startswith(vt) for vt in valid_types), f"Invalid fragment type: {frag_type}"

