"""
Test suite for glycofrag.io.tables module
Tests fragment table generation, deduplication, and Excel export functionality.
"""
import pytest
import pandas as pd
import numpy as np
import os
from glycofrag import Glycan, Peptide, Glycopeptide
from glycofrag.io.tables import (
    format_fragment_string,
    extract_fragment_composition,
    generate_fragment_table,
    generate_all_fragments_table,
    deduplicate_fragments_by_mass,
    export_fragment_table_to_excel
)


class TestFormatFragmentString:
    """Test fragment string formatting."""
    
    def test_basic_format(self):
        """Test basic fragment string formatting."""
        composition = {'HexNAc': 2, 'Hex': 3}
        result = format_fragment_string(composition, 'B')
        assert 'HexNAc' in result
        assert 'Hex' in result
        assert result.endswith('-B')
        
    def test_single_monosaccharide(self):
        """Test formatting with single monosaccharide type."""
        composition = {'Hex': 3}
        # Y-ions get PREFIX format: Y-<composition>-<suffix>
        result = format_fragment_string(composition, 'Y')
        assert 'Hex3' in result
        assert result.startswith('Y-')
        
    def test_empty_composition(self):
        """Test formatting with empty composition."""
        composition = {}
        result = format_fragment_string(composition, 'B')
        assert result == '-B'
        
    def test_multiple_monosaccharides(self):
        """Test formatting with many monosaccharide types."""
        composition = {'HexNAc': 2, 'Hex': 3, 'Fuc': 1, 'NeuAc': 2}
        result = format_fragment_string(composition, 'BY')
        assert 'HexNAc' in result
        assert 'Hex' in result
        assert 'Fuc' in result
        assert 'NeuAc' in result
        assert '-B' in result  # BY format uses Y-<comp>-B


class TestExtractFragmentComposition:
    """Test fragment composition extraction."""
    
    def test_basic_extraction(self):
        """Test basic composition extraction from fragment string."""
        fragment_str = "HexNAc2-Hex3-B"
        result = extract_fragment_composition(fragment_str)
        assert result == {'HexNAc': 2, 'Hex': 3}
        
    def test_single_monosaccharide_extraction(self):
        """Test extraction with single monosaccharide."""
        fragment_str = "Hex3-Y"
        result = extract_fragment_composition(fragment_str)
        assert result == {'Hex': 3}
        
    def test_complex_extraction(self):
        """Test extraction with multiple monosaccharides."""
        fragment_str = "HexNAc2-Hex3-Fuc1-NeuAc2-BY"
        result = extract_fragment_composition(fragment_str)
        assert result == {'HexNAc': 2, 'Hex': 3, 'Fuc': 1, 'NeuAc': 2}
        
    def test_empty_composition_extraction(self):
        """Test extraction from empty composition fragment."""
        fragment_str = "-B"
        result = extract_fragment_composition(fragment_str)
        assert result == {}


class TestGenerateFragmentTable:
    """Test fragment table generation."""
    
    def test_simple_glycan_table(self):
        """Test table generation for simple glycan."""
        glycan = Glycan('4503')  # HexNAc(4)Hex(5)Fuc(3)
        structures = glycan.predict_structures()
        assert len(structures) > 0, "No structures predicted"
        
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code='4503',
            modification_type=0,
            charge_states=[1]
        )
        
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert 'Fragment' in df.columns
        assert 'Fragment Type' in df.columns
        assert 'Mass(Da)' in df.columns
        assert 'm/z(z=1)' in df.columns
        
    def test_multiple_charge_states(self):
        """Test table generation with multiple charge states."""
        glycan = Glycan('4501')  # HexNAc(4)Hex(5)Fuc(1)
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code='4501',
            modification_type=0,
            charge_states=[1, 2, 3]
        )
        
        # Check for different charge state columns
        assert 'm/z(z=1)' in df.columns
        assert 'm/z(z=2)' in df.columns
        assert 'm/z(z=3)' in df.columns
        # All m/z columns should have values
        assert df['m/z(z=1)'].notna().any()
        assert df['m/z(z=2)'].notna().any()
        assert df['m/z(z=3)'].notna().any()
        
    def test_with_peptide_fragments(self):
        """Test table generation including peptide fragments."""
        peptide = Peptide('YPVLN')
        glycan = Glycan('3300')  # Small glycan HexNAc(3)Hex(3)
        structures = glycan.predict_structures()
        
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code='3300',
            modification_type=0,
            peptide=peptide.sequence,
            charge_states=[1]
        )
        
        # Should have glycan fragments
        assert len(df) > 0
        # Glycan-only table doesn't include Glycopeptide column
        # (that's added by Glycopeptide workflow)
        assert 'Fragment' in df.columns
        
    def test_dataframe_structure(self):
        """Test that generated DataFrame has correct structure."""
        glycan = Glycan('3300')  # HexNAc(3)Hex(3) - valid N-glycan
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=0
        )
        
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code='3300',
            modification_type=0,
            charge_states=[1]
        )
        
        required_columns = ['Fragment', 'Fragment Type', 'Mass(Da)', 'm/z(z=1)']
        for col in required_columns:
            assert col in df.columns, f"Missing column: {col}"
            
        # Check data types
        assert pd.api.types.is_numeric_dtype(df['m/z(z=1)'])


class TestGenerateAllFragmentsTable:
    """Test generation of combined fragment tables."""
    
    def test_multiple_structures(self):
        """Test table generation for multiple glycan structures."""
        glycan1 = Glycan('3300')  # HexNAc(3)Hex(3)
        structures1 = glycan1.predict_structures()
        fragments1, _ = glycan1.generate_fragments(structures1[0], modification_type=0)
        
        glycan2 = Glycan('4400')  # HexNAc(4)Hex(4)
        structures2 = glycan2.predict_structures()
        fragments2, _ = glycan2.generate_fragments(structures2[0], modification_type=0)
        
        results = {
            '3300': {'fragments': fragments1},
            '4400': {'fragments': fragments2}
        }
        
        df = generate_all_fragments_table(
            results=results,
            modification_type=0,
            charge_states=[1]
        )
        
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        # Should contain fragments from both structures
        assert len(df) > 0
        
    def test_empty_results(self):
        """Test handling of empty results."""
        results = {}
        df = generate_all_fragments_table(
            results=results,
            modification_type=0,
            charge_states=[1]
        )
        
        # Should return empty DataFrame with correct columns
        assert isinstance(df, pd.DataFrame)
        assert 'Fragment' in df.columns or 'm/z(z=1)' in df.columns


class TestDeduplicateFragmentsByMass:
    """Test fragment deduplication by mass."""
    
    def test_basic_deduplication(self):
        """Test basic deduplication with duplicate masses."""
        # Create test DataFrame with duplicates
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 100.0, 200.0, 100.0],
            'FragmentType': ['B', 'Y', 'B', 'BY'],
            'Ions': ['1H+', '1H+', '1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)', 'N(H)', 'N(H)']
        })
        
        result = deduplicate_fragments_by_mass(df)
        
        # Should have only unique masses
        assert len(result) == 2  # 100.0 and 200.0
        assert len(result['Fragment_mz'].unique()) == 2
        
    def test_priority_b_over_y(self):
        """Test that B fragments have priority over Y fragments."""
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 100.0],
            'FragmentType': ['Y', 'B'],
            'Ions': ['1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)']
        })
        
        result = deduplicate_fragments_by_mass(df)
        
        assert len(result) == 1
        assert result.iloc[0]['FragmentType'] == 'B'
        
    def test_priority_charge_1_over_2(self):
        """Test that charge 1+ has priority over 2+."""
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 100.0],
            'FragmentType': ['B', 'B'],
            'Ions': ['2H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)']
        })
        
        result = deduplicate_fragments_by_mass(df)
        
        assert len(result) == 1
        assert result.iloc[0]['Ions'] == '1H+'
        
    def test_no_duplicates(self):
        """Test deduplication with no duplicates."""
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 200.0, 300.0],
            'FragmentType': ['B', 'Y', 'BY'],
            'Ions': ['1H+', '1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)', 'N(H)']
        })
        
        result = deduplicate_fragments_by_mass(df)
        
        assert len(result) == 3
        assert len(result) == len(df)


class TestExportFragmentTableToExcel:
    """Test Excel export functionality."""
    
    def test_basic_export(self, tmp_path):
        """Test basic Excel export."""
        # Create test DataFrame
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 200.0],
            'FragmentType': ['B', 'Y'],
            'Ions': ['1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)']
        })
        
        filename = tmp_path / "test_fragments.xlsx"
        export_fragment_table_to_excel(df, None, str(filename))
        
        assert filename.exists()
        
        # Read back and verify
        df_read = pd.read_excel(filename, sheet_name='Fragments')
        assert len(df_read) == len(df)
        
    def test_export_with_all_fragments(self, tmp_path):
        """Test export with both individual and all fragments."""
        df = pd.DataFrame({
            'Fragment_mz': [100.0, 200.0],
            'FragmentType': ['B', 'Y'],
            'Ions': ['1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)']
        })
        
        df_all = pd.DataFrame({
            'Fragment_mz': [100.0, 200.0, 300.0],
            'FragmentType': ['B', 'Y', 'BY'],
            'Ions': ['1H+', '1H+', '1H+'],
            'Glycan': ['N(H)', 'N(H)', 'N(H)']
        })
        
        filename = tmp_path / "test_fragments_all.xlsx"
        export_fragment_table_to_excel(df, df_all, str(filename))
        
        assert filename.exists()
        
        # Read back and verify both sheets
        df_read = pd.read_excel(filename, sheet_name='Fragments')
        df_all_read = pd.read_excel(filename, sheet_name='All_Fragments')
        
        assert len(df_read) == len(df)
        assert len(df_all_read) == len(df_all)
        
    def test_export_preserves_data(self, tmp_path):
        """Test that export preserves data integrity."""
        df = pd.DataFrame({
            'Fragment_mz': [123.456, 789.012],
            'FragmentType': ['B', 'Y'],
            'Ions': ['1H+', '2H+'],
            'Glycan': ['N(H)', 'N(H)']
        })
        
        filename = tmp_path / "test_integrity.xlsx"
        export_fragment_table_to_excel(df, None, str(filename))
        
        df_read = pd.read_excel(filename, sheet_name='Fragments')
        
        # Check numeric precision
        assert abs(df_read['Fragment_mz'].iloc[0] - 123.456) < 0.001
        assert abs(df_read['Fragment_mz'].iloc[1] - 789.012) < 0.001


class TestIntegrationWorkflow:
    """Test complete workflow integration."""
    
    def test_end_to_end_workflow(self, tmp_path):
        """Test complete workflow from glycan to Excel export."""
        # 1. Create glycan and generate fragments
        glycan = Glycan('4501')  # HexNAc(4)Hex(5)Fuc(1)
        structures = glycan.predict_structures()
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=0
        )
        
        # 2. Generate fragment table
        df = generate_fragment_table(
            fragments=fragments,
            glycan_code='4501',
            modification_type=0,
            charge_states=[1, 2]
        )
        
        # 3. Deduplicate
        df_dedup = deduplicate_fragments_by_mass(df)
        
        # 4. Export to Excel
        filename = tmp_path / "workflow_test.xlsx"
        export_fragment_table_to_excel(df_dedup, df, str(filename))
        
        # Verify
        assert filename.exists()
        df_read = pd.read_excel(filename, sheet_name='Fragments')
        assert len(df_read) > 0
        assert len(df_read) <= len(df)  # Deduplication should reduce or keep same
        
    def test_multiple_glycans_workflow(self, tmp_path):
        """Test workflow with multiple glycan structures."""
        # Generate fragments for multiple glycans
        glycans = ['3300', '4501']  # HexNAc3Hex3, HexNAc4Hex5Fuc1
        results = {}
        
        for glycan_code in glycans:
            glycan = Glycan(glycan_code)
            structures = glycan.predict_structures()
            if len(structures) > 0:
                fragments, _ = glycan.generate_fragments(
                    structures[0],
                    modification_type=0
                )
                results[glycan_code] = {'fragments': fragments}
            
        # Generate combined table
        df_all = generate_all_fragments_table(
            results=results,
            modification_type=0,
            charge_states=[1]
        )
        
        # Deduplicate and export
        df_dedup = deduplicate_fragments_by_mass(df_all)
        filename = tmp_path / "multi_glycan_test.xlsx"
        export_fragment_table_to_excel(df_dedup, df_all, str(filename))
        
        # Verify
        assert filename.exists()
        df_read = pd.read_excel(filename, sheet_name='All_Fragments')
        assert len(df_read) > 0

