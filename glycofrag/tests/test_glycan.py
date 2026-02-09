"""
Tests for the Glycan class - structure prediction and analysis.
Tests use public API only (behavior testing, not implementation testing).
"""

import pytest
import networkx as nx

from glycofrag import Glycan


class TestGlycanInitialization:
    """Test Glycan object initialization and parsing."""
    
    def test_init_simple_nglycan(self):
        """Test initialization with simple N-glycan code produces valid structures."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Verify structures were created
        assert len(structures) > 0
        
        # Verify structure composition matches input
        counts = glycan.count_residues(structures[0])
        assert counts['HexNAc'] == 4
        assert counts['Hex'] == 5
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 1
    
    def test_init_oglycan(self):
        """Test initialization with O-glycan code."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        
        # Verify structures were created
        assert len(structures) > 0
        
        # Verify structure composition
        counts = glycan.count_residues(structures[0])
        assert counts['HexNAc'] == 1
        assert counts['Hex'] == 1
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 0
    
    def test_init_with_named_code(self):
        """Test initialization with named format code."""
        glycan = Glycan('HexNAc(2)Hex(3)', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Verify structure composition
        counts = glycan.count_residues(structures[0])
        assert counts['HexNAc'] == 2
        assert counts['Hex'] == 3
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 0
    
    def test_max_structures_parameter(self):
        """Test max_structures parameter is respected."""
        glycan = Glycan('2300', glycan_type='N', max_structures=50)
        assert glycan.max_structures == 50


class TestNGlycanCoreBuilding:
    """Test N-glycan core structure building through structure prediction."""
    
    def test_build_n_glycan_core(self):
        """Test N-glycan core structure has correct composition."""
        glycan = Glycan('2300', glycan_type='N')
        structures = glycan.predict_structures()
        
        # Verify structure created
        assert len(structures) > 0
        core = structures[0]
        
        # Check node count (2 HexNAc + 3 Hex)
        assert len(core.nodes()) == 5
        
        # Check node types using public API
        counts = glycan.count_residues(core)
        assert counts['HexNAc'] == 2
        assert counts['Hex'] == 3
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 0
    
    def test_n_glycan_core_structure(self):
        """Test N-glycan core has correct graph structure."""
        glycan = Glycan('2300', glycan_type='N')
        structures = glycan.predict_structures()
        core = structures[0]
        
        # Node 1 should be reducing end
        assert core.nodes[1]['type'] == 'HexNAc'
        assert core.nodes[1]['position'] == 'core_reducing'
        
        # Node 3 should be central mannose with 2 branches
        assert core.nodes[3]['type'] in ['Hex', 'Man']  # Could be labeled as Man
        assert core.out_degree(3) == 2  # Branches to nodes 4 and 5


class TestOGlycanCoreBuilding:
    """Test O-glycan core structure building through prediction."""
    
    def test_build_o_glycan_cores(self):
        """Test O-glycan generates structures with correct composition."""
        glycan = Glycan('3210', glycan_type='O')
        structures = glycan.predict_structures()
        
        # Should generate structures
        assert len(structures) > 0
        
        # Verify composition in generated structures
        for structure in structures:
            counts = glycan.count_residues(structure)
            # Total should match input
            total = counts['HexNAc'] + counts['Hex'] + counts['Fuc'] + counts['NeuAc']
            assert total == 3 + 2 + 1 + 0  # 3 HexNAc, 2 Hex, 1 Fuc, 0 NeuAc
    
    def test_o_glycan_core_type_1(self):
        """Test O-glycan Core 1 (HexNAc-Hex) T antigen."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        
        # Should generate at least one structure
        assert len(structures) >= 1
        
        # Check first structure is simple T antigen
        structure = structures[0]
        assert len(structure.nodes()) == 2
        counts = glycan.count_residues(structure)
        assert counts['HexNAc'] == 1
        assert counts['Hex'] == 1


class TestGlycanClassification:
    """Test glycan structure classification."""
    
    def test_classify_high_mannose(self):
        """Test high mannose N-glycan classification."""
        glycan = Glycan('2300', glycan_type='N')  # HexNAc(2)Hex(3)
        structures = glycan.predict_structures()
        classification = glycan.classify_structure(structures[0])
        assert classification == "High Mannose"
    
    def test_classify_complex(self):
        """Test complex N-glycan classification."""
        glycan = Glycan('4501', glycan_type='N')  # HexNAc(4)Hex(5)NeuAc(1)
        structures = glycan.predict_structures()
        classification = glycan.classify_structure(structures[0])
        assert classification == "Complex"
    
    def test_classify_hybrid(self):
        """Test hybrid N-glycan classification."""
        glycan = Glycan('3400', glycan_type='N')  # HexNAc(3)Hex(4)
        structures = glycan.predict_structures()
        classification = glycan.classify_structure(structures[0])
        assert classification == "Hybrid"
    
    def test_classify_o_glycan(self):
        """Test O-glycan classification."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        classification = glycan.classify_structure(structures[0])
        # O-glycans should be classified as such
        assert "Core" in classification or "O-GalNAc" in classification or classification == "O-glycan"


class TestStructurePrediction:
    """Test glycan structure prediction."""
    
    def test_predict_simple_nglycan(self):
        """Test prediction for simple N-glycan generates structures."""
        glycan = Glycan('2300', glycan_type='N', max_structures=10)
        structures = glycan.predict_structures()
        
        assert len(structures) > 0
        assert all(isinstance(s, nx.DiGraph) for s in structures)
    
    def test_predict_special_case_hexnac1(self):
        """Test special case: single HexNAc."""
        glycan = Glycan('1000', glycan_type='N')
        structures = glycan.predict_structures()
        
        assert len(structures) == 1
        assert len(structures[0].nodes()) == 1
    
    def test_predict_special_case_t_antigen(self):
        """Test special case: T antigen (HexNAc-Hex)."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        
        assert len(structures) == 1
        assert len(structures[0].nodes()) == 2
    
    def test_max_structures_limit(self):
        """Test that max_structures limit is respected."""
        glycan = Glycan('4501', glycan_type='N', max_structures=5)
        structures = glycan.predict_structures()
        
        assert len(structures) <= 5


class TestStructureFingerprinting:
    """Test structure uniqueness and deduplication behavior."""
    
    def test_fingerprint_uniqueness(self):
        """Test that different compositions produce different structures."""
        glycan1 = Glycan('2300', glycan_type='N')
        glycan2 = Glycan('2400', glycan_type='N')
        
        structures1 = glycan1.predict_structures()
        structures2 = glycan2.predict_structures()
        
        # Different compositions should produce different structure counts
        counts1 = glycan1.count_residues(structures1[0])
        counts2 = glycan2.count_residues(structures2[0])
        
        assert counts1 != counts2
    
    def test_fingerprint_consistency(self):
        """Test that prediction is consistent (same input = same output)."""
        glycan = Glycan('2300', glycan_type='N')
        
        structures1 = glycan.predict_structures()
        structures2 = glycan.predict_structures()
        
        # Same input should produce same number of structures
        assert len(structures1) == len(structures2)
        
        # Structures should have same composition
        counts1 = [glycan.count_residues(s) for s in structures1]
        counts2 = [glycan.count_residues(s) for s in structures2]
        assert counts1 == counts2


class TestResidueCount:
    """Test monosaccharide counting functions using public API."""
    
    def test_count_residues_simple(self):
        """Test counting residues in simple structure using public API."""
        glycan = Glycan('2300', glycan_type='N')
        structures = glycan.predict_structures()
        counts = glycan.count_residues(structures[0])
        
        assert counts['HexNAc'] == 2
        assert counts['Hex'] == 3
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 0
    
    def test_count_residues_complex(self):
        """Test counting residues in complex structure."""
        glycan = Glycan('4501', glycan_type='N')
        structures = glycan.predict_structures()
        counts = glycan.count_residues(structures[0])
        
        assert counts['HexNAc'] == 4
        assert counts['Hex'] == 5
        assert counts['Fuc'] == 0
        assert counts['NeuAc'] == 1


class TestOGlycanStructurePrediction:
    """Test O-glycan specific structure prediction."""
    
    def test_o_glycan_simple_prediction(self):
        """Test O-glycan structure prediction for simple case."""
        glycan = Glycan('1100', glycan_type='O')
        structures = glycan.predict_structures()
        
        assert len(structures) >= 1
        # Should generate simple HexNAc-Hex (T antigen)
        assert any(len(s.nodes()) == 2 for s in structures)
    
    def test_o_glycan_complex_prediction(self):
        """Test O-glycan structure prediction for complex case."""
        glycan = Glycan('2110', glycan_type='O', max_structures=20)
        structures = glycan.predict_structures()
        
        assert len(structures) > 0
        # All structures should use all sugars - verify with public API
        for structure in structures:
            counts = glycan.count_residues(structure)
            assert counts['HexNAc'] == 2
            assert counts['Hex'] == 1
            assert counts['Fuc'] == 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
