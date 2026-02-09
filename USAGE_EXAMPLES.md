# Glycofrag Usage Examples

Comprehensive examples demonstrating all features of the Glycofrag package.

## Table of Contents

1. [Basic Setup](#basic-setup)
2. [Glycan Analysis](#glycan-analysis)
3. [Peptide Fragmentation](#peptide-fragmentation)
4. [Glycopeptide Analysis](#glycopeptide-analysis)
5. [Real-World Workflows](#real-world-workflows)
6. [Advanced Techniques](#advanced-techniques)

---

## Basic Setup

### Installation and Import

```python
# Install from source
# pip install -e /path/to/glycofrag

# Import main classes
from glycofrag import Glycan, Peptide, Glycopeptide, GlycanMassCalculator
from glycofrag import MONOSACCHARIDE_MASSES, AMINO_ACID_MASSES, MODIFICATION_MASSES

# Verify installation
import glycofrag
print(f"Glycofrag version: {glycofrag.__version__}")
```

---

## Glycan Analysis

### Example 1: Basic Glycan Mass Calculation

```python
from glycofrag import GlycanMassCalculator

# Initialize calculator
calc = GlycanMassCalculator()

# Calculate mass from composition code
# Code format: HEXNAC, HEX, NEUAC, FUC (4 digits)
mass_4501 = calc.calculate_glycan_mass("4501")
print(f"Glycan 4501 mass: {mass_4501:.4f} Da")
# Output: Glycan 4501 mass: 1336.4779 Da (HexNAc-4, Hex-5, NeuAc-0, Fuc-1)

# Different compositions
compositions = ["3300", "4401", "5510", "2200"]
for comp in compositions:
    mass = calc.calculate_glycan_mass(comp)
    print(f"  {comp}: {mass:.2f} Da")
```

### Example 2: Predict N-Glycan Structures

```python
from glycofrag import Glycan

# Create glycan object
glycan = Glycan("4501", glycan_type="N")

# Predict all possible structures
structures = glycan.predict_structures()
print(f"Predicted {len(structures)} possible N-glycan structures from composition 4501")

# Classify each structure
classifications = {}
for i, structure in enumerate(structures):
    classification = glycan.classify_structure(structure)
    classifications[classification] = classifications.get(classification, 0) + 1
    if i < 3:
        print(f"  Structure {i+1}: {classification}")

print(f"\nStructure classifications:")
for classification, count in sorted(classifications.items(), key=lambda x: -x[1]):
    print(f"  {classification}: {count} structures")
```

### Example 3: Predict O-Glycan Structures

```python
from glycofrag import Glycan

# O-glycan example (Core 1)
o_glycan = Glycan("2201", glycan_type="O")
structures = o_glycan.predict_structures()

print(f"O-glycan 2201 (Core 1 structure):")
print(f"  Total structures: {len(structures)}")

# Get structure details
for i, structure in enumerate(structures[:2]):
    classification = o_glycan.classify_structure(structure)
    residues = o_glycan._count_residues(structure)
    print(f"\n  Structure {i+1}: {classification}")
    print(f"    Residues: {residues}")
```

### Example 4: Generate Glycan Fragment Ions

```python
from glycofrag import Glycan

# Create glycan and predict structures
glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()

if structures:
    # Generate fragment ions
    fragments, metadata = glycan.generate_fragments(structures[0])
    
    print(f"Generated {len(fragments)} fragment ions")
    print(f"\nFragment types distribution:")
    
    # Group by type
    by_type = {}
    for frag in fragments:
        ftype = frag['type']
        by_type[ftype] = by_type.get(ftype, 0) + 1
    
    for ftype in sorted(by_type.keys()):
        print(f"  {ftype}: {by_type[ftype]} ions")
    
    # Show some examples
    print(f"\nExample fragment ions:")
    for frag in fragments[:5]:
        print(f"  {frag['name']}: {frag['mass']:.4f} m/z (charge {frag['charge']})")
```

### Example 5: Glycan Fragment Modifications

```python
from glycofrag import Glycan

glycan = Glycan("4501", glycan_type="N")
structures = glycan.predict_structures()

# Generate fragments with different reducing end modifications
print("Fragment generation with different reducing end types:\n")

for mod_type in range(7):
    try:
        fragments, _ = glycan.generate_fragments(
            structures[0],
            modification_type=mod_type,
            fragment_types=['BY']
        )
        print(f"  Type {mod_type}: {len(fragments)} BY-series fragments")
    except:
        pass

# Permethylated glycan
print("\nPermethylated glycan:")
fragments_perm, _ = glycan.generate_fragments(
    structures[0],
    permethylated=True,
    fragment_types=['BY', 'CZ']
)
print(f"  Total fragments: {len(fragments_perm)}")
```

---

## Peptide Fragmentation

### Example 1: Basic Peptide Fragmentation

```python
from glycofrag import Peptide

# Create peptide (no modifications)
peptide = Peptide("PEPTIDE")

# Generate all fragments
fragments = peptide.generate_all_fragments()
print(f"Peptide: {peptide.sequence}")
print(f"Total fragments: {len(fragments)}")

# Count by type
fragment_counts = {}
for frag in fragments:
    ftype = frag['type']
    fragment_counts[ftype] = fragment_counts.get(ftype, 0) + 1

print("\nFragment counts by type:")
for ftype in sorted(fragment_counts.keys()):
    print(f"  {ftype}: {fragment_counts[ftype]}")

# Show example b-ions
print("\nExample b-ions:")
b_ions = [f for f in fragments if f['type'] == 'b']
for ion in sorted(b_ions, key=lambda x: x['mass'])[:5]:
    print(f"  {ion['name']}: {ion['mass']:.4f} m/z")
```

### Example 2: Peptide with CAM Modification

```python
from glycofrag import Peptide

# Peptide with cysteines (Carbamidomethylation) - positions are 1-based
peptide_seq = "LCPDCPLLAPLNDSR"
peptide = Peptide(peptide_seq, modifications={'CAM': [1, 4, 6]})  # L, D, P at positions 1, 4, 6

print(f"Peptide: {peptide_seq}")
print(f"Modifications: CAM at positions 1, 4, 6 (1-based indexing)")

fragments = peptide.generate_all_fragments()
print(f"Total fragments: {len(fragments)}")

# Check that CAM is applied
print("\nFragment mass increase from CAM (+57.021 Da each):")
# Compare with unmodified
peptide_unmod = Peptide(peptide_seq, modifications={})
fragments_unmod = peptide_unmod.generate_all_fragments()

# Get same fragment from both
for frag_mod in fragments[:3]:
    frag_name = frag_mod['fragment_name']
    unmod_frag = [f for f in fragments_unmod if f['fragment_name'] == frag_name]
    if unmod_frag:
        mass_diff = frag_mod['fragment_mass'] - unmod_frag[0]['fragment_mass']
        print(f"  {frag_name}: {mass_diff:.2f} Da difference")
```

### Example 3: Multiple Modification Types

```python
from glycofrag import Peptide

# Peptide with multiple modification types (1-based positions)
peptide = Peptide(
    "MEPEPTIDEC",
    modifications={
        'CAM': [10],       # Carbamidomethylation on cysteine at position 10
        'Oxidation': [1]   # Oxidation on methionine at position 1
    }
)

print("Peptide: MEPEPTIDEC")
print("Modifications: CAM at C (position 10), Oxidation at M (position 1) - 1-based indexing")

fragments = peptide.generate_all_fragments()
print(f"Total fragments: {len(fragments)}")

# Find fragments affected by each modification
print("\nFragments with modifications:")
for frag in fragments:
    if frag['fragment_mass'] > 0:  # Valid fragment
        metadata = frag.get('metadata', {})
        if metadata.get('modifications'):
            print(f"  {frag['fragment_name']}: {metadata['modifications']}")
```

### Example 4: Different Charge States

```python
from glycofrag import Peptide

peptide = Peptide("MEPEPTIDE")
fragments = peptide.generate_all_fragments(charge_states=[1, 2, 3])

print("Fragment ions with multiple charge states:")
print(f"Total fragments (all charges): {len(fragments)}")

# Group by charge state
by_charge = {}
for frag in fragments:
    charge = frag['charge']
    by_charge[charge] = by_charge.get(charge, 0) + 1

print("\nFragments by charge state:")
for charge in sorted(by_charge.keys()):
    print(f"  z={charge}: {by_charge[charge]} ions")

# Show m/z variation for single fragment
print("\nExample: mass variation across charge states")
name_to_check = "y3"
matching = [f for f in fragments if f['fragment_name'] == name_to_check]
for frag in sorted(matching, key=lambda x: x['charge']):
    print(f"  {frag['fragment_name']}^{frag['charge']}: {frag['fragment_mass']:.4f} m/z")
```

---

## Glycopeptide Analysis

### Example 1: Basic N-Glycopeptide

```python
from glycofrag import Glycopeptide

# Simple N-glycopeptide (1-based positions)
gp = Glycopeptide(
    peptide_sequence="LCPDCPLLAPLNDSR",
    glycan_code="4501",
    glycosylation_site=12,      # N at position 12 (1-based indexing)
    glycan_type="N"
)

print(f"Peptide: LCPDCPLLAPLNDSR")
print(f"Glycan: 4501 (HexNAc-4, Hex-5)")
print(f"Glycosylation site: N at position 12 (Asparagine is the 12th amino acid)")

# Generate all fragments
fragments = gp.generate_fragments()
print(f"\nTotal fragments: {len(fragments)}")

# Show fragment type distribution
frag_types = {}
for frag in fragments:
    ftype = frag['type']
    frag_types[ftype] = frag_types.get(ftype, 0) + 1

print("\nFragment types:")
for ftype in sorted(frag_types.keys()):
    print(f"  {ftype}: {frag_types[ftype]} ions")

# Show examples of each type
print("\nExample fragments:")
for ftype in ['Y0', 'Y1', 'Y-glycan', 'b', 'Intact']:
    frag = [f for f in fragments if f['type'] == ftype]
    if frag:
        print(f"  {ftype}: {frag[0]['name']} ({frag[0]['mass']:.4f} m/z)")
```

### Example 2: O-Glycopeptide Analysis

```python
from glycofrag import Glycopeptide

# O-glycopeptide (Core 1) - 1-based positions
o_gp = Glycopeptide(
    peptide_sequence="TPEPTIDES",
    glycan_code="2201",
    glycosylation_site=1,      # Threonine at position 1 (1-based indexing)
    glycan_type="O"
)

print("O-Glycopeptide Analysis:")
print(f"Peptide: TPEPTIDES")
print(f"Glycan: 2201 (Core 1)")
print(f"Glycosylation site: T at position 1 (first amino acid)\n")

fragments = o_gp.generate_fragments()

# Diagnostic ions (oxonium)
oxonium = [f for f in fragments if f['type'] == 'A']
print(f"Diagnostic oxonium ions (A-type): {len(oxonium)}")
for ion in oxonium[:3]:
    print(f"  {ion['name']}: {ion['mass']:.4f} m/z")

# Peptide backbone
y0 = [f for f in fragments if f['type'] == 'Y0']
print(f"\nPeptide backbone ions (Y0): {len(y0)}")

# Complete glycopeptide
intact = [f for f in fragments if f['type'] == 'Intact']
if intact:
    print(f"\nIntact glycopeptide: {intact[0]['mass']:.4f} m/z")
```

### Example 3: Multiple Glycan Structures

```python
from glycofrag import Glycopeptide, Glycan

# Analyze same peptide with different glycans
peptide_seq = "LCPDCPLLAPLNDSR"
position = 2

glycans = ["3300", "4401", "4501", "5510"]

results = []
for glycan_code in glycans:
    gp = Glycopeptide(
        peptide_sequence=peptide_seq,
        glycan_code=glycan_code,
        glycosylation_site=position,
        glycan_type="N"
    )
    
    fragments = gp.generate_fragments()
    
    # Get intact mass
    intact = [f for f in fragments if f['type'] == 'Intact']
    intact_mass = intact[0]['mass'] if intact else 0
    
    results.append({
        'glycan': glycan_code,
        'total_fragments': len(fragments),
        'intact_mass': intact_mass
    })

print("Glycopeptide Analysis - Glycan Comparison:")
print(f"{'Glycan':<8} {'Fragments':<12} {'Intact m/z':<12}")
print("-" * 35)
for r in results:
    print(f"{r['glycan']:<8} {r['total_fragments']:<12} {r['intact_mass']:<12.4f}")
```

### Example 4: Glycopeptide with Modifications

```python
from glycofrag import Glycopeptide

# N-glycopeptide with peptide modifications
gp = Glycopeptide(
    peptide_sequence="MEPTETTQSVKC",
    glycan_code="4501",
    glycosylation_site=4,
    glycan_type="N",
    modifications={
        'CAM': [11],       # Carbamidomethylation on C11
        'Oxidation': [0]   # Oxidation on M0
    }
)

print("Glycopeptide with modifications:")
print("Peptide: MEPTETTQSVKC")
print("Modifications: CAM at C11, Oxidation at M0")
print("Glycan: 4501 at N4\n")

fragments = gp.generate_fragments()
print(f"Total fragments: {len(fragments)}")

# Show fragments with modifications in metadata
print("\nFragments with modifications:")
count = 0
for frag in fragments:
    if frag.get('metadata', {}).get('modifications'):
        print(f"  {frag['name']}: {frag['metadata']['modifications']}")
        count += 1
        if count >= 5:
            print("  ...")
            break
```

---

## Real-World Workflows

### Workflow 1: Experimental Data Matching

```python
from glycofrag import Glycopeptide

# Your experimental m/z values (from MS data)
experimental_mz = [300.5, 400.3, 520.8, 654.2, 1200.4, 1500.3]

# Define suspected glycopeptide
gp = Glycopeptide(
    peptide_sequence="LCPDCPLLAPLNDSR",
    glycan_code="4501",
    glycosylation_site=2,
    glycan_type="N"
)

# Generate theoretical fragments
theoretical = gp.generate_fragments()

print("Experimental Data Matching:")
print(f"{'Exp m/z':<10} {'Fragment':<15} {'Theory m/z':<12} {'Delta (ppm)':<10}")
print("-" * 50)

matches = 0
for exp_mz in experimental_mz:
    # Find close matches (within 50 ppm)
    close = [f for f in theoretical if abs((f['mass'] - exp_mz)/exp_mz * 1e6) < 50]
    
    if close:
        for theory_frag in close[:1]:  # Show best match
            ppm = (theory_frag['mass'] - exp_mz) / exp_mz * 1e6
            print(f"{exp_mz:<10.2f} {theory_frag['name']:<15} {theory_frag['mass']:<12.4f} {ppm:<10.2f}")
            matches += 1
    else:
        print(f"{exp_mz:<10.2f} {'No match':<15}")

print(f"\nMatches found: {matches}/{len(experimental_mz)}")
```

### Workflow 2: Database Search Preparation

```python
from glycofrag import Glycopeptide, Glycan

# Prepare a database of theoretical masses for MS/MS search

peptides = [
    ("EPTETTQSVK", 2),
    ("TPEPTIDES", 0),
    ("MEPEPTIDEC", 3),
]

glycan_codes = ["3300", "4401", "4501"]

database = []

for peptide, position in peptides:
    for glycan_code in glycan_codes:
        gp = Glycopeptide(
            peptide_sequence=peptide,
            glycan_code=glycan_code,
            glycosylation_site=position,
            glycan_type="N"
        )
        
        fragments = gp.generate_fragments()
        
        for frag in fragments:
            if frag['type'] in ['Y0', 'Y1', 'Intact']:
                database.append({
                    'peptide': peptide,
                    'glycan': glycan_code,
                    'fragment_type': frag['type'],
                    'mass': frag['mass'],
                    'charge': frag['charge'],
                    'name': frag['name']
                })

print(f"Database entries: {len(database)}")
print(f"Unique peptides: {len(peptides)}")
print(f"Glycan compositions: {len(glycan_codes)}")

# Show sample entries
print("\nSample database entries (first 10):")
print(f"{'Peptide':<15} {'Glycan':<8} {'Type':<10} {'m/z':<12}")
print("-" * 50)
for entry in database[:10]:
    print(f"{entry['peptide']:<15} {entry['glycan']:<8} {entry['fragment_type']:<10} {entry['mass']:<12.4f}")
```

### Workflow 3: Comparative Glycosylation Analysis

```python
from glycofrag import Glycopeptide

# Compare different glycosylation patterns on same peptide

peptide = "LCPDCPLLAPLNDSR"
common_n_glycans = ["2300", "3300", "4401", "4501", "5510"]

print("Comparative Analysis: Same Peptide, Different Glycans\n")
print(f"{'Glycan':<10} {'Composition':<30} {'Y1 m/z':<12} {'Intact m/z':<12}")
print("-" * 70)

compositions_info = {
    "2300": "HexNAc-2, Hex-3",
    "3300": "HexNAc-3, Hex-3",
    "4401": "HexNAc-4, Hex-4, Fuc-1",
    "4501": "HexNAc-4, Hex-5",
    "5510": "HexNAc-5, Hex-5, Fuc-1"
}

for glycan_code in common_n_glycans:
    gp = Glycopeptide(
        peptide_sequence=peptide,
        glycan_code=glycan_code,
        glycosylation_site=2,
        glycan_type="N"
    )
    
    fragments = gp.generate_fragments()
    
    # Get Y1 and Intact masses
    y1 = [f for f in fragments if f['type'] == 'Y1']
    intact = [f for f in fragments if f['type'] == 'Intact']
    
    y1_mass = y1[0]['mass'] if y1 else 0
    intact_mass = intact[0]['mass'] if intact else 0
    
    print(f"{glycan_code:<10} {compositions_info[glycan_code]:<30} {y1_mass:<12.2f} {intact_mass:<12.2f}")
```

---

## Advanced Techniques

### Advanced 1: Custom Fragment Analysis

```python
from glycofrag import Glycopeptide

gp = Glycopeptide(
    peptide_sequence="EPTETTQSVK",
    glycan_code="4501",
    glycosylation_site=2,
    glycan_type="N"
)

fragments = gp.generate_fragments()

# Analyze by fragment mass range
print("Fragment Distribution by Mass Range:")
ranges = [(0, 200), (200, 400), (400, 600), (600, 1000), (1000, 2000)]

for low, high in ranges:
    in_range = [f for f in fragments if low <= f['mass'] < high]
    print(f"  {low}-{high} m/z: {len(in_range)} fragments")

# Find highest mass fragment
highest = max(fragments, key=lambda x: x['mass'])
print(f"\nHighest m/z: {highest['name']} = {highest['mass']:.4f} m/z")

# Find lowest mass fragment > 100 m/z
lowest = min([f for f in fragments if f['mass'] > 100], key=lambda x: x['mass'])
print(f"Lowest m/z (>100): {lowest['name']} = {lowest['mass']:.4f} m/z")
```

### Advanced 2: Fragment Type Specific Analysis

```python
from glycofrag import Glycopeptide

gp = Glycopeptide(
    peptide_sequence="LCPDCPLLAPLNDSR",
    glycan_code="4501",
    glycosylation_site=2,
    glycan_type="N"
)

fragments = gp.generate_fragments()

# Analyze Y-glycan fragments (peptide + individual glycan fragments)
y_glycan = [f for f in fragments if f['type'] == 'Y-glycan']

print(f"Y-Glycan Fragments (Peptide + Glycan Fragment):")
print(f"Total: {len(y_glycan)}\n")

# Group by glycan fragment type
glycan_types = {}
for frag in y_glycan:
    name_parts = frag['name'].split('-')
    if len(name_parts) > 1:
        glycan_part = name_parts[1]
        glycan_types[glycan_part] = glycan_types.get(glycan_part, 0) + 1

print("Glycan fragment types in Y-glycan ions:")
for gtype in sorted(glycan_types.keys()):
    print(f"  {gtype}: {glycan_types[gtype]} ions")
```

### Advanced 3: Charge State Optimization

```python
from glycofrag import Glycopeptide

gp = Glycopeptide(
    peptide_sequence="LCPDCPLLAPLNDSR",
    glycan_code="4501",
    glycosylation_site=2,
    glycan_type="N"
)

fragments = gp.generate_fragments()

# Find optimal charge state for different mass ranges
print("Optimal charge state for different mass ranges:\n")
print(f"{'Mass Range':<15} {'Preferred z':<12} {'Common m/z':<12}")
print("-" * 40)

ranges = [(0, 300), (300, 600), (600, 1000), (1000, 2000)]

for low, high in ranges:
    in_range = [f for f in fragments if low <= f['mass'] < high]
    if in_range:
        charge_dist = {}
        for f in in_range:
            z = f['charge']
            charge_dist[z] = charge_dist.get(z, 0) + 1
        
        best_z = max(charge_dist.keys(), key=lambda x: charge_dist[x])
        sample = [f for f in in_range if f['charge'] == best_z][0]
        
        print(f"{low}-{high} m/z{15-len(str(low))-len(str(high))} z={best_z:<11} {sample['mass']:<12.2f}")
```

---

## Tips and Best Practices

1. **Memory Efficiency**: For large batch processing, iterate rather than storing all fragments
2. **Fragment Filtering**: Use fragment type to focus analysis (e.g., only Y1 for glycan confirmation)
3. **Error Handling**: Always check fragment lists are non-empty before accessing elements
4. **Mass Tolerance**: Use ±50 ppm when matching experimental to theoretical
5. **Charge States**: Consider typical MS instrument charge distributions (usually z=1-3 for glycopeptides)
6. **Modifications**: Always specify expected modifications for accurate mass matching
7. **Glycan Position**: Verify glycan position matches actual sequencing (use 0-indexed)

---

**Last Updated**: February 1, 2026  
**Glycofrag Version**: 0.1.0
