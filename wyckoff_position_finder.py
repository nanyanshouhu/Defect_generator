from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from collections import Counter

# Load the structure from the POSCAR file
structure = Structure.from_file("POSCAR")

# Initialize SpacegroupAnalyzer to analyze symmetry
analyzer = SpacegroupAnalyzer(structure)

# Get space group symbol and number
space_group_symbol = analyzer.get_space_group_symbol()
space_group_number = analyzer.get_space_group_number()
print(f"Space Group Symbol: {space_group_symbol}")
print(f"Space Group Number: {space_group_number}")

# Get the list of symmetry operations
symmetry_operations = analyzer.get_symmetry_operations()
print("Symmetry Operations:")
for op in symmetry_operations:
    print(op)

# Get the Wyckoff position for each atom
wyckoff_positions = analyzer.get_symmetry_dataset()['wyckoffs']

# Print each atom's element symbol and its corresponding Wyckoff position
print("\nAtom Wyckoff Positions:")
for site, wyckoff in zip(structure, wyckoff_positions):
    print(f"Atom: {site.species_string}, Wyckoff Position: {wyckoff}")

# Get symmetry-equivalent positions
equivalent_sites = analyzer.get_symmetrized_structure().equivalent_sites
print("\nEquivalent Sites:")
for i, sites in enumerate(equivalent_sites):
    print(f"\nEquivalent Site Group {i + 1}:")
    for site in sites:
        print(f"  - {site.species_string} at {site.frac_coords}")

# Count Wyckoff positions for each species
wyckoff_count = {}
for site, wyckoff in zip(structure, wyckoff_positions):
    element = site.species_string
    if element not in wyckoff_count:
        wyckoff_count[element] = []
    wyckoff_count[element].append(wyckoff)

# Summarize Wyckoff positions with counts
print("\nElement-wise Wyckoff Positions:")
for element, positions in wyckoff_count.items():
    position_count = Counter(positions)
    wyckoff_summary = " ".join([f"{count}{wyckoff}" for wyckoff, count in position_count.items()])
    print(f"{element}: {wyckoff_summary}")
