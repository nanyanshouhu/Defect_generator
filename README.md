Automated Point Defect Generation Script
This repository contains a Python script designed to automate the generation of point defects in a crystal structure. The script uses a POSCAR file as input and outputs configurations with various defect types. It is particularly useful for materials scientists and computational chemists working on defect chemistry and related properties.

Features
Vacancy Defects: Automatically removes atoms from specified positions to create vacancies.
Doping Defects: Replaces atoms with dopant species to simulate doping.
Interstitial Defects: Identifies potential interstitial sites and adds atoms to these locations.
Antisite Defects: Swaps atoms between lattice sites to create antisite defects.
Requirements
Python 3.x
pymatgen library for structure manipulation and defect generation.


# Example usage:
# Remove one atom per element
remove_individual_wyckoff(structure, ["O", "La"])

# Replace one atom per element
replace_individual_wyckoff(structure, ["O", "La"])

# Insert interstitials
# insert one atom per element
insert_interstitials(structure, [("Li", [0.25, 0.25, 0.25]), ("Na", [0.5, 0.5, 0.5])])

# Example usage
# Specify multiple element pairs for generating antisite defects
# generate_antisite_defects_multiple_pairs(structure, [("Sr", "Ti"), ("La", "O")])
generate_antisite_defects_multiple_pairs(structure, [("O", "Sr")])
