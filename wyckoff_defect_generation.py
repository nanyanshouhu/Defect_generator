import os
from pymatgen.core import Structure, Element, Site
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from itertools import product
import re
from collections import defaultdict

# Load the structure from the POSCAR file
structure = Structure.from_file("POSCAR")

# Initialize SpacegroupAnalyzer to analyze symmetry
analyzer = SpacegroupAnalyzer(structure)

# Get the Wyckoff position for each atom
wyckoff_positions = analyzer.get_symmetry_dataset()['wyckoffs']

# Helper function to sanitize filenames by removing invalid characters
def sanitize_filename(filename):
    return re.sub(r'[\\/*?:"<>|]', "_", filename)

# Function to remove one atom per element at specified Wyckoff positions
def remove_individual_wyckoff(structure, element_symbols):
    element_indices = {element_symbol: [] for element_symbol in element_symbols}
    for i, (site, wyckoff) in enumerate(zip(structure, wyckoff_positions)):
        if site.species_string in element_symbols:
            element_indices[site.species_string].append(i)

    for combination in product(*element_indices.values()):
        modified_structure = structure.copy()
        modified_structure.remove_sites(combination)
        folder_suffix = "_".join([f"{structure[i].species_string}_{wyckoff_positions[i]}" for i in combination])
        sanitized_folder = sanitize_filename(f"Removed_{folder_suffix}")
        os.makedirs(sanitized_folder, exist_ok=True)
        poscar_path = os.path.join(sanitized_folder, "POSCAR")
        modified_structure.to(filename=poscar_path)
        print(f"Created folder and file: {sanitized_folder}/POSCAR")

# Function to replace one atom per element at specified Wyckoff positions
def replace_individual_wyckoff(structure, element_symbols, new_elements):
    element_indices = {element_symbol: [] for element_symbol in element_symbols}
    for i, (site, wyckoff) in enumerate(zip(structure, wyckoff_positions)):
        if site.species_string in element_symbols:
            element_indices[site.species_string].append(i)

    for combination in product(*element_indices.values()):
        modified_structure = structure.copy()
        for idx, i in enumerate(combination):
            modified_structure.replace(i, Element(new_elements[idx]))
        folder_suffix = "_".join([f"{structure[i].species_string}_{wyckoff_positions[i]}_to_{new_elements[idx]}" for idx, i in enumerate(combination)])
        sanitized_folder = sanitize_filename(f"Replaced_{folder_suffix}")
        os.makedirs(sanitized_folder, exist_ok=True)
        poscar_path = os.path.join(sanitized_folder, "POSCAR")
        modified_structure.to(filename=poscar_path)
        print(f"Created folder and file: {sanitized_folder}/POSCAR")

# Function to insert multiple interstitial atoms at specified coordinates
def insert_interstitials(structure, elements_coords):
    """
    Insert interstitial atoms at specified fractional coordinates.
    Args:
        structure: pymatgen.core.Structure object.
        elements_coords: List of tuples, each containing the element symbol and fractional coordinates
                         (e.g., [("Li", [0.25, 0.25, 0.25]), ("Na", [0.5, 0.5, 0.5])]).
    """
    modified_structure = structure.copy()
    for element_symbol, coords in elements_coords:
        # Insert interstitial atoms at the specified fractional coordinates
        modified_structure.append(Element(element_symbol), coords, coords_are_cartesian=False)

    # Generate unique folder and filename based on interstitials
    coords_str = "_".join([f"{element}_{'_'.join([f'{c:.3f}' for c in coords])}" for element, coords in elements_coords])
    sanitized_folder = sanitize_filename(f"Inserted_{coords_str}")
    os.makedirs(sanitized_folder, exist_ok=True)
    poscar_path = os.path.join(sanitized_folder, "POSCAR")
    modified_structure.to(filename=poscar_path)
    print(f"Created folder and file: {sanitized_folder}/POSCAR")
# Function to generate antisite defects considering Wyckoff positions automatically
def generate_antisite_defects_multiple_pairs(structure, element_pairs):
    """
    Generate antisite defects by moving multiple elements into unique Wyckoff positions of other elements,
    considering each group of equivalent Wyckoff positions only once.
    
    Args:
        structure: pymatgen.core.Structure object.
        element_pairs: List of tuples specifying pairs of elements, where the first element moves into 
                       the unique Wyckoff positions of the second element (e.g., [("Sr", "Ti"), ("La", "O")]).
    """
    element_indices = {element: [] for pair in element_pairs for element in pair}
    
    # Gather indices and Wyckoff positions for each specified element
    for i, (site, wyckoff) in enumerate(zip(structure, wyckoff_positions)):
        if site.species_string in element_indices:
            element_indices[site.species_string].append((i, wyckoff))
    
    for element1, element2 in element_pairs:
        # Process combinations of indices for the elements involved
        indices_element1 = element_indices[element1]
        indices_element2 = element_indices[element2]
        
        for idx1, idx2 in product(indices_element1, indices_element2):
            index1, wyckoff1 = idx1
            index2, wyckoff2 = idx2
            
            # Create a modified structure
            modified_structure = structure.copy()
            
            # Replace element2 site with element1
            modified_structure[index2] = Element(element1)
            
            # Sort the structure for consistency
            modified_structure.sort()
            
            # Generate folder name and save the modified structure
            folder_suffix = f"{element1}_into_{element2}_{wyckoff2}"
            sanitized_folder = sanitize_filename(f"Antisite_{folder_suffix}")
            os.makedirs(sanitized_folder, exist_ok=True)
            poscar_path = os.path.join(sanitized_folder, "POSCAR")
            modified_structure.to(filename=poscar_path)
            print(f"Created folder and file: {sanitized_folder}/POSCAR")


# Function to generate antisite defects by swapping two atoms

# Example usage:
# Remove one atom per element
# remove_individual_wyckoff(structure, ["O", "La"])

# Replace one atom per element
# replace_individual_wyckoff(structure, ["O", "La"], ["N", "Ca"])

# Insert interstitials
# insert one atom per element
# insert_interstitials(structure, [("Li", [0.25, 0.25, 0.25]), ("Na", [0.5, 0.5, 0.5])])

# Example usage
# Specify multiple element pairs for generating antisite defects
# generate_antisite_defects_multiple_pairs(structure, [("Sr", "Ti"), ("La", "O")])
generate_antisite_defects_multiple_pairs(structure, [("O", "Sr")])

