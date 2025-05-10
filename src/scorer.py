import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch

def calculate_rmsd(start_ca, end_ca, loop_coords):
    """Calculate RMSD of loop anchors to original coordinates."""
    if not loop_coords:
        return 5.0  # Arbitrary high RMSD for invalid loops
    start_rmsd = np.linalg.norm(loop_coords[0] - start_ca)
    end_rmsd = np.linalg.norm(loop_coords[-1] - end_ca)
    return np.sqrt((start_rmsd**2 + end_rmsd**2) / 2)

def validate_geometry(loop_coords):
    """Check C-alpha distances (~3.8 Å)."""
    if len(loop_coords) < 2:
        return 0.0
    distances = [np.linalg.norm(loop_coords[i] - loop_coords[i+1]) for i in range(len(loop_coords)-1)]
    valid = sum(1 for d in distances if 3.5 <= d <= 4.1) / len(distances)
    return valid

def check_steric_clashes(structure, loop_coords, threshold=2.0):
    """Check for steric clashes using NeighborSearch."""
    atoms = [atom for model in structure for chain in model for res in chain for atom in res]
    searcher = NeighborSearch(atoms)
    clash_count = sum(1 for coord in loop_coords if searcher.search(coord, threshold, level="A"))
    clash_score = 1.0 - min(clash_count / 10, 1.0)  # Normalize clash penalty
    return clash_score

def score_conformation(structure, start_ca, end_ca, loop_coords):
    """Combine RMSD, clash, and geometry scores."""
    rmsd = calculate_rmsd(start_ca, end_ca, loop_coords)
    clash_score = check_steric_clashes(structure, loop_coords)
    geometry_score = validate_geometry(loop_coords)
    # Composite score: Lower RMSD, fewer clashes, better geometry
    score = 0.4 * (1 - min(rmsd / 5, 1.0)) + 0.4 * clash_score + 0.2 * geometry_score
    return max(0.1, min(1.0, score))  # Clamp between 0.1 and 1.0
