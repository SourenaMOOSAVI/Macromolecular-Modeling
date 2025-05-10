import numpy as np
from Bio.PDB import Residue, Atom
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def get_residue_types(cif_file, chain_id, start_num, gap_size):
    """Parse residue types from _entity_poly_seq in mmCIF."""
    mmcif_dict = MMCIF2Dict(cif_file)
    seq = mmcif_dict.get("_entity_poly_seq", [])
    chain_map = {'A': '1', 'B': '2'}  # Adjust based on 3IDP's entity IDs
    entity_id = chain_map.get(chain_id, '1')
    residue_types = []
    for i in range(start_num, start_num + gap_size):
        res_num = str(i)
        res_type = "GLY"  # Default
        for entry in seq:
            if entry.get("entity_id") == entity_id and entry.get("num") == res_num:
                res_type = entry.get("mon_id", "GLY")
                break
        residue_types.append(res_type)
    return residue_types

def generate_loop_coordinates(start_coord, end_coord, gap_size, num_conformations=5):
    """Generate loop coordinates using simplified Ramachandran angles."""
    rng = np.random.default_rng()
    conformations = []
    # Simplified Ramachandran angles (phi, psi) for alpha-helix-like loops
    ramachandran_angles = [(-60, -40)] * gap_size  # Placeholder; can use a library
    ca_distance = 3.8  # Approximate C-alpha distance

    for _ in range(num_conformations):
        coords = []
        current_coord = start_coord
        for i in range(gap_size):
            # Simulate peptide bond geometry (simplified)
            phi, psi = ramachandran_angles[i]
            # Approximate next C-alpha position
            displacement = np.array([ca_distance * np.cos(np.radians(phi)),
                                   ca_distance * np.sin(np.radians(phi)),
                                   0.0])
            perturbation = rng.normal(0, 0.5, 3)
            next_coord = current_coord + displacement + perturbation
            coords.append(next_coord)
            current_coord = next_coord
        # Adjust final coordinate to approach end_coord
        if coords:
            coords[-1] = 0.5 * (coords[-1] + end_coord)  # Smooth connection
        conformations.append(coords)
    print(f"Generated {len(coords)} coordinates for gap_size {gap_size}")
    return conformations

def create_loop_residues(chain, start_num, gap_size, loop_coords, residue_types):
    """Create residues with backbone atoms for the modeled loop."""
    residues = []
    for i in range(gap_size):
        if i >= len(loop_coords):
            print(f"Warning: Not enough coordinates ({len(loop_coords)}) for gap_size {gap_size}")
            break
        res_id = (' ', start_num + i, ' ')
        res_name = residue_types[i]
        residue = Residue(res_id, res_name, ' ')
        ca_coord = loop_coords[i]
        # Add backbone atoms with simplified geometry
        n_coord = ca_coord + np.array([1.45, 0, 0])  # N-Ca bond ~1.45 Å
        c_coord = ca_coord + np.array([-1.33, 0, 0])  # Ca-C bond ~1.53 Å
        o_coord = c_coord + np.array([0, 1.24, 0])   # C=O bond ~1.24 Å
        residue.add(Atom('CA', ca_coord, 10.0, 1.0, ' ', 'CA', serial_number=1, element='C'))
        residue.add(Atom('N', n_coord, 10.0, 1.0, ' ', 'N', serial_number=2, element='N'))
        residue.add(Atom('C', c_coord, 10.0, 1.0, ' ', 'C', serial_number=3, element='C'))
        residue.add(Atom('O', o_coord, 10.0, 1.0, ' ', 'O', serial_number=4, element='O'))
        residues.append(residue)
    print(f"Created {len(residues)} residues")
    return residues
