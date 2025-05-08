import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import MMCIFIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import os
import Bio

def parse_structure(pdb_id, cif_file):
    """Parse a PDBx/mmCIF file and return the structure."""
    print(f"BioPython version: {Bio.__version__}")
    print(f"MMCIFParser: {MMCIFParser}")
    if not os.path.exists(cif_file):
        raise FileNotFoundError(f"File {cif_file} not found.")
    parser = MMCIFParser()
    print(f"Parser instance: {parser}")
    print(f"get_structure method: {parser.get_structure}")
    structure = parser.get_structure(pdb_id, cif_file)
    print(f"Structure parsed: {structure}")
    return structure

def find_missing_loops(structure):
    """Identify missing loop regions in the structure."""
    missing_regions = []
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            for i in range(len(residues) - 1):
                current_res = residues[i]
                next_res = residues[i + 1]
                current_num = current_res.get_id()[1]
                next_num = next_res.get_id()[1]
                if next_num > current_num + 1 and is_aa(current_res) and is_aa(next_res):
                    missing_regions.append({
                        'chain': chain.id,
                        'start_res': current_num,
                        'end_res': next_num,
                        'gap_size': next_num - current_num - 1
                    })
    return missing_regions

def generate_loop_coordinates(start_coord, end_coord, gap_size, num_conformations=5):
    """Generate simple loop coordinates for a missing region."""
    print(f"NumPy version: {np.__version__}")
    rng = np.random.default_rng()
    print(f"rng.normal: {rng.normal}")
    ca_distance = 3.8
    conformations = []
    
    for _ in range(num_conformations):
        coords = []
        step = (end_coord - start_coord) / (gap_size + 1)
        for i in range(1, gap_size + 1):
            coord = start_coord + step * i
            perturbation = rng.normal(0, 0.5, 3)
            print(f"Perturbation: {perturbation}")
            coord += perturbation
            coords.append(coord)
        conformations.append(coords)
    print(f"Generated {len(coords)} coordinates for gap_size {gap_size}")
    return conformations

def check_steric_clashes(structure, loop_coords, threshold=2.0):
    """Check for steric clashes (simplified: distance-based)."""
    print(f"Entering check_steric_clashes")
    print(f"np.array: {np.array}")
    print(f"np.linalg.norm: {np.linalg.norm}")
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res:
                    atom_coord = np.array(atom.get_coord())
                    for loop_coord in loop_coords:
                        distance = np.linalg.norm(atom_coord - loop_coord)
                        if distance < threshold:
                            return False
    return True

def create_loop_residues(chain, start_num, gap_size, loop_coords, res_name="GLY"):
    """Create residues for the modeled loop."""
    print(f"Entering create_loop_residues with start_num={start_num}, gap_size={gap_size}, len(loop_coords)={len(loop_coords)}")
    residues = []
    for i in range(gap_size):
        if i >= len(loop_coords):
            print(f"Warning: Not enough coordinates ({len(loop_coords)}) for gap_size {gap_size}")
            break
        res_id = (' ', start_num + i, ' ')
        print(f"Creating residue with ID: {res_id}")
        residue = Residue(res_id, res_name, ' ')
        atom = Atom('CA', loop_coords[i], 10.0, 1.0, ' ', 'CA', serial_number=1, element='C')
        residue.add(atom)
        residues.append(residue)
    print(f"Created {len(residues)} residues")
    return residues

def save_structure(structure, output_file):
    """Save structure to PDBx/mmCIF file."""
    print(f"Saving structure to {output_file}")
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)

def main():
    pdb_id = "3IDP"
    cif_file = "3IDP.cif"
    num_conformations = 5
    
    try:
        structure = parse_structure(pdb_id, cif_file)
        missing_loops = find_missing_loops(structure)
        if not missing_loops:
            print(f"No missing loops found in {pdb_id}.")
            return
        
        print(f"Found {len(missing_loops)} missing loop(s) in {pdb_id}:")
        for region in missing_loops:
            print(f"Chain {region['chain']}: Missing residues {region['start_res']+1} to {region['end_res']-1} "
                  f"(size: {region['gap_size']})")
        
        output_structures = []
        metadata = []
        for region in missing_loops:
            chain_id = region['chain']
            start_res = region['start_res']
            end_res = region['end_res']
            gap_size = region['gap_size']
            
            chain = structure[0][chain_id]
            start_ca = None
            end_ca = None
            for res in chain:
                if res.get_id()[1] == start_res and 'CA' in res:
                    start_ca = np.array(res['CA'].get_coord())
                if res.get_id()[1] == end_res and 'CA' in res:
                    end_ca = np.array(res['CA'].get_coord())
            
            if start_ca is None or end_ca is None:
                print(f"Skipping loop in Chain {chain_id}: Anchor atoms not found.")
                continue
            
            loop_conformations = generate_loop_coordinates(start_ca, end_ca, gap_size, num_conformations)
            
            for i, loop_coords in enumerate(loop_conformations):
                print(f"Creating new structure for conformation {i}")
                print(f"Structure class: {Structure}")
                print(f"Structure module: {Structure.__module__}")
                new_structure = Structure(pdb_id + f"_loop_{i}")
                new_model = Model(0)
                new_structure.add(new_model)
                
                # Copy chains, inserting loop residues at the correct position
                for chain in structure[0]:
                    new_chain = Chain(chain.id)
                    new_model.add(new_chain)
                    if chain.id == chain_id:
                        # Copy residues strictly before the loop
                        print(f"Copying residues before loop for Chain {chain_id} (res <= {start_res})")
                        for res in chain:
                            if res.get_id()[1] <= start_res:
                                new_chain.add(res.copy())
                        # Debug: List residues in new_chain before adding loop
                        print(f"Residues in new_chain before loop: {[res.get_id() for res in new_chain]}")
                        # Add loop residues
                        loop_residues = create_loop_residues(new_chain, start_res + 1, gap_size, loop_coords)
                        for res in loop_residues:
                            print(f"Adding loop residue: {res.get_id()}")
                            new_chain.add(res)
                        # Debug: List residues in new_chain after adding loop
                        print(f"Residues in new_chain after loop: {[res.get_id() for res in new_chain]}")
                        # Copy residues after the loop
                        print(f"Copying residues after loop for Chain {chain_id} (res >= {end_res})")
                        for res in chain:
                            if res.get_id()[1] >= end_res:
                                print(f"Adding residue after loop: {res.get_id()}")
                                new_chain.add(res.copy())
                        # Debug: List final residues in new_chain
                        print(f"Final residues in new_chain: {[res.get_id() for res in new_chain]}")
                    else:
                        # Copy all residues for other chains
                        for res in chain:
                            new_chain.add(res.copy())
                
                clash_free = check_steric_clashes(structure, loop_coords)
                score = 0.8 if clash_free else 0.2
                
                output_structures.append((new_structure, score))
                metadata.append({
                    'chain': chain_id,
                    'start_res': start_res + 1,
                    'end_res': end_res - 1,
                    'conformation': i,
                    'score': score
                })
        
        for i, (struct, score) in enumerate(output_structures):
            output_file = f"3IDP_modeled_{i}.cif"
            save_structure(struct, output_file)
            print(f"Saved structure {output_file} with score {score}")
        
        with open("3IDP_metadata.txt", "w") as f:
            for entry in metadata:
                f.write(f"Chain {entry['chain']}, Residues {entry['start_res']}-{entry['end_res']}, "
                        f"Conformation {entry['conformation']}, Score {entry['score']}\n")
        print("Saved metadata to 3IDP_metadata.txt")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()