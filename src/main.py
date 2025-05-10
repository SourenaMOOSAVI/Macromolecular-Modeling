import argparse
import numpy as np
import os
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from loop_finder import parse_structure, find_missing_loops
from loop_modeler import generate_loop_coordinates, create_loop_residues, get_residue_types
from scorer import score_conformation
from io_utils import save_structure, save_metadata

def main():
    parser = argparse.ArgumentParser(description="Reconstruct missing loops in protein structures.")
    parser.add_argument("--pdb_id", default="3IDP", help="PDB ID of the structure")
    parser.add_argument("--cif_file", default="data/3IDP.cif", help="Path to mmCIF file")
    parser.add_argument("--max_structures", type=int, default=15, help="Maximum number of output structures")
    parser.add_argument("--num_conformations", type=int, default=5, help="Conformations per loop")
    args = parser.parse_args()

    try:
        structure = parse_structure(args.pdb_id, args.cif_file)
        missing_loops = find_missing_loops(structure)
        if not missing_loops:
            print(f"No missing loops found in {args.pdb_id}.")
            return

        print(f"Found {len(missing_loops)} missing loop(s) in {args.pdb_id}:")
        for region in missing_loops:
            print(f"Chain {region['chain']}: Missing residues {region['start_res']+1} to {region['end_res']-1} "
                  f"(size: {region['gap_size']})")

        output_structures = []
        metadata = []
        structures_generated = 0
        output_dir = "../data"
        os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

        for region in missing_loops:
            if structures_generated >= args.max_structures:
                break
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

            residue_types = get_residue_types(args.cif_file, chain_id, start_res + 1, gap_size)
            loop_conformations = generate_loop_coordinates(start_ca, end_ca, gap_size, args.num_conformations)

            for i, loop_coords in enumerate(loop_conformations):
                if structures_generated >= args.max_structures:
                    break
                print(f"Creating new structure for conformation {i}")
                new_structure = Structure(args.pdb_id + f"_loop_{structures_generated}")
                new_model = Model(0)
                new_structure.add(new_model)

                for chain in structure[0]:
                    new_chain = Chain(chain.id)
                    new_model.add(new_chain)
                    if chain.id == chain_id:
                        for res in chain:
                            if res.get_id()[1] <= start_res:
                                new_chain.add(res.copy())
                        loop_residues = create_loop_residues(new_chain, start_res + 1, gap_size, loop_coords, residue_types)
                        for res in loop_residues:
                            new_chain.add(res)
                        for res in chain:
                            if res.get_id()[1] >= end_res:
                                new_chain.add(res.copy())
                    else:
                        for res in chain:
                            new_chain.add(res.copy())

                score = score_conformation(structure, start_ca, end_ca, loop_coords)
                output_structures.append((new_structure, score))
                metadata.append({
                    'chain': chain_id,
                    'start_res': start_res + 1,
                    'end_res': end_res - 1,
                    'conformation': i,
                    'score': score
                })
                structures_generated += 1

        for i, (struct, score) in enumerate(output_structures):
            output_file = os.path.join(output_dir, f"{args.pdb_id}_modeled_{i}.cif")
            save_structure(struct, output_file)
            print(f"Saved structure {output_file} with score {score:.2f}")

        metadata_file = os.path.join(output_dir, f"{args.pdb_id}_metadata.txt")
        save_metadata(metadata, metadata_file)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
