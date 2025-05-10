from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
import os

def parse_structure(pdb_id, cif_file):
    """Parse a PDBx/mmCIF file and return the structure."""
    if not os.path.exists(cif_file):
        raise FileNotFoundError(f"File {cif_file} not found.")
    parser = MMCIFParser()
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
    print(f"Found {len(missing_regions)} missing loop(s): {missing_regions}")
    return missing_regions
