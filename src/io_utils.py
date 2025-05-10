from Bio.PDB.mmcifio import MMCIFIO

def save_structure(structure, output_file):
    """Save structure to PDBx/mmCIF file."""
    print(f"Saving structure to {output_file}")
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)

def save_metadata(metadata, output_file):
    """Save metadata to a text file."""
    with open(output_file, "w") as f:
        for entry in metadata:
            f.write(f"Chain {entry['chain']}, Residues {entry['start_res']}-{entry['end_res']}, "
                    f"Conformation {entry['conformation']}, Score {entry['score']:.2f}\n")
    print(f"Saved metadata to {output_file}")
