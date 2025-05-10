# Loop Reconstruction for Protein Structures

This repository contains a Python toolkit to identify and model missing loop regions in protein structures from PDBx/mmCIF files, addressing the assignment "Reconstruction of Missing Loops." Version 2 enhances the original script with modular code, realistic loop modeling, and improved scoring.

## Overview

Protein structures in the Protein Data Bank (PDB) often have missing loop regions due to experimental limitations. This toolkit:

- Parses a PDBx/mmCIF file (e.g., `3IDP.cif`) using BioPython.
- Identifies inner missing loops by detecting gaps in residue sequences.
- Models loops with backbone atoms (N, C, O) and correct residue types, using simplified Ramachandran angles.
- Generates multiple conformations per loop, scored by RMSD, steric clashes, and geometry.
- Saves modeled structures as `.cif` files and metadata in a text file.

For 3IDP, it models three loops:

- Chain B: Residues 598–613 (16 residues).
- Chain A: Residues 600–614 (15 residues).
- Chain A: Residues 629–630 (2 residues).

## Requirements

- Python 3.12
- BioPython 1.85
- NumPy 1.26.4

## Setup

1. **Clone the Repository**:

   ```bash
   git clone git@github.com:SourenaMOOSAVI/Macromolecular-Modeling.git
   cd Macromolecular-Modeling
   git checkout enhancement
   ```

2. **Create a Virtual Environment**:

   ```bash
   python3 -m venv Myenv
   source Myenv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install Dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

4. **Obtain Input File**:
   - Place `3IDP.cif` in the `data/` directory.
   - Alternatively, download it from [RCSB PDB](https://www.rcsb.org/structure/3IDP).

## Usage

Run the main script with command-line arguments:

```bash
cd src
python3 main.py --pdb_id 3IDP --cif_file ../data/3IDP.cif --max_structures 15 --num_conformations 5
```

**Arguments**:

- `--pdb_id`: PDB ID (default: `3IDP`).
- `--cif_file`: Path to mmCIF file (default: `data/3IDP.cif`).
- `--max_structures`: Maximum number of output structures (default: 15).
- `--num_conformations`: Conformations per loop (default: 5).

The script will:

- Read the specified `.cif` file.
- Generate up to `max_structures` modeled structures.
- Save outputs in `data/`:
  - `<pdb_id>_modeled_0.cif` to `<pdb_id>_modeled_N.cif`.
  - `<pdb_id>_metadata.txt` with details of modeled regions and scores.

## Output

- **Structural Files**: `.cif` files in `data/` (e.g., `3IDP_modeled_0.cif`).
  - Each file contains the original structure with one modeled loop conformation (backbone atoms, correct residues).
- **Metadata File**: `data/3IDP_metadata.txt` with entries like:

  ```plaintext
  Chain B, Residues 598-613, Conformation 0, Score 0.75
  Chain A, Residues 600-614, Conformation 0, Score 0.82
  ...
  ```

## Improvements in enhancement branch

- **Modular Structure**: Split into `loop_finder`, `loop_modeler`, `scorer`, `io_utils`, and `main`.
- **Physical Plausibility**:
  - Backbone atoms (N, C, O) for modeled residues.
  - Correct residue types from `_entity_poly_seq`.
  - Simplified Ramachandran angles for loop geometry.
- **Scoring**: Combines RMSD, precise clash detection (`NeighborSearch`), and C-alpha distance validation.
- **Input Flexibility**: Command-line arguments for PDB ID, file, and structure limits.

## Limitations

- Ramachandran angles are simplified (fixed alpha-helix-like); could use a full distribution.
- No side-chain atoms (focus on backbone).
- Scoring could include energy-based metrics or MolProbity validation.

## Future Work

- Use full Ramachandran distributions or loop libraries.
- Add side-chain atoms.
- Validate with MolProbity or molecular dynamics.
- Test on additional structures (e.g., 6X18, 8RX0).

## License

MIT License. See [LICENSE](LICENSE) for details.
