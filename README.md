# Loop Reconstruction for Protein Structures

This repository contains a Python script to identify and model missing loop regions in protein structures from `PDBx/mmCIF` files. The script processes the PDB entry 3IDP, identifies missing loops, generates multiple conformations for each loop, and saves the results as new `.cif` files with metadata.

## Overview

Protein structures in the Protein Data Bank (PDB) often have missing loop regions due to experimental limitations. This script:

- Parses a `PDBx/mmCIF` file (e.g., `3IDP.cif`) using `BioPython`.
- Identifies inner missing loops by detecting gaps in residue sequences.
- Models missing loops using linear interpolation between anchor C-alpha atoms with random perturbations.
- Generates 5 conformations per loop, checks for steric clashes, and assigns a quality score (0.8 for clash-free, 0.2 for clashes).
- Saves modeled structures as `.cif` files and metadata in a text file.

For `3IDP` included in the `data` directory of the current repository, it models three loops:

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
   ```

2. **Create a Virtual Environment**:

   ```bash
   python3 -m venv Myenv
   source Myenv/bin/activate  # On Windows: Myenv\Scripts\activate
   ```

3. **Install Dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

4. **Obtain Input File**:
   - `3IDP.cif` is placed in the `data/` directory.
   - Alternatively, download it from [RCSB PDB](https://files.rcsb.org/download/3IDP.cif).

## Usage

Run the script from the repository root:

```bash
cd src
python3 loop_reconstruction.py
```

The script will:

- Read `data/3IDP.cif`.
- Generate 15 modeled structures (5 conformations per loop).
- Save outputs in `data/`:
  - `3IDP_modeled_0.cif` to `3IDP_modeled_14.cif`.
  - `3IDP_metadata.txt` with details of modeled regions and scores.

## Output

- **Structural Files**: 15 `.cif` files in `data/` (e.g., `3IDP_modeled_0.cif` to `3IDP_modeled_14.cif`).
  - Each file contains the original structure with one modeled loop conformation.
  - Conformations 0–4: Chain B (598–613), 5–9: Chain A (600–614), 10–14: Chain A (629–630).
- **Metadata File**: `data/3IDP_metadata.txt` with entries like:

   ```bash
  Chain B, Residues 598-613, Conformation 0, Score 0.2
  Chain A, Residues 600-614, Conformation 0, Score 0.2
  ...
   ```

## Limitations

- Models only C-alpha atoms (no backbone N, C, O or side chains).
- Uses generic GLY residues instead of the true sequence.
- All conformations score 0.2 due to steric clashes (simplistic coordinate generation).
- Binary scoring (0.2/0.8) lacks nuance (e.g., no RMSD or geometry validation).
- Hardcoded for 3IDP; lacks command-line input for other PDB files.

## Future Improvements

- Add backbone atoms and correct residue types from `_entity_poly_seq`.
- Use Ramachandran angles for realistic loop geometry.
- Enhance scoring with RMSD, precise clash detection (`NeighborSearch`), and geometry checks.
- Support command-line arguments for PDB ID and max structures.
- Validate models with MolProbity or molecular dynamics.
- Test on additional structures (e.g., 6X18, 8RX0).

## License

MIT License. See [LICENSE](LICENSE) for details.
