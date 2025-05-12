# MCS-Align
## Overview

The codebase consists of two main components:

1. **JobManager.ipynb**: A Jupyter notebook that:
   - Processes Gaussian output files to identify OH vibrations in molecule-water clusters
   - Performs PCA analysis on water displacement matrices
   - Visualizes the results using RDKit and Bokeh
   - Generates structure files for further analysis

2. **MCSAlign**: A Python package that:
   - Finds common molecular substructures using Maximum Common Substructure (MCS) alignment
   - Generates new molecular configurations based on QM9 database matches
   - Prepares input files for extended dataset generation

## Key Features

- **Vibrational Analysis**: Identifies and analyzes OH vibrations in the frequency range 3400-3780 cm⁻¹
- **Structure Generation**: Uses RDKit and MCS alignment to generate new molecular configurations
- **Dataset Enrichment**: Creates an extended dataset of molecule-water clusters
- **Visualization**: Provides PCA projections of water displacement matrices vs. frequencies
- **Quality Control**: Includes tools for structure selection and dataset validation

## Usage

### JobManager.ipynb

The notebook processes Gaussian output files and generates visualizations:

1. Reads and processes `.out` files from various directories
2. Analyzes vibrational frequencies and displacements
3. Performs PCA on water displacement matrices
4. Generates interactive visualizations of the results


## Dependencies

- RDKit
- NumPy
- Pandas
- Bokeh
- cclib
- py3Dmol
- scikit-learn
- Jupyter


## Installation
```bash
pip install -e .
```

## Usage

```python
from src.MCSAlign import MCSAlign
Initialize with output file and motif
mcsa = MCSAlign("path/to/output.out", "motif")
Generate new structures
mcsa.find_interaction_mol()
mcsa.set_smiles(smiles)
mcsa.set_mols(mols)
mcsa.set_similarities()
mcsa.find_matches(cutoff, n_values=20)
mcsa.find_MCSs()
```

## Directory Structure
```bash
.
├── JobManager.ipynb # Main analysis notebook
├── MCSAlign/ # Structure generation package
│ ├── generate.py # Generation script
│ ├── main.py # Main execution script
│ └── src/ # Source code
├── jobs/ # generated Gaussian input files
├── hydra/ # example dataset 
    ├── train/ # Training set
    └── test/ # Test set
```



## Reference

```bibtex
@article{HyDRA2023,
title = {The first HyDRA challenge for computational vibrational spectroscopy},
journal = {Phys. Chem. Chem. Phys.},
year = {2023},
volume = {25},
pages = {22089-22102},
doi = {10.1039/D3CP01216F}
}
```


## License

[Add your license information here]