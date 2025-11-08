# LigX
LigX is a research tool designed for ligand generation and extension within the binding site of a target molecule, following a fragment-based approach.

## Highlights of the LigX

### 1. Synthetic Accessibility
During the fragment-linking process, LigX incorporates the concept of retrosynthetic analysis to ensure the synthetic accessibility of the generated molecules.

### 2. Structural Diversity
LigX employs a breadth-first search strategy to systematically explore a wider range of chemical space.

### 3. Generating Molecules as 3D Structures
LigX generate multiple compounds as a *xyz file* format. This ensures that the generated structure satisfies the structural requirements of the target protein. Furthermore, users can conduct subsequent more reliable evaluations (e.g. docking simulation and molecular dynamics simulation) based on this structure, allowing them to identify compounds that show both high affinity and synthetic feasibility from a large pool of candidates.


## Requirements

We tested on the python 3.9.23. It requires to install several modules. 
- **Required**
  - `numpy (2.0.2)`
- **Recommended**
  - `rdkit` — **must be a build that includes the `rdkit.Chem.rdDetermineBonds` module**.  
    Example of a verified build: `rdkit=2025.03.5=py39h364ec62_0` (from `conda-forge`).  
    > If a version of RDkit that satisfies the requirements is not installed, SMILES filter will be unavailable. LigX automatically switches the SMILES filter depending on the environment. Even if the SMILES filter is unavailable, compound generation in xyz format can be performed.

> **Note on portability**  
> Pinning to `rdkit=2025.03.5` ensures reproducibility on the tested macOS setup.  
> For cross-platform use, you may relax the build pin to any version that ships `rdkit.Chem.rdDetermineBonds`.

---

## Installation

### pip and conda
```
conda create -n ligx python=3.9 numpy pip -c conda-forge
conda activate ligx
pip install git+https://github.com/akatsuyama/LigX.git
```
Installation of RDkit (recommended)
```
conda install -c conda-forge rdkit
```


## Preparation of Input Files
LigX uses xyz coordinate of protein and ligand (initiator) molecules. Additionally, a text file describing integer values of charge (charge.txt) must be provided in the same directory.

### Input file layout 

Directorys containing input files are placed at the same level as `run_ligx.py`. 

<pre>
work_dir/
├── run_ligx.py
├── ligand_xyz/
│   └── ligand_name/
│       ├── ligand_name.xyz
│       └── charge.txt
├── protein_xyz/
│   └── protein_name.xyz
├── fragment/ (optional)
│   ├── ***
│   ├── .
│   └── .    
└── rollout_fragment/ (optional)
    ├── ***
    ├── .
    └── .
</pre>

The directory names `ligand_xyz/` and `protein_xyz/` are specified as LigX default values. The results will be output to `work_dir`.


## Quick Start

Please prepare a python script to run LigX program (example: `run_ligx.py`)

```run_ligx.py
#example

#import LigX package and class LigXRun
from LigX import LigXRun

input_file = LigXRun(linknum=4)
input_file.exec_ligx()
```
### Running

```
python run_ligx.py
```

### Output Files

Output files in LigX are stored in `result/`. Log files are stored in current directory. `result_all_total_conf_(number).xyz` is suitable for the docking simulations. `image_ligx/` stores the chemical structure of the generated compound (SMIELES filter required). This file can be visualized using programs such as PyMOL. `result_score.csv` contains SMILES and interaction energies roughlly estimated by LigX. *(Using this value for evaluation of compounds is not recommended. Please use other reliable docking tools for affinity evaluation.)*

<pre>
work_dir/
├── run_ligx.py
├── .
├── .
├── .
├── logfile.log
├── paramfile.log
└── result/
    ├── image_ligx/
    │   ├── chemical_structure1.png
    │   └── .
    ├── result_all_total_conf_(number).xyz  #suitable for the following docking simulations.
    ├── result_conf_(number)/
    │   ├── sorted_structures1/
    │   ├── .
    │   ├── .
    │   ├── sorted_structures.xyz
    │   ├── .
    │   └── .
    └── result_score.csv  #contains interaction energy and SMILES
</pre>

## Documentation

[LigX User Guide](https://akatsuyama.github.io/LigX)

## License

MIT License — see LICENSE for details.<br>
Third-party libraries are governed by their respective licenses (RDKit is BSD).

## Author

- Akira Katsuyama

## Contributor
- Kazuki Hoshi

