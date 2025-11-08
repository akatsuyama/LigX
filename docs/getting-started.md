# Getting Started

## Requirements

This program is cross-platform and can be used on MacOS, Windows, and Ubuntu operating systems.
We tested on the python 3.9.23. It requires to install several modules. 

- **Required**
  - `numpy (2.0.2)`
- **Recommended**
  - `rdkit` — **must be a build that includes the `rdkit.Chem.rdDetermineBonds` module**.  
    Example of a verified build: `rdkit=2025.03.5=py39h364ec62_0` (from `conda-forge`).  
    If a version of RDkit that satisfies the requirements is not installed, SMILES filter will be unavailable. LigX automatically switches the SMILES filter depending on the environment. Even if the SMILES filter is unavailable, compound generation in xyz format can be performed.

> **Notes**  
> Pinning to `rdkit=2025.03.5` ensures reproducibility.  
> You may relax the build pin to any version that ships `rdkit.Chem.rdDetermineBonds`.

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

## Quick Start

Please prepare a python script to run LigX program (example: `run_ligx.py`).<br> 
Key paratemers and attributes are summarized in [LigX Functions-LigXRun](resources.md).
It is recommended that key parameters **linknum** and **outnum** are specified in the python script.

```run_ligx.py
#example (linknum = 4 and outnum = 3000)

#import LigX package and class LigXRun
from LigX import LigXRun

input_file = LigXRun(linknum=4,outnum=3000)
input_file.exec_ligx()
```

Please refer to [Data Layout](data-layout.md) and place the input files in their designated directories accordingly.

### Running
```
python run_ligx.py
```
For output structure and file locations, please refer to [Data Layout](data-layout.md).