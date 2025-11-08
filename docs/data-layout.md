# Data Layout

## Input file layout 

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

The directory names `ligand_xyz/` and `protein_xyz/` are specified as LigX default values. `ligand_name.xyz` is a xyz file of initiator prepared according to the instructions on [Initiator Preparation](input_details/initiator.md). The results will be output to `work_dir`. These name can be specified by Attributes of LigXRun.<br>

*Note*: The `fragment` and `rollout_fragment` in the above example illustrate the layout when using external fragments. These two are optional because LigX contains default fragment sets required for the execution. For details, refer to [Using User-Specified Fragments](user_fragment.md)].

## Output file layout

Output files in LigX are stored in `result/`. Log files are stored in current directory. `result_all_total_conf_(number).xyz` is suitable for the docking simulations. This file can be visualized using programs such as PyMOL. `result_score.csv` contains SMILES and interaction energies roughlly estimated by LigX. *(Using this value for evaluation of compounds is not recommended. Please use other reliable docking tools for affinity evaluation.)*<br>

After the termination, `result_conf_(number)/` will be created that strores all xyz files and charge parameters.

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

