# Using User-Specified Fragments

## Additional fragment set
Your can obtain additional fragments that can be used in LigX from the Github repository ([akatsuyama/LigX-data](https://github.com/akatsuyama/LigX-data)).

Each fragment dataset (an xyz file and two txt files) is stored in a directory, categorized as either C2L or F2L. Both fragment types have two connectors. In the case of C2L fragments, the SC connector is the one that binds to the initiator. In contrast, F2L fragments bind to the initiator through the FC connector.

Place any user-specified fragment directories into either the C2L or F2L directory, depending on their category, following the layout shown below.

The following is an example layout when using six fragments: xxx, yyy, zzz, aaa, bbb, and ccc.

<pre>
work_dir/
├── run_ligx.py
├── .
├── .
└── external_fragment/
    ├── C1L
    ├── C2L/
    │   ├── xxx/        #obtained fragment (directory)
    │   ├── yyy/        #obtained fragment (directory)
    │   └── zzz/        #obtained fragment (directory)
    ├── F1L
    └── F2L/
        ├── aaa/        #obtained fragment (directory)
        ├── bbb/        #obtained fragment (directory)
        └── ccc/        #obtained fragment (directory)

</pre>

An external_fragment directory containing empty C2L and F2L subfolders can be downloaded from the GitHub repository. Other directories, such as C1L, may already contain required data set. It is recommended to use this template when organizing your custom fragments according to their respective categories.

*Note*: When an external fragment set is specified, the default fragment set normally loaded by LigX will not be included. This means that if you want to add fragments, you must combine your additional fragments with those from the default set. The default fragment set can be obtained from the GitHub repository ([akatsuyama/LigX-data](https://github.com/akatsuyama/LigX-data)).


## Default fragment set in LigX

The following fragment set is used by default.

### Default fragment set
- C2L fragments

| fragment name | chemical structure |
|---|---|
| B1L_0000 (SC)  | hydrogen (SC) |
| C2L_0001  | para-phenyl (SC to SC) |
| C2L_0002  | ethylene |
| C2L_0003  | methylene |
| C2L_0004  | meta-phenyl (SC to SC) |
| C2L_0005  | ortho-phenyl (SC to SC) |
| C2L_0032  | cyclopropane |
| C2L_0035  | 1,2,3-triazole (C to N) |
| C2L_0039  | pyridine (2,4) |
| C2L_0040  | pyridine (2,4) |
| C2L_0041  | pyridine (2,5) |
| C2L_0042  | pyridine (2,5) |
| C2L_0044  | para-phenyl (SC to FC) |
| C2L_0045  | meta-phenyl (SC to FC) |
| C2L_0046  | ortho-phenyl (SC to FC) |

- F2L fragments

| fragment name | chemical structure |
|---|---|
| AMIL2_0001  | ammonium nitrogen |
| AML2_0001_A1  | trans-amide (C to N) |
| AML2_0001_B1  | trans-amide (N to C) |
| ANIL2_0001  | aniline nitrogen |
| C2L_F_0035  | 1,2,3-triazole (N to C) |
| C2L_F_0044  | para-phenyl (FC to SC) |
| C2L_F_0045  | meta-phenyl (FC to SC) |
| C2L_F_0046  | ortho-phenyl (FC to SC) |
| F2L_0004  | piperazine (aniline N to ammonium N) |
| F2L_0005  | N-Me-amide (C to N) |
| F2L_0006  | N-Me-amide (N to C) |
| SAM_YN2_A  | sulfonamide (YN2, S to N) |
| SAM_YN2_B  | sulfonamide (YN2, N to S) |
| SAM_YN3_A  | sulfonamide (YN3, S to N) |
| SAM_YN3_B  | sulfonamide (YN3, N to S) |
| YO  | oxygene (ether) |

### Default rollout fragment set
- C2L fragments

| fragment name | chemical structure |
|---|---|
| B1L_0000 (SC)  | hydrogen (SC) |
| C2L_0001  | para-phenyl (SC to SC) |
| C2L_0003  | methylene |
| C2L_0005  | ortho-phenyl (SC to SC) |
| C2L_0044  | para-phenyl (SC to FC) |
| C2L_0046  | ortho-phenyl (SC to FC) |

- F2L fragments

| fragment name | chemical structure |
|---|---|
| AMIL2_0001  | ammonium nitrogen |
| AML2_0001_A1  | trans-amide (C to N) |
| AML2_0001_B1  | trans-amide (N to C) |
| C2L_F_0044  | para-phenyl (FC to SC) |
| C2L_F_0046  | ortho-phenyl (FC to SC) |
| YO  | oxygene (ether) |


### Fragment set at terminal
Regardless of the contents of the fragment set, the following fragments are used at the terminal.

| fragment name | chemical structure |
|---|---|
| B1L_0000 (SC)  | hydrogen (SC) |
| B1L_0000 (FC)  | hydrogen (FC) |
| carboxylate  | carboxylate |
