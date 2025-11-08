# Special atom names used in LigX

Special atom names are used in an initiator and fragments to ensure synthetic accessibility and remove unsuitable fragment linking.
Atom names are categorized based on **skeletal connector (SC)** and **functional group connector (FC)**. All of the atoms bound to the SC are carbon atoms.

## 1. Atom names and properties bound for SCs.

| atom name | description |
|---|---|
| XC3  | A sp3 carbon atom except for methylene type (XC3g) carbon atom. |
| XC2  | A sp2 carbon atom except for carbonyl carbon atom. |
| XC1  | A sp carbon atom of an alkyne. |
| XC3g | A sp3 carbon atom containing vicinal connectors such as methylene group. |
| XCA  | A sp2 carbon atom of aromatic ring. |

*Note*: The above atom names are actually used with an underscore and numbers appended. For details, refer to [Initiator Preparation](initiator.md)

## 2. Atom names and properties bound for FCs.

| atom name | description |
|---|---|
| YCA | A carbonyl carbon atom. This atom accepts many fragments, allowing flexible utilization not only in carbonyl carbon but also in various substituents. |
| YCO | A carbon atom in the carboxylate fragment. |
| YO  | An oxygen atom in an ether moiety. |
| YN3 | An ammonia type of nitrogen atom (positively charged). |
| YN2 | An aniline type of nitrogen atom (neutral nitrogen atom). |
| YNA | An amide type of nitrogen atom (neutral nitrogen atom). This atom accepts many fragments, allowing flexible utilization not only in amide nitrogen but also in various substituents.|
| YS  | A sulfur atom in a thioether moiety. |
| YSA | A sulfur atom in a sulfonamide moiety. |

*Note*: The above atom names are used *without* an underscore and numbers.

## 3. Rule of the linking using special atom names

In LigX program, the fragment linking which would generate chemically unstable or less synthetically accessible compounds can be removed based on the atom names between the fragments. During compound generation, the special atom names of the initiator and the used fragments are recorded. The linking permission can be determined by the combination between the atom name of each fragment and previous one or two atom recorded. Examples of combinations excluded by AAA are shown below.

### a. Two-atoms combination
Two-atom combination listed in `LigX.LigX_core.distance_dict` defines the bond length between these atoms. Additionally, If the combination is not defined, the linking between this combination is not permitted. For example, a fragment linking between a sp2 carbon atom (XC2) and a nitrogen atom (YN2) generates enamine moiety that can be hydrolyzed in the presence of water. To eliminate the chance for yielding molecule with enamine moiety, the combination of XC2 and YN2 is not defined in the LigX program.

### b. Three-atoms combination
In certain cases, chemically unstable structures can be generated through the two specific linking. For example, a sequential connection of a nitrogen atom (YN2), methylene group (XC3g), and an oxygen atom (YO) generates an hemiaminal moiety that can be hydrolyzed easily. To remove such structures, LigX defines the unacceptable sequence of three-atoms combination in `LigX.LigX_core.not_permitted_comb_dict`.

## 4. Special atom names only used for the positional reference
Fragment linking requires another special atom. The atom name is as shown in the table below.<br>
For details, refer to [Initiator Preparation](initiator.md)

| atom name | description |
|---|---|
| CX_0  | A carbon atom for the positional reference |
| NX_0  | A carbon atom for the positional reference |
| RH_2  | A hydrogen atom for the positional reference |
| RC_2  | A carbon atom for the positional reference |
| RO_2  | A oxygene atom for the positional reference |
| RN_2  | A nitrogen atom for the positional reference |
| RS_2  | A sulfur atom for the positional reference |
| RCl_2  | A chlorine atom for the positional reference |
| RBr_2  | A bromine atom for the positional reference |
| RI_2  | A iodine atom for the positional reference |
| RF_2  | A fluorine atom for the positional reference |
| RP_2  | A phosphorus atom for the positional reference |
