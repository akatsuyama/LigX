# Score Calculation for Breadth-First Search Based on Interaction Energy

Physically based interaction between each fragment and proteins is calculated by Lennard-Jones potential and Coulomb potential as a following equation.

$$
E = \frac{ \left\{ \displaystyle \sum 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]
\;+\; \alpha \displaystyle \sum \frac{Q_1 \cdot Q_2}{4\pi \varepsilon_0\, r} \right\} }{L}
$$

- $\alpha$: hyper parameter  
- $r$: distance  
- $L$: number of heavy atoms

These interactions are calculated only for atoms within the cutoff value (default value is 10 Å). Except for hydrogen atoms, the uniform LJ parameters and charge parameters are applied. For hydrogen atoms, two types of parameters are applied by classifying them into polar hydrogen and nonpolar hydrogen. The values of hyperparameters and the availability of L can be modified in the input file. Fragments with high E-values undergo rollout, and their E-values are updated. The score for selection is calculated using the following equation, based on this E value and the number N, which is the number of times the fragment has been selected for a given initiator.

$$
\mathrm{score}
= \sqrt{2}\,\sqrt{\frac{\log n}{N} + \frac{1}{E_{\min}-\mathrm{cutoff}}}\,(E-\mathrm{cutoff})
$$

- \(E_{\min}\): The highest interaction among the possible fragments in this stage.
- \(\mathrm{cutoff}\): Cutoff value for the standardization.
- \(n\): Total selection number at each initiator.

## LJ and charge parameters used in LigX

 | Atom | σ (Å)        | ε (kcal/mol) | Q (e)   |
|:----:|:-------------|:-------------|:-------:|
| H (non polar)   | 2.38760856462 | 0.1882800    | 0.09    |
| H (polar)    | 0.40001352445 | 0.1924640    | 0.43    |
| C    | 3.58141284692 | 0.2343040    | 0.07    |
| N    | 3.52795892384 | 0.2928800    | -0.40   |
| O    | 3.02905564168 | 0.5020800    | -0.50   |
| F    | 3.02905564168 | 0.5020800    | -0.19   |
| P    | 3.83086448800 | 2.4476400    | 1.00    |
| S    | 3.56359487256 | 1.8828000    | 0.12    |
| Cl   | 3.31414323148 | 0.9623200    | -0.17   |
| Br   | 3.52795892384 | 1.3388800    | -0.14   |
| I    | 3.99122625727 | 2.1756800    | -0.08   |
| Mg   | 2.11142996199 | 0.0627600    | 2       |
| Zn   | 1.94215920555 | 1.0460000    | 2       |

*Note*: The values of &sigma; and &sigma; were determined based on the charmm36-jul2022.ff force field (atom types: HGA1, HGP1, CG321, NG331, OG2D1, FGR1, PG0, SG311, CLGR1, BRGR1, IGR1, MG, ZN), as provided by the [MacKerell lab website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs).

