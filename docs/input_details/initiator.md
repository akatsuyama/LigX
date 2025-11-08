# Initiator Preparation
An initiator can be created by the manual modification of the xyz file. An initiator must contain the following special element, which can be prepared by replacing the atom names in the original xyz file. It is recommended that a hydrogen atom is replaced to a skeletal connector (SC) or functional group connector (FC). 

## 1. skeletal connector (SC) or functional group connector (FC)

To define a SC, the atom name must be replaced to **VC_2**, and similarly, atom name **VLB** can be used to define a FC.

## 2. Special atom names bonded to SC or FC

A sequence of three connected atoms—*comprising an atom directly bonded to the SC or FC and the two atoms adjacent to it*—must be renamed.<br>
The atoms requiring renaming are schematically highlighted in the following chemical structures.<br>
<pre>
        H   H
        |   |   |
SC/FC – X - Y - C - H
        |   |   |
        H   O
            |
            H
</pre>
The atoms labeled **SC, FC, X and Y** in the above figure should be renamed as specified in the following table.
[Special atom names used in LigX](atomname.md) summarizes the requirements for each atomic name.


### Table for initiator with SC

| SC | X | Y | remarks |
|---|---|---|---|
| VC_2 | XC3_2 | R'Z'_2 | 'Z' = H, C, O, N, S, Cl, Br, I, P, F | 
| VC_2 | XC2_2 | R'Z'_2 | 'Z' = H, C, O, N, S, Cl, Br, I, P, F | 
| VC_2 | XC1_2 | R'Z'_2 | 'Z' = H, C, O, N, S, Cl, Br, I, P, F | 
| VC_2 | XC3g_2 | R'Z'_2 | 'Z' = H, C, O, N, S, Cl, Br, I, P, F | 
| VC_2 | XCA | R'Z'_2 | 'Z' = H, C, O, N, S, Cl, Br, I, P, F | 

*Note* 'Z' used to define Y is a variable. The atom name for Y must be given based on the element.

### Table for initiator with FC

| FC | X | Y | remarks |
|---|---|---|---|
| VLB | YCA | CX_0 or NX_0 | CX_0: FC-C-C, NX_0: FC-C-N | 
| VLB  | YO | CX_0 or NX_0 | CX_0: FC-C-C, NX_0: FC-C-N | 
| VLB  | YN3 | CX_0 | fragment linking begins with ammonium nitrogen | 
| VLB  | YN2 | CX_0 | fragment linking begins with aniline nitrogen | 
| VLB  | YNA | CX_0 | fragment linking begins with other nitrogens | 
| VLB  | YS | CX_0 | fragment linking begins with sulger (thioether) | 
| VLB  | YSA | YN3 | fragment linking begins with sulger (sulfonamide) | 

Other elements attached on FC/SC are not supported.

## 3. Terminal atom (Q atom)

To inform LigX program of the end of initiator, a terminal atom (Q atom) must be defined. This can be achieved by simply added “Q” after the element listed in the last line of the xyz file.
*example*: H → **HQ**

*Note*: If any of the atoms corresponding to SC, X, or Y appear in the last line of the file, please swap its line with that of another atom so that it does not remain on the final line before setting the Q atom.

## 4. Example of xyz files

### a. Phenyl group with SC

```initiator_phenyl_SC.xyz
12
initiator_phenyl_SC
XCA_2 31.430000 2.169000 29.525000
RC_2 31.000999 0.823000 29.065001
C 31.391001 0.264000 27.724001
C 32.221001 1.083000 26.849001
C 32.676998 2.428000 27.245001
C 32.269001 2.945000 28.582001
VC_2 31.106001 2.693000 30.707001
H 30.385000 0.206000 29.702000
H 31.056999 -0.720000 27.431000
H 32.519001 0.708000 25.881001
H 33.290001 3.010000 26.573000
HQ 32.612999 3.932000 28.851999
```

*Structural representation*
<pre>
       HQ       H
        \      /
         C == C
        /      \
VC_2 -XCA_2     C - H
        \\    //
        RC_2- C
        /      \
       H        H
</pre>

### b. Phenyl group with FC

```initiator_phenyl_FC.xyz
12        
initiator_phenyl_FC
YCA 31.430000 2.169000 29.525000
CX_0 31.000999 0.823000 29.065001
C 31.391001 0.264000 27.724001
C 32.221001 1.083000 26.849001
C 32.676998 2.428000 27.245001
C 32.269001 2.945000 28.582001
VLB 31.106001 2.693000 30.707001
H 30.385000 0.206000 29.702000
H 31.056999 -0.720000 27.431000
H 32.519001 0.708000 25.881001
H 33.290001 3.010000 26.573000
HQ 32.612999 3.932000 28.851999
```
*Structural representation*
<pre>
       HQ       H
        \      /
         C == C
        /      \
VLB - YCA      C - H
        \\    //
        CX_0- C
        /      \
       H        H
</pre>


## 5. charge.txt

It is a single-line text file containing only integer values for charge (1, 0, -1, etc.).<br>
*example*
> 1
