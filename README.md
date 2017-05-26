# eAromatic

Analysis of aromatic interactions in protein-ligand complexes

eAromatic calculates two geometric parameters for aromatic contacts between ligands and proteins. The first parameter is the Cartesian distance between the geometric centers of two rings, referred to as the distance. The second parameter is an angle between normal vectors of two aromatic rings, referred to as the angle.

I. Prerequisite Perl modules:

1. File::Slurp
2. Chemistry::Mol
3. Chemistry::File::MDLMol
4. Chemistry::Ring::Find
5. Math::Vector::Real
6. Math::Trig

II. Input files:

1. Protein in PDB format (mandatory)
2. Ligand in MOL format (mandatory)
3. Output from LPC (optional)

III. Output:

A list of interacting pairs of aromatic rings, each annotated with 
the distance and the angle.

IV. Note on file formats:

We recommend OpenBabel (http://openbabel.org) to convert between 
different file formats.

V. Note on contacts:

You can define contacts either by providing the output from Ligand 
Protein Contacts (LPC) (https://www.ncbi.nlm.nih.gov/pubmed/10320401) 
or with an atomic distance threshold.

VI. Example:

1. Run eAromatic without arguments to get the list of available options:

perl eAromatic-1.0.pl 

2. Run eAromatic with distance-based contacts at a default 4.5 A 
threashold:

perl eAromatic-1.0.pl -p 3lc6A.pdb -l 3lc6A00.mol

3. Run eAromatic with distance-based contacts at a different 
threshold, e.g. 8.0 A:

perl eAromatic-1.0.pl -p 3lc6A.pdb -l 3lc6A00.mol -t 8.0

4. Run eAromatic with LPC-based contacts:

perl eAromatic-1.0.pl -p 3lc6A.pdb -l 3lc6A00.mol -c 3lc6A00.lpc

5. Output should look like this:

AROMATIC TYR  467  1  4.0371    6.32  6 25:23:15:13:16:24

AROMATIC TYR  467  0  3.9689    7.38  5 25:17:14:22:24

2nd column is the residue name

3rd column is the residue index

4th column is the index of ligand aromatic ring (these indices start 
from 0)

5th column is the distance (as defined above)

6th column is the angle (as defined above)

7th column is the number of atoms in ligand ring

8th column is the list of atom indices for that ligand ring

For tryptophan residues, TR5 and TR6 denote 5- and 6-member rings, 
respectively.

The first detected aromatic contact is between aromatic ring of 
tyrosine 467 and a 6-member ring in the ligand comprising atoms 13, 
15, 16, 23, 24, and 25. The Cartesian distance between ring geometric 
centers is 4.0371 Angstrom, whereas the angle between normal vectors 
of these rings is 6.32 deg.
