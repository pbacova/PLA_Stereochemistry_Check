Intra- and intermolecular radial distribution functions for specific types of atoms in PLA. The labelling of the atoms is in agreement with the labelling in ../stereochemistry_check/100LA.itp. 
The ditributions are calculated in radial layers centered in each specific atom, taking into account either the atoms belonging to the same chain (intramolecular) or atoms of neighbouring chains (intermolecular).
Procedure:
1) change the name of the gro file at line 92
2) gcc -lm gr_atom_specific_intra_github.c -o gr_atom_specific_intra_github.out
3) ./gr_atom_specific_intra_github.out to see the flags. Flags -x -y nd -z refer to the maximum box sizes in respected direction reached during the simulation (for NPT simulations). The flag -g defines the width of the grid for the histogram of the distances between atoms (i.e., the binning of the x-axis).
4) run the code with the corresponding flags
5) the output is written in a column-like style, where each column corresponds to the radial distribution function for different atom (see the first commented line of the output)

The procedure is identical for gr_atom_specific_inter_github.c file.     
