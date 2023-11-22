Code for calculating the internal distances R_n^2/n of PLA chains.
The input file contains unwrapped coordinates (x,y,z) in Gromacs format (gro) without velocities. 
It includes "nchain" linear chains, each with "lchain" atoms. The itp file ../stereochemistry_check/100LA.itp illustrates the order of the atoms in the gro file (i.e., 100-mer PLA contains 903 atoms) as also reported in doi:10.1021/ct200251x.
The output is a dependence of n (index difference between two backbone atoms) vs. Rn^2/n, where Rn is the distance between the given two atoms i and i+n on the backbone.

Procedure:
1) change the name of the gro file at line 78
2) gcc -lm inner_dist_github.c -o inner_dist_github.out
3) see the flags with ./inner_dist_github.out
4) run ./inner_dist_github.out -f 100 -p 70 -k 903 -o 10 for a system of 70 100-mer chains with 100 frames, from which first 10 will be skipped.
5) use "gnuplot inner_dist_github.g" to plot the results 
