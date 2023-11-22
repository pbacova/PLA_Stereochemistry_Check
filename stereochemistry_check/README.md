# PLA_Stereochemistry_Check
C codes which change the stereochemistry of poly(lactic acid) (PLA) monomers.  

The input file contains unwrapped coordinates (x,y,z) in Gromacs format (gro) without velocities. It includes "nchain" linear chains, each with "lchain" atoms. The itp file 100LA.itp illustrates the order of the atoms in the gro file (i.e., 100-mer PLA contains 903 atoms) as also reported in doi:10.1021/ct200251x.
Procedure:
1) change the name of your gro file on line 74 in find_sequence_after_run_github.c
2) gcc -lm find_sequence_after_run_github.c -o find_sequence_after_run_github.out
3) see the flags with ./find_sequence_after_run_github.out
4) run ./find_sequence_after_run_github.out -p numerofchains -k numberofatomsperchain (replace with your own numbers)
5) you get a file called "sequence". By definition, "1" stands for D- and "0" for L- stereoisomer of PLA. The length of the sequence (i.e., the number of lines in the file) should correspond to the number of monomers in the system (each inner monomer of PLA contain 9 atoms, therefore number of lines in sequence file=numerofchains*int(numberofatomsperchain/9)). The labeling starts from monomer 1 in the first chain (line 1) until the last monomer in the last chain (last line).
6) compare this sequence with your desired sequence and create a file called "sequence_wrong" with only 0 or 1. In "sequence_wrong" 0 is used to label monomers, which maintain their stereochemistry and 1 for monomers, whose stereochemistry will be changed by the code.
7) change the name of your gro file on line 78 in correct_stereochemistry_github.c and copy "sequence_wrong" file into the same directory
8) gcc -lm correct_stereochemistry_github.c -o correct_stereochemistry_github.out
9) see the flags with ./correct_stereochemistry_github.out
10) run ./correct_stereochemistry_github.out -p numerofchains -k numberofatomsperchain
11) the output modified2.gro contains the corrected sequence. This sequence can be checked again by repeating steps 1-4. This configuration may contain overlaps, therefore energy minimization needs to be applied before an MD run.
