

CODE SAMPLES FROM PHD WORK
------------------------------

First committed to git by Max D. Porter on 8/29/19

In this directory are some code samples that show my coding ability. The 
samples are taken from throughout the five years of my physics PhD (2014-2019) at the 
University of Texas at Austin, all written by myself. The code is all Mathematica (Wolfram 
Language), which is difficult to combine with control version systems. As such the files 
have been converted from their original notebook format (.nb) that possessed output, 
figures, and better formatting to a simpler markdown (.md) file that displays well on 
GitHub. Separations between 'cells' (chunks of input) have been lost by this, and some 
formatting such as renaming Greek characters has been necessary. Hopefully the results are 
nonetheless very readable. The original notebook (.nb) files are included for download if 
desired.


FILES:

1. ET_1D-Make_VV1_git.md

This file was for my final project, concerning the 'elliptic-theta one-dimensional 
lattice' model. This project became the paper "Bound states in the Brillouin zone 
continuum" (https://doi.org/10.1016/j.physb.2019.06.057). In this file I construct a 
Hamiltonian matrix corresponding to the potential term in our Hamiltonian, which my group
called 'VV1'. This file outputs the resulting VV1 (which can get very large, up to the 
gigabyte scale) to be read in by another Mathematica notebook and used to construct the 
quantum eigenstates of the system.

2. Sinai-Classical_simulate.md

This comes from my second project, concerning a 'soft Sinai lattice' model, based on the classic Sinai billiard. This project led to two papers: "Chaos in the band structure of a soft Sinai lattice" (https://doi.org/10.1103/PhysRevE.95.052213) and "Signatures of chaos in the Brillouin zone" (https://doi.org/10.1063/1.5001186). In this file that I've heavily abbreviated from the .nb file I simulate a classical point-particle traveling around the lattice. I repeat the simulation for many different initial conditions, and each time harvest the 'surface of section' data from the simulation. This data is then plotted to create a surface of section for analyzing chaos. Care is taken to preserve sufficient precision throughout the simulation based on a specified desired precision. Abbreviated from this .md file are many further analyses of this simulation data.
