
CODE SAMPLES FROM PHD WORK
------------------------------
NOTE: Look at .md files, not .nb!

First committed to git by Max D. Porter on 8/29/19

Last updated on 9/19/19

In this directory are some code samples that show my coding ability. The samples are taken 
from throughout the five years of my physics PhD (2014-2019) at the University of Texas at 
Austin. All code is written 100% on my own. The code is all Mathematica (Wolfram 
Language), as I don't have code in other languages from previous employers which I can
share. I would be happy to demonstrate my coding ability in another language, and can
learn new languages very quickly.

Mathematica is difficult to combine with control version systems. As such the files 
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
Hamiltonian matrix corresponding to the potential term in the Schrodinger equation, which 
my group called 'VV1'. This file outputs the resulting VV1 (which can get very large, up 
to the gigabyte scale) to be read in by another Mathematica notebook and used to construct 
the quantum eigenstates of the system. Constructing these states was the starting point for Wigner-Eisenbud scattering analysis, for which I also did much of the coding.

2. Sinai-Classical_simulate.md

This comes from my second project, concerning a 'soft Sinai lattice' model, based on the classic Sinai billiard. This project led to two papers: "Chaos in the band structure of a soft Sinai lattice" (https://doi.org/10.1103/PhysRevE.95.052213) and "Signatures of chaos in the Brillouin zone" (https://doi.org/10.1063/1.5001186). In this file I simulate a classical point-particle traveling around two versions of our model lattice: a symmetric case formed from one periodic Gaussian (one bump) and an asymmetric case formed from three overlapping periodic Gaussians (three bumps). The code also allows for other cases we considered, like inverting the bumps to make holes. I repeat the simulation for many different initial conditions, and each time harvest the 'surface of section' data from the simulation. This data is then plotted to create a surface of section for analyzing chaos. Care is taken to preserve sufficient precision throughout the simulation based on a specified desired precision. I've heavily abbreviated this file from the .nb file, which contains further analysis of this simulation data.
