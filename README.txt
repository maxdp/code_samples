

RESEARCH SAMPLES FROM PHD WORK
------------------------------

First committed to git by Max D. Porter on 8/29/19

In this directory are some code samples for the purpose of showing my coding ability. The samples are taken from throughout my five years of physics PhD work (2014-2019) at the University of Texas at Austin. The code is all Mathematica (Wolfram Language), which does not generally work well with control version systems. As such the files have been converted from their original notebook format (.nb) that possessed output, figures, and better formatting to a simple Wolfram Language file (.wl) that is easier to display on Github. This also helps focus the attention on my code writing, though it may obfuscate the purpose of some lines where the output (figures, numbers) were only for clarity and bug prevention for myself while coding. Some unclear lines have been removed for display purposes.


FILES:

1. ET 1D Lattice - Make VV1 (for git).wl

This file was for my final project, concerning the "elliptic-theta one-dimensional lattice" model. This project became the paper "Bound states in the Brillouin zone continuum." (https://doi.org/10.1016/j.physb.2019.06.057) In this file I construct a Hamiltonian matrix corresponding to the potential term in our Hamiltonian, which among my group's code we called 'VV1'. This file outputs the resulting VV1 (which can get very large, up to the gigabyte scale) to be read in by another file and used to construct the quantum eigenstates of the system.