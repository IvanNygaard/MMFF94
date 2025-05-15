# MMFF94
An attemnt at a partial implementation of the MMFF94 chemistry force field in C++, developed as a final project for the course "CHEM179: Numerical Algorithms Applied to Computational Quantum Chemistry" at UC Berkeley (Spring 2025). The implementation is partial in the sense that it only supports saturated hydrocarbons. Most of the work is based on the equations presented in the article:
"Merck Molecular Force Field. I. Basis, Form, Scope, Parameterization, and Performance of MMFF94" by T. A. Halgren, J. Comp. Chem., 17, 490–519 (1996).

MMFF parameters were retrieved from:
https://ftp.wiley.com/public/journals-back/jcc/suppmat/17/490.

Molecular geometries for sample inputs were either retrieved from:
https://cccbdb.nist.gov/expgeom1x.asp
or built using Avogadro.

Analytical definitions of angles can be found in Chapter 4 of the book:
https://archive.org/details/E.BrightWilsonJR.J.C.DeciusPaulC.CrossMolecularVibrationsTheTheoryOfInfraredAndRama/page/n59/mode/2up?view=theater. 

Ignore the note "My code: 0.807905 kcal/mol" in the presentation — this result was based on a geometry I built manually, which likely differs from the one Halgren used. I later ran the same calculation (barrier height of gauche–anti n-butane) using RDKit’s MMFF94 implementation and obtained 0.7822 kcal/mol, while my own code (using the same geometry) yielded 0.7852941 kcal/mol. I recall reading that MMFF94 parameters may have been updated at some point, so it's unclear to me whether RDKit uses the exact same parameter set I do or if there's still a pesky bug hiding somewhere in my code.

## Installation & Usage
```bash
# Clone
git clone https://github.com/IvanNygaard/MMFF94.git
cd MMFF94

# Examples: 
./compute_energy ../sample_input/butane/butane_gauche.txt 
./compute_energy ../sample_input/butane/butane_anti.txt
