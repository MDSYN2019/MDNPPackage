[![Build Status](https://app.travis-ci.com/MDSYN2019/MDNPPackage.svg?branch=master)](https://app.travis-ci.com/MDSYN2019/MDNPPackage)

# Martini-PyNP - the constructor for Coarse-grained/All-atomic Molecular Dynamics simulations

Last Updated: 08/05/2022
------------------------

This package unites the many heuristic NP constructor methods and constructs bniomolecular and physical systems with the NPs, specifically with the Martini forcefield. The types
of NPs supported with this package are:

- Ligand-functionalized NPs
- Carbon nanotubes 
- Permutations of the C60 buckyball carbon nanoparticle
- Insertion into systems by leveraging the polyply program - https://github.com/marrink-lab/polyply_1.0

At the point of writing this (which is shown in the 'Last Updated' part of this documentation), the Ligand-functionalized 
part of this project is nearly complete, and the next step would be to integrate this with polyply. The way I am trying to add 
this project as a 'plugin' to that project is something I need to consider, and currently the main bottleneck I am facing. 

## Requirements
Martini-PyNP requires:

* Python3
* [NumPy](http://www.numpy.org/)
* [simpletraj](https://github.com/arose/simpletraj)
* [MDAnalysis](https://www.mdanalysis.org/)
* [rdkit](https://www.rdkit.org/)
* [parmed](https://parmed.github.io/ParmEd/)
