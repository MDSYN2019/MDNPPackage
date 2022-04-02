[![Build Status](https://app.travis-ci.com/MDSYN2019/MDNPPackage.svg?branch=master)](https://app.travis-ci.com/MDSYN2019/MDNPPackage)

# Martini-PyNP - the constructor for Coarse-grained/All-atomic Molecular Dynamics simulations

Last Updated: 02/04/2022
------------------------


This package unites the many heuristic NP constructor methods and constructs bniomolecular and physical systems with the NPs, specifically with the Martini forcefield. The types
of NPs supported with this package are:

- Ligand-functionalized NPs

- Carbon nanotubes 

- Permutations of the C60 buckyball carbon nanoparticle

- Insertion into systems by leveraging the polyply program - https://github.com/marrink-lab/polyply_1.0


## Requirements
Martini-PyNP requires:
* Python3
* [NumPy](http://www.numpy.org/)
* [simpletraj](https://github.com/arose/simpletraj)
* [MDAnalysis](https://www.mdanalysis.org/)
* [rdkit](https://www.rdkit.org/)
