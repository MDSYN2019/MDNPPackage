[![Build Status](https://app.travis-ci.com/MDSYN2019/MDNPPackage.svg?branch=master)](https://app.travis-ci.com/MDSYN2019/MDNPPackage)

# Martini-PyNP - the constructor for Coarse-grained/All-atomic Molecular Dynamics simulations

Last Updated: 02/09/2022
------------------------

This package unites the many heuristic NP constructor methods and constructs bniomolecular and physical systems with the NPs, specifically with the Martini forcefield. The types
of NPs supported with this package are:

- Ligand-functionalized NPs
- Carbon nanotubes 
- Permutations of the C60 buckyball carbon nanoparticle
- Insertion into systems by leveraging the polyply program - https://github.com/marrink-lab/polyply_1.0 - not yet implemented

At the point of writing this (which is shown in the 'Last Updated' part of this documentation), the Ligand-functionalized 
part of this project is nearly complete, and the next step would be to integrate this with polyply. The way I am trying to add 
this project as a 'plugin' to that project is something I need to consider, and currently the main bottleneck I am facing. 

Carbon nanotubes is part of the roadmap, but will most likely follow after the long work of integrating this work into 
polyply and refactoring this code.  

## Requirements

Martini-PyNP requires:

* Python3
* [NumPy](http://www.numpy.org/)
* [simpletraj](https://github.com/arose/simpletraj)
* [MDAnalysis](https://www.mdanalysis.org/)
* [rdkit](https://www.rdkit.org/)
* [parmed](https://parmed.github.io/ParmEd/)
* [vermouth-martini](https://github.com/marrink-lab/vermouth-martinize)

### Most recent features to be updated:

- The Martini version of the NP builder can now account for building the itp file 
  from ligand to NP. Confirmed that minization works well. 
  
- All-atomic version can construct a NP, but require additional work to make the right itp file  

### Immediate future issues to be resolved:

- - The status as of 02/09/2022 is that I am working on further integration with polyply, but require some work on understanding how
parts of this code can be refactored. Integration with Martini Vermouth and polyply - the documentation on vermouth is scarce, meaning at the moment I have to look into the polyply code to see where exactly I can add this code as a 'plugin'. 
  
  
