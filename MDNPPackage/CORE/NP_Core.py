"""
-----------------------
Last updated: 13/03/2021 
------------------------

Prerequisite packages reuqired: OpenBabel (and the API python package, pybel, as linked here: https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) 


Summary
-------

This package builds the central sphere of the NP in question, using a voronoi approach. The input is as follows:

Useful Links: 
-------------



-> https://py-sphere-voronoi.readthedocs.io/en/latest/voronoi_utility.html

-> https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/

-> https://cedric.bouysset.net/blog/2020/08/07/rdkit-interoperability

-> https://stackoverflow.com/questions/47319238/python-plot-3d-vectors
 
-> https://en.wikipedia.org/wiki/XYZ_file_format - information on the xyz coordination file 

-> https://mattermodeling.stackexchange.com/questions/3961/recalculate-atom-positions-to-account-for-periodic-boundary-conditions/3970#3970  - Atomsk

-> https://stackoverflow.com/questions/49064611/how-to-find-different-groups-in-networkx-using-python - grouping networks 

-> https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/

-> http://cgmartini.nl/index.php/component/kunena/8-martini-philosophy/5776-mapping-of-benzene-ring

-> https://docs.mdanalysis.org/1.0.0/documentation_pages/lib/NeighborSearch.html

-> https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

-> https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html

-> https://nanotube.msu.edu/fullerene/fullerene-isomers.html

-> https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753

-> https://plotly.com/python/3d-scatter-plots/


Quote by Riccardo: 
------------------

In the open beta of Martini 3, benzene is indeed mapped with three TC4 beads, that that's again a 2-to-1 mapping. 
The bond length is changed to 0.29 nm in 3.0 because this allows to represent more closely the volume of a benzene 
molecule, taking into account also the smaller size of T-beads as compared to S -beads

---------------------------
How the modules are divided
---------------------------

-> We want to be given free reign into making NPs of multiple types. The primary types we are concerned with at the 
   moment is 

   1. Functionalized AuNP type ones.  
   2. Carbon Nanotube like ones. 
   3. Large spherical buckyball type structures - Like C70. 

-> TODO 
    
   - Am not able to get the box dimenisons correct for the xyz file...
   - Need to write Striped/Janus functionality to the NP maker  


Trying to translate the atomic C70 structure into the coarse-grained c70 structure 

TODO - isolate per 2 beads - There is a two to one mapping, and find the center of mass 
       for each of these beads. Then construct the itp file for the bead. 


"""

# Boilerplate libraries                                                                                                                                                                                    
import sys                                                                                                                                                                                                
import re                                                                                                                                                                                                 
import math 

# Rdkit libraries 
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem                                                                                                                                                                            
from rdkit.Chem import ChemicalFeatures                                                                                                                                                                   
from rdkit.Chem import rdchem                                                                                                                                                                             
from rdkit.Chem import rdMolDescriptors                                                                                                                                                                   
from rdkit import RDConfig  

# Alignment libraries in MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, PDB, XTC

# Scipy libraries
import scipy                                                                                                                                                                                                                                                                                                                                                                                                          
from scipy.sparse import csr_matrix                                                                                                                                                                        
from scipy.sparse.csgraph import floyd_warshall                                                                                                                                                            
from scipy.spatial import ConvexHull, convex_hull_plot_2d 
from scipy.linalg import solve
from scipy.spatial import distance

# Matplotlib libraries 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d    
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objs as go

# Pandas 
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import math 
from operator import itemgetter
import itertools                                                                                                                                                                                           
import requests                                                                                                                                                                                           
import collections                                                                                                                                                                                        
import random



#import numpy as np
#import csv
#import scipy

#import matplotlib
#import matplotlib.pyplot as plt
#import matplotlib.colors as colors

#from mpl_toolkits.mplot3d import proj3d 
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#import scipy as sp
#from scipy.spatial import SphericalVoronoi, geometric_slerp
#from MDAnalysis.analysis.distances import distance_array


class CentralCoreGenerator:
    """
    Generation of the structure of the center core of the NP. 
    Initially, we will be using the f
    """
    def __init__(self, filename, points, R, outputPDB, center):
        self.points = points
        self.R = R
        self.outputPDB = outputPDB
        self.center = center 
        self.output =  open(filename+".itp", 'w')
      
    def Nanoparticle_Base_Fibonacci_Sphere(samples=1):
        """ 
        Function to create even points on a sphere for the base of a Nanoparticle.
        """
        points = []
        phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2 # y goes from 1 to -1 radius = math.sqrt(1 - y * y) # radius at y
            theta = phi * i # golden angle increment
            x = math.cos(theta) * radius
            z = math.sin(theta) * radius points.append((x, y, z))
        # Return the surface points on the sphere  
        self.points = points

    def Nanoparticle_Base_Nanotube():
        """
        Function to create even points on a tube for the base of a Nanoparticle.
        
        This function would be based on Martin Vogele's code referenced. 
        """
        pass
        
    def Connect_Ligand():
        """
        Find the connecting vector to each of the core atoms on the surface of the core, and 
        find the closest atom to the anchor atom of the ligand. From this, we can 'attach' the 
        relevant ligand. 
        """
        x = [p[0] for p in self.points]
        y = [p[1] for p in self.points]
        z = [p[2] for p in self.points]
        centroid = (sum(x) / len(self.points), sum(y) / len(self.points), sum(z) / len(self.points)) # Find the center
        # Now need to compute the vectors that connect the center of the core
        # to the surface of the core
        VectorStore = [] 
        for point in self.points:
            pass
        
    def Write_Coordinates:
        """
        This writes the coordinates generated from either the nanotube or 
        the spherical core and converts the xyz file into the pdb file, or leaves the 
        xyz file. 
        """
        with open(self.outputPDB, "w") as f:
            writer = csv.writer(f)
            writer.writerows()
        pass
    
class CoarseGrainer:
    """
    This class creates a MARTINI3 compatible mapping over the 
    smiles pdb of the ligands that has been constructed from the RDKit code. 

    """
    def __init__(self):
        pass
    
