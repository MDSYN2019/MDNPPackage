# MDAnalysis and Alignment libraries
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np 
import pandas as pd 

# Scipy libaries                                                                                                                             from scipy.sparse import csr_matrix                                                                                                          from scipy.sparse.csgraph import floyd_warshall                                                                                              from scipy.spatial import ConvexHull, convex_hull_plot_2d  

# Boilerplate libraries                                                                                                                      import sys                                                                                                                                   import re                                                                                                                                                                                                 
import math                                                                                                                                  import scipy

# scipy libaries                                                                                                                                                                                           
from scipy.sparse import csr_matrix                                                                                                                                                                       
from scipy.sparse.csgraph import floyd_warshall                                                                                                                                                           
from scipy.spatial import ConvexHull, convex_hull_plot_2d  

import numpy as np
import scipy 
import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda


"""
------------------------
Last Updated: 02/07/2022
------------------------
"""

class MolecularConverter:
    """
    """    
    def __init__(self, gro, firstatoms, lastatoms, spherelist, option = 'plain'):
        """
        """
        self.gro = gro
        self.firstatoms = firstatoms
        self.lastatoms = lastatoms
        self.spherelist = spherelist
        self.option = option
    
