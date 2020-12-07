"""
-----------------------
Last updated: 06/12/2020 
------------------------

Prerequisite packages reuqired: OpenBabel (and the API python package, pybel, as linked here: https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) 


Summary
-------

This package builds the central sphere of the NP in question, using a voronoi approach. The input is as follows:


"""

import numpy as np
import csv

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


from mpl_toolkits.mplot3d import proj3d 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import scipy as sp
from scipy.spatial import SphericalVoronoi, geometric_slerp

class CentralCoreGenerator:
    """
    NEW 
    """
    def __init__(self, points, R, outputPDB, center):
        self.points = points
        self.R = R
        self.outputPDB = outputPDB
        self.center = center 
    def Sphere:
        """
        Write custom spheres - hypothetical structures
        """
        VoroniBaseSphere = SphericalVoronoi(self.points, self.R, self.center)
        
    def Write:
        with open(self.outputPDB, "w") as f:
            writer = csv.writer(f)
            writer.writerows()
        
      
    
