"""
-----------------------
Last updated: 20/10/2021 
------------------------

Prerequisite packages reuqired: OpenBabel (and the API python package, pybel, as linked here: https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) 


Summary
-------

This package builds the central sphere of the NP in question, using a voronoi approach. The input is as follows:


Useful Links: 

-https://py-sphere-voronoi.readthedocs.io/en/latest/voronoi_utility.html

-https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/


"""

import numpy as np
import csv
import scipy

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
    
