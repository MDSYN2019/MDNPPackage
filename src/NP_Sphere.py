"""
-----------------------
Last updated: 06/12/2020 
------------------------

Prerequisite packages reuqired: () 

Summary
-------

This package builds the central sphere of the NP in question, using a voronoi approach. The input is as follows:


"""

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


from mpl_toolkits.mplot3d import proj3d 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import scipy as sp
from scipy.spatial import SphericalVoronoi, geometric_slerp
