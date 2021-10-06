"""
Last Updated: 30/03/2021
------------------------

To create the Nanoparticle, the OpenBabel project is used (https://openbabel.org/docs/dev/Introduction/goals.html) 
to change formats (e.g. sdf etc) 




"""


import math
import numpy as np
import scipy

# Importing OpenBabel Libraries
import openbabel
from openbabel import pybel # Proper way to import pybel 

def Nanoparticle_Base_Fibonacci_Sphere(samples=1): """ Function to
    create even points on a sphere for the base of a Nanoparticle.

    Initially most likely I will have to create a coarse-grained
    version first, then convert to """ points = [] phi = math.pi *
    (3. - math.sqrt(5.))  # golden angle in radians for i in
    range(samples): y = 1 - (i / float(samples - 1)) * 2 # y goes from
    1 to -1 radius = math.sqrt(1 - y * y) # radius at y theta = phi *
    i # golden angle increment x = math.cos(theta) * radius z =
    math.sin(theta) * radius points.append((x, y, z)) return points

def Check_Plausibility(Ligand, Surface): """ Checks for whether there
    are steric clashes between the
    
    """

