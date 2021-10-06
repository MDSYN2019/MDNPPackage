"""

NP Construction Library 
-----------------------

Last Updated: 30/9/2021
-----------------------

# These need to be implemented as part of the smiles functionality where we can add identify what can be martinized 
# and added onto the surface of the NPs. 

-> Useful links:

- https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html - help on the numpy functionality 


-> https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/ - smiles functionaliry

"""
import MDAnalysis as mda
from rdkit import Chem
import math
import numpy as np
import scipy
import argparse

def Nanoparticle_Base_Fibonacci_Sphere(samples=1):
    """ 
    Function to create even points on a sphere for the base 
    of a Nanoparticle - in the form of a fibanocci sphere 

    Initially most likely I will have to create a coarse-grained
    version first, then convert to 


    Format follows the numpy convention: 
    
    https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html 

    Parameters 
    ----------

    Returns 
    -------
    """
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2 # y goes from1 to -1
        radius = math.sqrt(1 - y * y) # radius at y
        theta = phi * i # golden angle increment
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append((x, y, z))
    return points

def Check_Plausibility(Ligand, Surface):
    """ 
    Checks for whether there are steric clashes between the beads within the made NP. 
    """
    pass

def Create_Smiles_Ligand(SmileString):
    """
    Write tutorial here 
    """
    u1 = mda.Universe.from_smiles(SmileString)
    mol1 = u1.atoms.convert_to("RDKIT")
    u2 = mda.Universe(mol1)
    u2.atoms.write("mol.pdb")
    mol2 = Chem.MolFromPDBFile("mol.pdb")
    

parser = argparse.ArgumentParser()
#parser.parse_args()
parser.add_argument("square", help = "display the square of a given number", type = int)
parser.add_argument("coreNumber", help = "Number of beads to constitute the core of a ligand-attached NP", type = int)
# Optional arguements - we have positional arguments
parser.add_argument("--verbosity", help = "increase output verbosity")

parser.add_argument( "-nr", "--numrings",   type=int,   default=12,     help='Number of rings (~height)'                           )

args = parser.parse_args()

if args.verbosity:
    print("verbosity turned on")
    
#print(args.square)
print(args.square**2)
print(type(args.square))

if condition:
    print( "------------------------------------------------------------------------------" )
    print( "Generating a Martini model for an open CNT using "+str(numrings)+" rings with "+str(ringsize)+" each." )
    print( "------------------------------------------------------------------------------" )

else:
    pass



