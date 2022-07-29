import pandas as pd
import numpy as np 
from MDNPPackage.core.NP_Core import CentralCoreGenerator

gropath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/gros" 
MIMILigandPath = gropath + "/" + "CAFF.gro" # First ligand 
NITLLigandPath = gropath + "/" + "2NITL.gro" # Second ligand 

itppath1 = '/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.top'
itppath2 = '/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.top'

CoreGen = CentralCoreGenerator(10, 50, [MIMILigandPath, NITLLigandPath], ['R7', 'N1'] , ['O2', 'R4'],  itppath1, itppath2, option = 'Striped')

CoreGen.generate_coordinates('NP')
CoreGen.generate_itp('NP')


