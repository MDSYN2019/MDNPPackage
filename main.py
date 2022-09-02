import pandas as pd
import numpy as np 
from MDNPPackage.core.NP_Core import CentralCoreGenerator
from MDNPPackage.vis.NP_Vis import *

# sample NP parameters 
gro_path = "/home/sang/Desktop/git/Martini3-small-molecules/models/gros" 
ligand_1_path = gro_path + "/" + "1MIMI.gro" # First ligand 
ligand_2_path = gro_path + "/" + "2NITL.gro" # Second ligand 
itp_1_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.itp'
itp_2_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.itp'

CoreGen = CentralCoreGenerator(10, 50, [ligand_1_path, ligand_2_path], ['N1', 'N1'] , ['R3', 'R4'],  itp_1_path, itp_2_path, option = 'Janus')
return_streamlit_NP(CoreGen.return_ordered_coordinates()) # show as streamlit dashboard 
CoreGen.generate_coordinates('NP')
CoreGen.generate_itp('NP')


