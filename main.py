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


# --

ligand_3_path = gro_path + "/" + "ANIL.gro" # First ligand 
ligand_4_path = gro_path + "/" + "4NIAN.gro" # Second ligand 
itp_3_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/ANIL_cog.itp'
itp_4_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/4NIAN_cog.itp'

CoreGen_I = CentralCoreGenerator(10, 50, [ligand_1_path, ligand_2_path], ['N1', 'N1'] , ['R3', 'R4'],  itp_1_path, itp_2_path, option = 'Striped')
##return_streamlit_NP(CoreGen_I.return_ordered_coordinates()) # show as streamlit dashboard 
CoreGen_I.generate_coordinates('NP_MIMI_NITL')
CoreGen_I.generate_itp('NP_MIMI_NITL')

CoreGen_II = CentralCoreGenerator(10, 50, [ligand_3_path, ligand_4_path], ['N1', 'N1'] , ['R3', 'O4'],  itp_3_path, itp_4_path, option = 'Striped')
##return_streamlit_NP(CoreGen_II.return_ordered_coordinates()) # show as streamlit dashboard 
CoreGen_II.generate_coordinates('NP_ANIL_4NIAN')
CoreGen_II.generate_itp('NP_ANIL_4NIAN')


