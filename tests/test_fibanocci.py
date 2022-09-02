import logging
import sys
import pandas as pd
import numpy as np
import pytest

sys.path.append("..")

# testing the following modules
from MDNPPackage.core.NP_Core import CentralCoreGenerator
from MDNPPackage.utils.NP_Utils import generate_core

test_core_len = 50
gro_path = "/home/sang/Desktop/git/Martini3-small-molecules/models/gros" 
ligand_1_path = gro_path + "/" + "1MIMI.gro" # First ligand 
ligand_2_path = gro_path + "/" + "2NITL.gro" # Second ligand 
itp_1_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.itp'
itp_2_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.itp'

def test_fibanocci_sphere():
    """
    """
    CoreGen = CentralCoreGenerator(10, test_core_len, [ligand_1_path , ligand_2_path], ['R7', 'N1'] , ['O2', 'R4'],  itp_1_path, itp_2_path, option = 'Striped')
    fib_points = CoreGen._fibanocci_sphere()
    assert len(fib_points) == test_core_len 
    
def test_core():
    """
    """
    CoreGen = CentralCoreGenerator(10, test_core_len, [ligand_1_path , ligand_2_path ], ['R7', 'N1'] , ['O2', 'R4'],  itp_1_path, itp_2_path, option = 'Striped')
    assert len(CoreGen.sphere_list) == test_core_len 

def test_something():
    """
    """
    
