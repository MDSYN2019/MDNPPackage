import logging
import sys
sys.path.append("..")

import pytest
import pandas as pd
import numpy as np

# Testing the following modules
from MDNPPackage.core.NP_Core import CentralCoreGenerator
from MDNPPackage.utils.NP_UTILS import generate_core

testcorelen = 50
gropath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/gros" 
MIMILigandPath = gropath + "/" + "CAFF.gro" # First ligand 
NITLLigandPath = gropath + "/" + "2NITL.gro" # Second ligand 
itppath1 = '/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.top'
itppath2 = '/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.top'

def test_fibanocci_sphere():
    """
    """
    CoreGen = CentralCoreGenerator(10, testcorelen, [MIMILigandPath, NITLLigandPath], ['R7', 'N1'] , ['O2', 'R4'],  itppath1, itppath2, option = 'Striped')
    fibpoints = CoreGen._fibanocci_sphere()
    assert len(fibpoints) == testcorelen 
    
def test_core():
    """
    """
    CoreGen = CentralCoreGenerator(10, testcorelen, [MIMILigandPath, NITLLigandPath], ['R7', 'N1'] , ['O2', 'R4'],  itppath1, itppath2, option = 'Striped')
    assert len(CoreGen.spherelist) == testcorelen 
