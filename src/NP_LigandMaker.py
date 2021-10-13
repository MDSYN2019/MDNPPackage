import numpy as np
import scipy 
import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda

"""
Last Updated: 17/12/2020 

Author: Sang Young Noh 

TODO 

"""
class Converter():
    """
    Written Manual here 
    """
    def __init__(self, option):
        """
        What does this constructor do? 
        """
        if option not in ['mol', 'sdf']:
            except ValueError:
                print("Need a mol or sdf input") 
        else:
            self.option = option
            
    def GetRingSystems(mol, includeSpiro=False):
        """
        What is this function doing?
        """
        ri = mol.GetRingInfo() # Sets out the indices of the rings structures within the mol file 
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
                nSystems.append(ringAts)
                systems = nSystems
        return systems
