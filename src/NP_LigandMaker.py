import numpy as np
import scipy 
import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda

"""
------------------------
Last Updated: 20/10/2021 
------------------------

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

        # Adding Martini beads here 
        self._SmilesToMartiniDictionary = {}
        SmilesToMartiniDictionary["CC(=O)CO"] = 'P2' # P2 Bead 
        SmilesToMartiniDictionary["CC(=O)O"] = 'SP2' # SP2 Bead 
        SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 
        SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 
    
    def GetRingSystems(self, mol, includeSpiro=False):
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

    def FunctionPlaceHolder(self):
        pass
