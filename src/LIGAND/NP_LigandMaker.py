import numpy as np
import scipy 
import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda

"""
------------------------
Last Updated: 22/10/2021 
------------------------

Breakdown: 
----------

The class MolecularConverter translates the smiles string of the ligand(s)
into a cartesian coordinate xyz file. From this, we use the indices of the
smiles sting and match it to the SmilesToMartiniDictionary. From this, we
can use the center of geometry of the indices and replace with the relevant 
Martini bead.   


"""

class MolecularConverter:
    """
    Written Manual here 
    """    
    def __init__(self, option, smilesString):
        """
        The constructor must detect whether the smiles string is valid or not  
        """

        


        # Adding Martini Beads dictionary - at the moment, I am using
        # Martini 3 references
        
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

    def ComputeCoordinatesLigands(self):
        """
        Get the generic basis of the xyz coorinates of the 
        ligands 
        """
        u = mda.Universe.from_smiles(smilesString)
        # new feature
        Molecule = u1.select_atoms('all')
        MoleculeAtomPositions = Molecule.positions # Finds 
        
        
