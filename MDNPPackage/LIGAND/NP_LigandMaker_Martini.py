# MDAnalysis and Alignment libraries
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np 
import pandas as pd 

# Scipy libaries                                                                                                                             from scipy.sparse import csr_matrix                                                                                                          from scipy.sparse.csgraph import floyd_warshall                                                                                              from scipy.spatial import ConvexHull, convex_hull_plot_2d  

# Boilerplate libraries                                                                                                                      import sys                                                                                                                                   import re                                                                                                                                                                                                 
import math                                                                                                                                  import scipy

# scipy libaries                                                                                                                                                                                           
from scipy.sparse import csr_matrix                                                                                                                                                                       
from scipy.sparse.csgraph import floyd_warshall                                                                                                                                                           
from scipy.spatial import ConvexHull, convex_hull_plot_2d  

import numpy as np
import scipy 
import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda


"""
------------------------
Last Updated: 29/05/2022
------------------------
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
        
        #self._SmilesToMartiniDictionary = {}
        #SmilesToMartiniDictionary["CC(=O)CO"] = 'P2' # P2 Bead 
        #SmilesToMartiniDictionary["CC(=O)O"] = 'SP2' # SP2 Bead 
        #SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 
        #SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 
        
    def ReadMartiniMolecules(GroFile, First, Last):
        """ Generate the normalized coordinates, name, and vector of the Martini molecule 
    
        Access the Martini3 small molecules library and reads the parameterized coordinates from it, 
        with future view of looking at generating automatically generating Martini 3 representations 
        from smiles strings 
    
        One needs to describe the attaching bead to the main core and the atom furthest away from the 
        core, to create the directional vector to which the struture will be placed on the surface of the NP 
        core. 
    
        Args:
            GroFile:
              path the gromacs file of the ligand
        Returns: 
            Placeholder
        Raises: 
            Placeholder 
        
        """
        TransformationList= []
        MartiniUniverse = mda.Universe(GroFile) # Load the Martini gro file in as a universe 
        ids = [i.name for i in MartiniUniverse.atoms]
    
        Molecule = MartiniUniverse.select_atoms('all')
        # In this case, the atoms will be N1 and R3 
        FirstAtom = Molecule.select_atoms('name {}'.format(First)) 
        LastAtom = Molecule.select_atoms('name {}'.format(Last))
        LigandAlignmentVector = (FirstAtom.positions - LastAtom.positions)[0] # Get the alignment vector created from the first and COM  
    
        # Loop over the positions 
        for i,j in enumerate(Molecule.positions):
            vector = (j - FirstAtom.positions)[0]
            vector[0] = LigandAlignmentVector[0] - vector[0]
            vector[1] = LigandAlignmentVector[1] - vector[1]    
            vector[2] = LigandAlignmentVector[2] - vector[2]
            if vector[0] == -math.inf:
                pass
            if vector[0] == 0.0:
                pass
            else:
                TransformationList.append([vector, Molecule.atoms[i].type])   
            
        # Return the universe, the transformed (normalized) coordinate list of the ligand molecule, and the 
        # alignment vector that shows the arrow of direction of the vector, which we will be able to reorientate
        return Molecule, TransformationList, LigandAlignmentVector

        def AttachLigandsMartini(GroFiles, FirstAtoms, LastAtoms, SphereList, option = 'Plain'):
        """ Placeholder 

        Here, we follow the same logic as the the AttachLigands to create a Martini version 
        of it. We are currently only using the Martini3 small molecules dictionary to create 
        the martini ligands 
    
        Args: 
            Placeholder
        Returns: 
            Placeholder
        Raises:
            Placeholder
        """
        if option == 'Plain':
            # Information here 
            Molecule, TransformationList, LigandAlignmentVector = ReadMartiniMolecules(
                GroFiles[0], FirstAtoms[0], LastAtoms[0])
        
            MolLen =  len(Molecule)
            Coordinates = ReturnPandasNPMartini(Molecule, LigandAlignmentVector, 
                                              TransformationList, SphereList[0])
            return Coordinates
    
        if option == 'Janus' or option == 'Striped':
            # Information here 
        
            # First ligand 
            Molecule_I, TransformationList_I, LigandAlignmentVector_I = ReadMartiniMolecules(
                GroFiles[0], FirstAtoms[0], LastAtoms[0])
        
            # Get length of first ligand
            MolLen_I = len(Molecule_I)
            # Second Ligand 
            Molecule_II, TransformationList_II, LigandAlignmentVector_II = ReadMartiniMolecules(
                GroFiles[1], FirstAtoms[1], LastAtoms[1])
        
            # Get length of second ligand
            MolLen_II = len(Molecule_II)
        
            Coordinates_I = ReturnPandasNPMartini(Molecule_I, LigandAlignmentVector_I, 
                                                  TransformationList_I, SphereList[0], 'Lig1', 'Core').reset_index()
            # Add indices for Coordinates_I 
            Coordinates_I['index'] = Coordinates_I.index
            Coordinates_II = ReturnPandasNPMartini(Molecule_II, LigandAlignmentVector_II, 
                                              TransformationList_II, SphereList[1], 'Lig2', 'Core').reset_index()
            # Add indices for Coordaintes_II 
            Coordinates_II['index'] = Coordinates_II.index
            # This bit will be migrated to later on 
            Coordinates = Coordinates_I.append(Coordinates_II)
            logging.info(f"{Coordinates}")
        
            # We return the coordinates of the first batch of NP atoms attached with ligands, 
            # and then attach the second batch of NP atoms attached with the second type of ligands 
            # Add in the indices of the Coordinates first
        
            SkippedLigand_I_from_SphereList = Coordinates_I[:-len(SphereList[0])].iloc[::MolLen_I, :]
            SkippedLigand_II_from_SphereList = Coordinates_II[:-len(SphereList[1])].iloc[::MolLen_II, :]
        
            #logging.info(f"{SkippedLigand_I_from_SphereList}, {SkippedLigand_II_from_SphereList}")
        
            # Assert that the filtered coordinates of the ligands is identical to the 
            # NP atoms for the 
            assert len(SkippedLigand_I_from_SphereList) == len(SphereList[0])
            assert len(SkippedLigand_II_from_SphereList) == len(SphereList[1])
        
            AttachementBonds = []    
            Dist_LigI = [[SkippedLigand_I_from_SphereList['index'].iloc[i] + 1, 
                          Coordinates_I[-len(SphereList[0]):]['index'].iloc[i] + 1,
                          distance.euclidean(SkippedLigand_I_from_SphereList.iloc[i][['X', 'Y', 'Z']].to_numpy(),
                                             SphereList[0][i])] for i in range(0, len(SkippedLigand_I_from_SphereList))]
        
            Dist_LigII = [[SkippedLigand_II_from_SphereList['index'].iloc[i] + len(Coordinates_I) + 1, 
                           Coordinates_II[-len(SphereList[1]):]['index'].iloc[i] + len(Coordinates_I) + 1,
                          distance.euclidean(SkippedLigand_II_from_SphereList.iloc[i][['X', 'Y', 'Z']].to_numpy(), 
                                             SphereList[1][i])] for i in range(0, len(SkippedLigand_II_from_SphereList))]
        
            AttachementBonds.append('; Ligand - NP bond')
        
            for bond in Dist_LigI:
                bondstring = f"{bond[0]} {bond[1]} 1 {bond[2]} 10000"
                AttachementBonds.append(bondstring)
            for bond in Dist_LigII:
                bondstring = f"{bond[0]} {bond[1]} 1 {bond[2]} 10000"
                AttachementBonds.append(bondstring)
            
            logging.info(f"{Dist_LigI}, {Dist_LigII}")
            Coordinates = Coordinates.reset_index()
            # Now add into a list with indices, and return 
            return Coordinates, AttachementBonds
