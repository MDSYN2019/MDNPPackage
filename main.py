"""
Last Updated: 29/05/2022
------------------------
Author: Sang Young Noh 
----------------------
"""
import sys                                                                                                                                   import re                                                                                                                                    import pandas as pd
import numpy as np
import plotly.graph_objs as go
import math 
from operator import itemgetter
import itertools                                                                                                                             import requests                                                                                                                              import collections                                                                                                                           import random            

# Scipy libraries 
import scipy                                                                                                                                 from scipy.sparse import csr_matrix                           
from scipy.sparse.csgraph import floyd_warshall               
from scipy.spatial import ConvexHull, convex_hull_plot_2d 
from scipy.linalg import solve
from scipy.spatial import distance

# Rdkit libraries 
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem                               
from rdkit.Chem import ChemicalFeatures                       
from rdkit.Chem import rdchem                                 
from rdkit.Chem import rdMolDescriptors                       
from rdkit import RDConfig  

# Alignment libraries in MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
#from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, PDB, XTC

# Matplotlib libraries 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d    
from mpl_toolkits.mplot3d import Axes3D

# plotly functionalities 
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

# Abstract classes and methodology 
from abc import ABC, abstractmethod

# Importing parmed to be able to read the appropriate Martini ligand itp files
# for the ligands

import parmed as pmd
from parmed.gromacs.gromacstop import GromacsTopologyFile

# logging module 
import logging 
logging.basicConfig(level = logging.INFO)



class MDNPPackage(object):
    def __init__(self):
        """
        """
        self.FibSphere = []
        self.BatchSphereCoordinates = None
        
    def Fibanocci_Sphere(samples=1):
        """ Return a Fibanocci sphere with N number of points on the surface. 

        This will act as the template for the nanoparticle core. 
    
        Args:
            Placeholder
        Returns:
            Placeholder
        Raises:
            Placeholder
        """
        #points = []
        phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
            radius = math.sqrt(1 - y * y)  # radius at y
            theta = phi * i  # golden angle increment
            x = math.cos(theta) * radius
            z = math.sin(theta) * radius
            self.FibSphere.append((x, y, z))
            
    def rotation_matrix_from_vectors(vec1, vec2):
        """ Find the rotation matrix that aligns vec1 to vec2
        Args:
            vec1: 
                A 3d "source" vector
            vec2: 
                A 3d "destination" vector
        Returns:
            rotation_matrix:
                A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
        Raises:
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        self.rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

        
    def LabelStripedNP(Core, Type = 'Janus'):
        """Placeholder
    
        Depending on the type of NP we want in the input, we can try to generate 
        different patterns on the surface of the spehre, which will help us generate the 
        lists correponding to the anisotropic nature of the NP. 
    
        The types of NPs we can currently have are:
        
        - Janus 
        - Striped
    
        More options will be added. As in its current iteraton:
    
        1. The Janus type divides the NP into two hemispheres.
    
        2. The Striped type divides the NP into three hemispheres, typically used with a hydrophobic middle 
        especially when it comes to using with biosimulations. 
        
        Args:
            Core: 
                Placeholder
            Type:
                Placeholder
        Returns:
    
        Raises:

        """
        XCoordinates = [i[0] for i in Core] # Find x coordinates
        YCoordinates = [i[1] for i in Core] # Find y coordinates
        ZCoordinates = [i[2] for i in Core] # Find z coordinates 
        Length = 2 * abs(max(ZCoordinates)) # From 2 * the radius, we know the total length of the NP 
    
        if Type == 'Striped':
            # As we have a spherical structure, we just need to find the minimum/maximum in 
            # one of the axes to find that for the rest 
            # define the threshold for how you wish to generate the NP with striped pattern 
            Threshold = Length / 3 
            # Find the central band of the sphere where you wish to put 
            # different ligands 
            StripedValues = [i for i in Core if i[2] > (min(ZCoordinates) + Threshold) and i[2] < (max(ZCoordinates) - Threshold)] # Return middle strip values
            CeilingValues = [i for i in Core if i not in StripedValues] # Return top 'ceiling' strip values 
            return StripedValues, CeilingValues
        elif Type == 'Janus':
            Threshold = Length / 2 
            TopValues = [i for i in Core if i[2] > (min(ZCoordinates) + Threshold)] # Return top hemipshere 
            BotValues = [i for i in Core if i not in TopValues] # Return bottom hemisphere 
            return TopValues, BotValues
        
    def GenerateCore(Radius, N, Option = 'Plain'):
        """ Creates a Fibanocci sphere that represents the NP core 
        and allocates the radius. Using the radius, 
    
        The core is scaled down/up to the size that one wishes to have. 
        We can generate arrays corresponding  to a plain core, or a tuple with 
        two entries with different parts of the NP core that corresponds to positions 
        with striped or janus type positions.
    
        Args:
            Radius: 
                The radius of the core of the NP 
            N: 
                The number of core 'atoms' we want create that constitutes the core 
            Option (default = 'Plain'):
                If the option is plain, we simply return a single list. Otherwise, 
                we generate a tuple which separates out core atoms.
        Returns:
    
        Raises: 
    
        """
        Sphere = FibanocciSphere(N) # Create the fibanocci sphere representing the NP core 
        XSphere, YSphere, ZSphere  = [], [], []
    
        for entry in Sphere:
            XSphere.append(entry[0])
            YSphere.append(entry[1])
            ZSphere.append(entry[2])
        
        # Append as 2d list
        SphereList = [] 
        for index in range(0, len(XSphere)):
            SphereList.append([XSphere[index], YSphere[index], ZSphere[index]])
        # Take the radius value, and then multiply the unit vector in each 
        # Direction by that radius value to increase the total volume of the 
        # NP core.
        for index in range(0, len(SphereList) -1):
            SphereList[index][0] = SphereList[index][0] * Radius
            SphereList[index][1] = SphereList[index][1] * Radius
            SphereList[index][2] = SphereList[index][2] * Radius
        # Return just the whole list without any further modifications
        if Option == 'Plain':
            return [SphereList[1:-1]]
        # Separate out the anisotropy for the striped variant 
        elif Option == 'Striped':
            StripedValues, CeilingValues = LabelStripedNP(SphereList[1:-1], Option)
            return StripedValues, CeilingValues
        # Separate out the anisotropy for the Janus variant 
        elif Option == 'Janus':
            TopValues, BottomValues = LabelStripedNP(SphereList[1:-1], Option)
            return TopValues, BottomValues  
        
        
    def ReturnPandasNP(LigandString, FirstAtom, LastAtom, SphereList, 
                   LigandName, CoreName, Length = 1.0):
        """Placeholder
    
        Also needs to return the indices i,j of the sulfur 
        ligands, with the 
    
        , and the distance between the 
    
        Args:
        Returns:
        Raises:    
        """
        TransformationList, NameList = [], [] # 
        LigandList = [] 
        Sphere = []
        Xplot, Yplot, Zplot = [], [], []
        XplotSphere, YplotSphere, ZplotSphere = [], [], []
        
        u = mda.Universe.from_smiles(LigandString)
        Ligand = u.select_atoms('all')
        logging.info(f"The length of the ligand is {len(Ligand)}")
        FirstAtomGroup = u.select_atoms('name {}'.format(FirstAtom))
        LastAtomGroup = u.select_atoms('name {}'.format(LastAtom))
        LigandAlignmentVector = (FirstAtomGroup.positions- LastAtomGroup.positions)[0]
    
        for i,j in enumerate(Ligand.positions):
            vector = (j - FirstAtomGroup.positions)[0]
            vector[0] = LigandAlignmentVector[0] - vector[0]
            vector[1] = LigandAlignmentVector[1] - vector[1]    
            vector[2] = LigandAlignmentVector[2] - vector[2]
            if vector[0] == -math.inf:
                pass
            if vector[0] == 0.0:
                pass
            else:
                TransformationList.append([vector, Ligand.atoms[i].type])        
        
        unitVector = np.linalg.norm(LigandAlignmentVector)
        vecLigand = LigandAlignmentVector.tolist()
    
        # Loop over the sphere and find the 
        for index in range(0, len(SphereList)):
            vec2 = SphereList[index]
            # Find the transformationvector for the ligand vector to vec2, which is the position of the point on sphere
            TransformationVector = rotation_matrix_from_vectors(vecLigand, vec2)
            # Rotate the vector 
            vec1_rot = TransformationVector.dot(vecLigand) # Rotate the vector to match the surface point on the sphere 
            # Get the absolute length of the unit vector 
            unitVectorabs = np.linalg.norm(LigandAlignmentVector)
            # Change the rotation vector in unit vector, then multiply by the absolute 
            # length of the sphere 
            vecMultiplier = vec1_rot/unitVectorabs * (np.linalg.norm(np.array(vec2))) + (vec1_rot/unitVectorabs * Length)
            # Find the difference in length 
        
            Sphere.append(vec2)
        
            # Translate the vector further out 
            for trans in TransformationList:
                LigandAtomcoordinate = TransformationVector.dot(trans[0])
                LigandAtomcoordinate[0] = LigandAtomcoordinate[0] + vecMultiplier[0]
                LigandAtomcoordinate[1] = LigandAtomcoordinate[1] + vecMultiplier[1]
                LigandAtomcoordinate[2] = LigandAtomcoordinate[2] + vecMultiplier[2]
                LigandList.append(LigandAtomcoordinate.tolist()) # Append coordinates of the 
                NameList.append(trans[1]) # Append the names of the atoms
    
        # Append the coordinates of the ligands 
        for index, entry in enumerate(LigandList):
            Xplot.append(entry[0])
            Yplot.append(entry[1])
            Zplot.append(entry[2])  
    
        LigandConstituent = [atom.name for atom in Ligand]
        Ligands = [] # TODO 
        for index in range(0, len(Sphere)): 
            Ligands = Ligands + LigandConstituent

        SphereName = [] 
        # Append the coordinates of the sphere 
        for entry in Sphere:
            XplotSphere.append(entry[0])
            YplotSphere.append(entry[1])
            ZplotSphere.append(entry[2])
            SphereName.append('Au')
    
        dfLigand = pd.DataFrame(list(zip(Xplot, Yplot, Zplot, Ligands)), columns =['X', 'Y', 'Z', 'NAME'])
        dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
        dfLigand['name'] = LigandName
        dfCore['name'] = CoreName
        
        return dfLigand, dfCore
    
    
    def ReturnPandasNPMartini(Molecule, LigandAlignmentVector, TransformationList, SphereList, LigandName, 
                          CoreName, Length = 1.0):
        """ Function to read Martini molecule information and orientate on NP surface 
        
        Args:
            Placeholder
        Returns: 
            Placeholder
        Raises:
            Placeholder 
        """
        LigandList, NameList = [], []
        Sphere = []
        Xplot, Yplot, Zplot = [], [], []
        XplotSphere, YplotSphere, ZplotSphere = [], [], []

        # Sulfur/ligand vector 
        unitVector = np.linalg.norm(LigandAlignmentVector)
        vec1 = LigandAlignmentVector.tolist()
        for index in range(0, len(SphereList)):
            vec2 = SphereList[index] 
            # Find the rotation matrix that aligns ligand vector representation to the NP surface vector point representation 
            TransformationVector = rotation_matrix_from_vectors(vec1, vec2)  
            vec1_rot = TransformationVector.dot(vec1) # Rotate the vector to match the surface point on the sphere 
            # TODO 
            unitVectorabs = np.linalg.norm(LigandAlignmentVector)  
            #vecMultiplier = vec1_rot * (np.linalg.norm(np.array(vec2))) # Controls how far we want the ligands to be placed away from
                                               # the NP surface
            #vecMultiplier = [1.0, 1.0, 1.0]
            vecMultiplier = vec1_rot/unitVectorabs * (np.linalg.norm(np.array(vec2))) + (vec1_rot/unitVectorabs * Length)
            #Sphere.append(vec1_rot.tolist())
            Sphere.append(vec2)
            # Get the factors to translate the vector 
            for trans in TransformationList:
                #if vec1_rot[0] > 0 and vec1_rot[1] > 0 and vec1_rot[2] > 0:   # (+ , +, +)
                LigandAtomcoordinate = TransformationVector.dot(trans[0])
                LigandAtomcoordinate[0] = LigandAtomcoordinate[0] + vecMultiplier[0]
                LigandAtomcoordinate[1] = LigandAtomcoordinate[1] + vecMultiplier[1]
                LigandAtomcoordinate[2] = LigandAtomcoordinate[2] + vecMultiplier[2]
                LigandList.append(LigandAtomcoordinate.tolist()) # Append coordinates of the 
                NameList.append(trans[1]) # Append the names of the atoms 

        # Append the coordinates of the ligands 
        for index, entry in enumerate(LigandList):
            #ax.plot3D(entry[0], entry[1], entry[2], 'red')
            Xplot.append(entry[0])
            Yplot.append(entry[1])
            Zplot.append(entry[2])
        
        # Add in the ligand index 
        LigandConstituent = [atom.name for atom in Molecule] # Molecule is utilized here 
        Ligands = []
        for index in range(0, len(Sphere)): 
            Ligands = Ligands + LigandConstituent
    
        SphereName = [] 
        # Append the coordinates of the sphere 
        for entry in Sphere:
            #ax.plot3D(entry[0], entry[1], entry[2], 'red')
            XplotSphere.append(entry[0])
            YplotSphere.append(entry[1])
            ZplotSphere.append(entry[2])
            SphereName.append('Au')
    
        dfLigand = pd.DataFrame(list(zip(Xplot, Yplot, Zplot, Ligands)), columns =['X', 'Y', 'Z', 'NAME'])
        dfLigand['name'] = LigandName
        dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
        dfCore['name'] = CoreName
        Total = dfLigand.append(dfCore)
        return Total
    
    def AttachLigands(LigandSmilesString, FirstAtomList, LastAtomList, SphereList, Length = 1.0, option = 'Plain'):
        """ This function is utilized to place the smiles description of the ligand(s). 
    
        e.g use for the case of a plain ligand functionalized NP: 
    
        LigandSmilesStrings = ['C1=C(C=CC=C1)CS[H]']
        PandasNPDataframe = AttachLigands(LigandSmilesStrings, ['S7'], ['C4'], SphereList)

        e.g use for the case of a Janus or Striped NP:
    
        LigandSmilesStrings = ['C1=C(C=CC=C1)CS[H]']
        PandasNPDataframe = AttachLigands(LisphegandSmilesStrings, ['S7', 'S0'], ['C4', 'O3'], SphereList)
       
        Args:
            LigandSmilesString: 
                smiles of the ligand we want to attach. We want to give a list 
                as we want to have the options of attaching the ligands onto 
                a Janus or Striped NP 
            FirstAtomList:
                Name of the first atoms in the ligands
            LastAtomList: 
                Name of the last atoms in the ligands 
            SphereList:
                The core which has been generated by GenerateCore function 
            option (default = 'plain'):
                 The type of NP we want to generate. 
         
        Returns:
        
            
        Raises:
    
        """
        # If the option is 'Plain' - we create a simple Nanoparticle without any patterns included  
        if option == 'Plain':    
            Ligand_I, Core_I = ReturnPandasNP(LigandSmilesString[0], FirstAtomList[0], 
                                              LastAtomList[0], SphereList[0], 'Ligand1', 'Core', Length)
            Total = Core_I.append(Ligand_I)
            return Total 
    
        # If the option is 'Janus' or 'Striped', then we have to include the different 
        # ligands.  We have two entries of ligands we need to take into account 
        elif option == 'Janus' or option == 'Striped':
        
            # If the option is 'Janus' or 'Striped', then we have to include the different 
            # ligands.  We have two entries of ligands we need to take into account 
        
            Ligand_I, Core_I = ReturnPandasNP(LigandSmilesString[0], FirstAtomList[0], 
                                              LastAtomList[0], SphereList[0], 'Ligand1', 'Core', Length)
        
            Ligand_II, Core_II = ReturnPandasNP(LigandSmilesString[1], FirstAtomList[1], 
                                              LastAtomList[1], SphereList[1], 'Ligand2', 'Core', Length)
            # Append Core with Ligands
            MainCore = Core_I.append(Core_II)
            # Add index to the core 
            Ligands = Ligand_I.append(Ligand_II)
            # Add index to the ligands 
            #Ligands['index'] = range(1, len(Ligands) + 1)
            Total = MainCore.append(Ligands)
            Total = Total.reset_index()
            Total['index_col'] = Total.index
            return Total
        
        
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
        
    def itp():
        """
        """
        MartiniAltItp_I = GromacsTopologyFile('/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.top')
        MartiniAltItp_II = GromacsTopologyFile('/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.top')
    
        indices_1 = [i for i in PandasNPDataframe[PandasNPDataframe['name'] == 'Lig1']['index'].iloc[::3]]
        # These indices need a variable to show the length of each ligand
        indices_2 = [i for i in PandasNPDataframe[PandasNPDataframe['name'] == 'Lig2']['index'].iloc[::4]]
    
        LigandStringBonds = [] 
        LigandStringImpropers = [] 
        # get bond parameters for first type of ligand on the NP
        LigandStringBonds.append('; ST - P5 ')
    
        for index in indices_1: 
            for bond in MartiniAltItp_I.bonds:
                bondstring = f"{bond.atom1.idx + (index + 1)} {bond.atom2.idx + (index + 1)} {bond.funct} {bond.type.req} {bond.type.k}"
                LigandStringBonds.append(bondstring)
            # get improper dihedral parameters 
            for improper in MartiniAltItp_I.impropers:
                dihedralstring = f"{improper.atom1.idx + (index + 1)} {improper.atom2.idx + (index + 1)} {improper.atom3.idx + (index + 1)} {improper.atom4.idx + (index + 1)} {improper.funct} {improper.type.psi_eq} {improper.type.psi_k}"
                LigandStringImpropers.append(dihedralstring)
        
        # get bond parameters for the second type of ligand on the NP 
        for index in indices_2:
            for bond in MartiniAltItp_II.bonds:
                bondstring = f"{bond.atom1.idx + (index + 1)} {bond.atom2.idx + (index + 1)} {bond.funct} {bond.type.req} {bond.type.k}"
                LigandStringBonds.append(bondstring)
            # get improper dihedral parameters   
            for improper in MartiniAltItp_II.impropers:
                dihedralstring = f"{improper.atom1.idx + (index + 1)} {improper.atom2.idx + (index + 1)} {improper.atom3.idx + (index + 1)} {improper.atom4.idx + (index + 1)} {improper.funct} {improper.type.psi_eq} {improper.type.psi_k}"
                LigandStringImpropers.append(dihedralstring)
        
        return LigandStringBonds, LigandStringImpropers
