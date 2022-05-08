#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
------------------------
Last Updated: 30/04/2022
------------------------
----------------------------------------
THIS NOTEBOOK NEEDS TO BE MORE ORGANIZED 
-----------------------------------------

The purpose of this notebook is to get a idea of how the translation of the smiles to 
the coarse-grained beads of Martini can work. 

-> https://cedric.bouysset.net/blog/2020/08/07/rdkit-interoperability
-> https://stackoverflow.com/questions/47319238/python-plot-3d-vectors
-> https://en.wikipedia.org/wiki/XYZ_file_format - information on the xyz coordination file 
-> https://mattermodeling.stackexchange.com/questions/3961/recalculate-atom-positions-to-account-for-periodic-boundary-conditions/3970#3970  - Atomsk
-> https://stackoverflow.com/questions/49064611/how-to-find-different-groups-in-networkx-using-python - grouping networks 
-> https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/
-> http://cgmartini.nl/index.php/component/kunena/8-martini-philosophy/5776-mapping-of-benzene-ring
-> https://docs.mdanalysis.org/1.0.0/documentation_pages/lib/NeighborSearch.html
-> https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
-> https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
-> https://nanotube.msu.edu/fullerene/fullerene-isomers.html
-> https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753
-> https://plotly.com/python/3d-scatter-plots/
-> https://undergroundmathematics.org/circles/cutting-spheres/solution
-> https://www.ossila.com/products/pc70bm - pcbm 
-> https://github.com/mosdef-hub/nanoparticle_optimization
-> https://pubs.acs.org/doi/10.1021/acs.jctc.8b01269

Quote by Riccardo: 
------------------
In the open beta of Martini 3, benzene is indeed mapped with three TC4 beads, that that's again a 2-to-1 mapping. 
The bond length is changed to 0.29 nm in 3.0 because this allows to represent more closely the volume of a benzene 
molecule, taking into account also the smaller size of T-beads as compared to S -beads

---------------------------
How the modules are divided
---------------------------

-> We want to be given free reign into making NPs of multiple types. The primary types we are concerned with at the 
   moment is 
   
   1. Functionalized AuNP type ones.  
   2. Carbon Nanotube like ones. 
   3. Large spherical buckyball type structures - Like C70. 

-> TODO  
   - Am not able to get the box dimenisons correct for the xyz file...
   - Need to write Striped/Janus functionality to the NP maker  

Trying to translate the atomic C70 structure into the coarse-grained c70 structure 

TODO - isolate per 2 beads - There is a two to one mapping, and find the center of mass 
       for each of these beads. Then construct the itp file for the bead. 

"""


# In[2]:


# Boilerplate libraries                                                                                                                                                                                    
import sys                                                                                                                                                                                                 
import re                                                                                                                                                                                                  
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import math 
from operator import itemgetter
import itertools                                                                                                                                                                                           
import requests                                                                                                                                                                                            
import collections                                                                                                                                                                                         
import random            

# Scipy libraries 
import scipy                                                                                                                                                                                                                                                                                                                                                                                                          
from scipy.sparse import csr_matrix                                                                                                                                                                        
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

# logging module 
import logging 

# Change logging configurations so that we are printing out 

pd.set_option('display.max_colwidth', None)
IPythonConsole.ipython_useSVG=True  #< set this to False if you want PNGs instead of SVGs


# In[6]:



"""
--------------------
Misc Functionalities
--------------------

"""
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

"""
Plotting functions 
"""

def mol_with_atom_index(mol):
    """
    
    """
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def vector_plot(tvects,is_vect=True,orig=[0,0,0]):
    """Plot vectors using plotly
    
    Args:
        Placeholder
    Returns:
        Placeholder
    Raises:
        Placeholder
    """
    if is_vect:
        if not hasattr(orig[0],"__iter__"):
            coords = [[orig,np.sum([orig,v],axis=0)] for v in tvects]
        else:
            coords = [[o,np.sum([o,v],axis=0)] for o,v in zip(orig,tvects)]
    else:
        coords = tvects
    data = []
    for i,c in enumerate(coords):
        X1, Y1, Z1 = zip(c[0])
        X2, Y2, Z2 = zip(c[1])
        vector = go.Scatter3d(x = [X1[0],X2[0]],
                              y = [Y1[0],Y2[0]],
                              z = [Z1[0],Z2[0]],
                              marker = dict(size = [0,5],
                                            color = ['blue'],
                                            line=dict(width=5,
                                                      color='DarkSlateGrey')),
                              name = 'Vector'+str(i+1))
        data.append(vector)

    layout = go.Layout(
             margin = dict(l = 4,
                           r = 4,
                           b = 4,
                           t = 4)
                  )
    fig = go.Figure(data=data,layout=layout)
    fig.show()
    

def FibanocciSphere(samples=1):
    """ Return a Fibanocci sphere with N number of points on the surface. 

    This will act as the template for the nanoparticle core. 
    
    Args:
        Placeholder
    Returns:
        Placeholder
    Raises:
        Placeholder
    """
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

"""
---------------------------------
Actual NP construction functions 
---------------------------------
"""
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
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

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
    if Option == 'Plain':
        return [SphereList[1:-1]]
    elif Option == 'Striped':
        StripedValues, CeilingValues = LabelStripedNP(SphereList[1:-1], Option)
        return StripedValues, CeilingValues
    elif Option == 'Janus':
        TopValues, BottomValues = LabelStripedNP(SphereList[1:-1], Option)
        return TopValues, BottomValues 
    
def ReturnPandasNP(LigandString, FirstAtom, LastAtom, SphereList, LigandName, CoreName):
    """Placeholder
    
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
    
    for index in range(0, len(SphereList)):
        vec2 = SphereList[index]
        TransformationVector = rotation_matrix_from_vectors(vecLigand, vec2)
        vec1_rot = TransformationVector.dot(vecLigand) # Rotate the vector to match the surface point on the sphere 
        unitVectorabs = np.linalg.norm(LigandAlignmentVector)
        vecMultiplier = vec1_rot/unitVectorabs * (np.linalg.norm(np.array(vec2))) 
        Sphere.append(vec2)
        
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
        #ax.plot3D(entry[0], entry[1], entry[2], 'red')
        XplotSphere.append(entry[0])
        YplotSphere.append(entry[1])
        ZplotSphere.append(entry[2])
        SphereName.append('Au')
    
    dfLigand = pd.DataFrame(list(zip(Xplot, Yplot, Zplot, Ligands)), columns =['X', 'Y', 'Z', 'NAME'])
    dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
    dfLigand['name'] = LigandName
    dfCore['name'] = CoreName
    return dfLigand, dfCore


def ReturnPandasNPMartini(Molecule, LigandAlignmentVector, TransformationList, SphereList, LigandName, CoreName):
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
        
        vecMultiplier = vec1_rot/unitVectorabs * (np.linalg.norm(np.array(vec2)))
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
    



def AttachLigands(LigandSmilesString, FirstAtomList, LastAtomList, SphereList, option = 'Plain'):
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
                                          LastAtomList[0], SphereList[0], 'Ligand1', 'Core')
        Total = Core_I.append(Ligand_I)
        return Total 
    # If the option is 'Janus' or 'Striped', then we have to include the different 
    # ligands.  We have two entries of ligands we need to take into account 
    elif option == 'Janus' or option == 'Striped':
        
        # If the option is 'Janus' or 'Striped', then we have to include the different 
        # ligands.  We have two entries of ligands we need to take into account 
        
        #ReturnPandasNP(LigandString, FirstAtom, LastAtom, SphereList, LigandName, CoreName):
        Ligand_I, Core_I = ReturnPandasNP(LigandSmilesString[0], FirstAtomList[0], 
                                          LastAtomList[0], SphereList[0], 'Ligand1', 'Core')
        Ligand_II, Core_II = ReturnPandasNP(LigandSmilesString[1], FirstAtomList[1], 
                                          LastAtomList[1], SphereList[1], 'Ligand2', 'Core')
        # Append Core with Ligands
        MainCore = Core_I.append(Core_II)
        Ligands = Ligand_I.append(Ligand_II)
        Total = MainCore.append(Ligands)
        
        return Total
        
def AttachLigandsOriginal(LigandSmilesString, FirstAtomList, LastAtomList, SphereList):
    """ The original Attach ligands function 
    
    This function takes a smiles strings, orientates 
    the ligand of that of the smiles string, and then positions them onto the 
    surface of the NP

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
    Returns:
        Placeholder:
    Raises:
        Placeholder: 
    """
    Molecule = u.select_atoms('all')
    # Select Atom attached to the core 
    FirstAtom = u.select_atoms('name {}'.format(FirstAtomList[0])) # Pick out the atoms attached
    # Select end atom on the ligand 
    LastAtom = u.select_atoms('name {}'.format(LastAtomList[0])) # Pick out the last atom at the tip of the ligand 
    TransformationList, NameList = [], []
    # Find the vector representing the direction from the sulfur to the tip of the ligand. 
    LigandAlignmentVector = (FirstAtom.positions- LastAtom.positions)[0]  
     
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
            
    LigandList = [] 
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
        vecMultiplier = vec1_rot/unitVectorabs * 3 # Controls how far we want the ligands to be placed away from
                                               # the NP surface
        # TODO
        Sphere.append(vec1_rot.tolist())
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
    LigandConstituent = [atom.name for atom in Molecule]
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
    dfLigand['name'] = 'Ligand'
    dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
    dfCore['name'] = 'Core'
    Total = dfLigand.append(dfCore)
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
    #GroPath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/gros"
    #ItpPath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono"
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

        Coordinates = ReturnPandasNPMartini(Molecule, LigandAlignmentVector, 
                                              TransformationList, SphereList[0])
        return Coordinates
    
    if option == 'Janus' or option == 'Striped':
        # Information here 
        # First ligand 
        Molecule_I, TransformationList_I, LigandAlignmentVector_I = ReadMartiniMolecules(
            GroFiles[0], FirstAtoms[0], LastAtoms[0])
        # Second Ligand 
        Molecule_II, TransformationList_II, LigandAlignmentVector_II = ReadMartiniMolecules(
            GroFiles[1], FirstAtoms[1], LastAtoms[1])
    
        Coordinates_I = ReturnPandasNPMartini(Molecule_I, LigandAlignmentVector_I, 
                                              TransformationList_I, SphereList[0], 'Lig1', 'Core')
        Coordinates_II = ReturnPandasNPMartini(Molecule_II, LigandAlignmentVector_II, 
                                              TransformationList_II, SphereList[1], 'Lig2', 'Core')
        
        Coordinates = Coordinates_I.append(Coordinates_II)
        return Coordinates 
    
    
def AttachLigandsMartiniOriginal(GroFiles, FirstAtoms, LastAtoms, SphereList, option = 'Plain'):
    """
    
    """
    Molecule, TransformationList, LigandAlignmentVector = ReadMartiniMolecules(
        GroFiles[0], FirstAtom[0], LastAtom[0])
    
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
        vecMultiplier = vec1_rot/unitVectorabs * 2 # Controls how far we want the ligands to be placed away from
                                               # the NP surface
        # TODO
        Sphere.append(vec1_rot.tolist())
        #LigandList.append(vec1_rot.tolist())
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
    LigandConstituent = [atom.name for atom in Molecule]
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
    dfLigand['name'] = 'Ligand'
    dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
    dfCore['name'] = 'Core'
    Total = dfLigand.append(dfCore)
    return Total
    
def ComputeC70Distances(coordinate):
    """    
    
    Args:
        coordinate: 
            path and file where the pdb file is.
    Returns: 
        
    Raises:
    
    
    """
    u = mda.Universe(coordinate)
    
    CorePositions, CoreIndex = [], []
    C70CGBeadPositions = []
    DistanceDict = {}
    CorePositions = [[index, atoms.position] for index, atoms in enumerate(u.atoms)]
    
    for index, atoms in enumerate(u.atoms): # iterate through items 
        DistanceDict[index] = []
        # Get the distance between atoms in the core positions  
        for items in CorePositions:
            # Find the distance between the index atoms in the index, atom.. 
            # and the core 
            dist = distance.euclidean(items[1], atoms.position)
            Entry = [index, items[0], dist]
            # Sort entry by dist 
            DistanceDict[index].append(Entry)
    ClosestAtomsDistance = [] # List to store the atoms with the closest distances
    # Sort each entry by distance 
    
    for key in DistanceDict.keys():
        DistanceDict[key] = sorted(DistanceDict[key], key=itemgetter(2)) # sort entries by closest distance 
        DistanceDict[key][1][0:2] = sorted(DistanceDict[key][1][0:2]) # 
        ClosestAtomsDistance.append(DistanceDict[key][1])
    # Remove duplicate entries 
    UniqueSetsCoordinates = [list(x) for x in set(tuple(x) for x in ClosestAtomsDistance)]       
    # Take the indices in the unique data and compute the averge coordinates
    
    for entry in UniqueSetsCoordinates:
        data = [list(CorePositions[entry[0]][1]), list(CorePositions[entry[1]][1])]
        averaged = np.average(data, axis=0)
        C70CGBeadPositions.append(averaged)
        
    return C70CGBeadPositions


def GenerateXYZ(NPDataframe, name):
    """ Generate an XYZ file from the pandas dataframe 
        we generated from the previous functions 
    Args:
        NPDataframe:
            Placeholder
        name: 
            Placeholder 
    Returns:
        Placeholder
    Raises::
        Placeholder 
        
    """
    with open(name, 'w') as xyz_file:
        xyz_file.write("%d\n%s\n" % (len(NPDataframe), 'NP'))
        for index, row in NPDataframe.iterrows():
            xyz_file.write("{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                row['NAME'], row['X'], row['Y'], row['Z']))
    


# In[7]:


# Generating Striped configuration 
PlainSphereList = GenerateCore(10, 40)
StripedLigandSmileStrings = ["SCCO[H]", "C1=C(C=CC=C1)CS[H]"]
StripedPandasNPDataframe = AttachLigands(StripedLigandSmileStrings, ['S0', 'S7'], ['O3', 'C4'], PlainSphereList)
fig = px.scatter_3d(StripedPandasNPDataframe, x='X', y='Y', z='Z', color='name')
fig.show()


# In[8]:


StripedPandasNPDataframe.to_pickle('stripedPickle')


# In[10]:


# Generating Striped configuration 
StripedSphereList = GenerateCore(10, 40, 'Striped')
StripedLigandSmileStrings = ["SCCO[H]", "C1=C(C=CC=C1)CS[H]"]
StripedPandasNPDataframe = AttachLigands(StripedLigandSmileStrings, ['S0', 'S7'], ['O3', 'C4'], StripedSphereList, option = 'Striped')
fig = px.scatter_3d(StripedPandasNPDataframe, x='X', y='Y', z='Z', color='name')
fig.show()


# In[11]:


# Generating Janus configuration 
JanusSphereList = GenerateCore(10, 30, 'Janus')
JanusLigandSmileStrings = ["SCCO[H]", "C1=C(C=CC=C1)CS[H]"]
JanusPandasNPDataframe = AttachLigands(JanusLigandSmileStrings, ['S0', 'S7'], ['O3', 'C4'], JanusSphereList, option = 'Janus')
fig = px.scatter_3d(JanusPandasNPDataframe, x='X', y='Y', z='Z', color='name')
fig.show()


# In[13]:


""" 
Annotating and drawing an example sulfur ligand - The first ligand we want to add 
"""
u1 = mda.Universe.from_smiles("C1=C(C=CC=C1)CS[H]")
AromaticSulfurSmilesString = 'C1=C(C=CC=C1)CS[H]'
m_aromatic = Chem.MolFromSmiles(AromaticSulfurSmilesString)
#mol_with_atom_index(m_aromatic)
m_aromatic


# In[17]:


""" 
Annotating and drawing an example sulfur ligand - The The second ligand we want to add in case we are adding in 
a Janus/Striped Nanoparticle 


"""
u2 = mda.Universe.from_smiles("SCCO[H]")
NonAromaticSulfurSmilesString_II = 'SCCO[H]'
m_nonaromatic = Chem.MolFromSmiles(NonAromaticSulfurSmilesString_II)
#mol_with_atom_index(m_aromatic)
m_nonaromatic


# In[20]:


for i in u2.atoms:
    print(i)


# In[14]:


""" 
Working on the same ligand.. 
"""
# Through RDKIT
AromaticSulfurSmilesString = 'C1=C(C=CC=C1)CS[H]'
m_aromatic = Chem.MolFromSmiles(AromaticSulfurSmilesString)
mol_with_atom_index(m_aromatic)

# Through MDAnalysis 
u1 = mda.Universe.from_smiles(AromaticSulfurSmilesString)
Molecule = u1.select_atoms('all')
Molecule.positions # Finds the cartesian coordinates of the ligands 


# In[11]:


C70 = mda.Universe.from_smiles("c1(c2c3c4c15)c6c7c8c2c9c%10c3c%11c%12c4c%13c%14c5c%15c6c%16c7c%17c%18c%19c%20c%21c%22c%23c%24c%21c%25c%26c%20c%18c%16c%27c%15c%14c%28c(c%25c%29c%24c%30c%31c%23c%32c%33c%22c%19c%34c%33c(c9c8c%34%17)c%35c%10c%11c(c%31c%32%35)c%36c%12c%13c%28c%29c%30%36)c%26%27")


# In[14]:


C70


# In[32]:


C70.atoms


# In[33]:


# Highlight a Substructure in a Molecule - Can we identify the benzene part and sulfur part?
ConvertedSmiles = Chem.MolToSmiles(Chem.MolFromSmiles('C1=CC=CN=C1'))
ConvertedSmiles
ExampleMartiniString = "CC(=O)CO" # Example smiles string that is compatible with MARTINI 
ConvertedSmiles2 = Chem.MolToSmiles(Chem.MolFromSmiles(ExampleMartiniString))
ConvertedSmiles2
m2  = Chem.MolFromSmiles(ConvertedSmiles2)
# Hence, we need to create a dictionary that catalogues the string with the relevant martini bead 
SmilesToMartiniDictionary = {}
SmilesToMartiniDictionary["CC(=O)CO"] = 'P2' # P2 Bead 
SmilesToMartiniDictionary["CC(=O)O"] = 'SP2' # SP2 Bead 
SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 
SmilesToMartiniDictionary["CC(C)O"] = 'P1' # P1 Bead 


# In[34]:


"""
This part links the lignad coordinates with the 
"""
# new feature
u1 = mda.Universe.from_smiles("c1ncncc1C(=O)[O-]")
# new feature
Molecule = u1.select_atoms('all')
Molecule.positions # Finds 
# Need to label each of the xyz coordinates with the relevant indices within the 


# In[35]:


m = Chem.MolFromSmiles('c1cc(C(=O)O)c(OC(=O)C)cc1')
substructure = Chem.MolFromSmarts('CC(=O)O')
print(m.GetSubstructMatches(substructure)) # Shows the indices of the molcule that matches the MARTINI bead 

# Find part of the smiles which is part of a ring 

ri = m.GetRingInfo().AtomRings()
for ring in ri:
    print(ring)


# In[36]:


"""
In the above situation, we can identify the 3 parts - the thiolated part, the saturated hydrocarbon chain attached 
to the benzene group, and the benzene group itself. 
"""
substructure = Chem.MolFromSmarts('SC')
print(m.GetSubstructMatches(substructure))


# In[15]:


vector_plot(GenerateCore(10, 100))


# In[44]:


"""
Generate an all-atomic based NP 
"""
SphereList = GenerateCore(30, 20)
AromaticSulfurSmilesString = 'C1=C(C=CC=C1)CS[H]'
PandasNPDataframe = AttachLigands('C1=C(C=CC=C1)CS[H]', ['S7'], ['C4'], SphereList)
#PandasNPDataframe = AttachLigandsMartini(MIMILigandPath, 'N1', 'R3', SphereList)
PandasNPDataframe


# In[71]:


"""

Generating a Martini based NP - works both for Striped and Janus

"""
GroPath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/gros" # self-explanatory 
# ItpPath = "/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono"  # not quite as self-explanatory
MIMILigandPath = GroPath + "/" + "1MIMI.gro" # First ligand 
NITLLigandPath = GroPath + "/" + "2NITL.gro" # Second ligand 
SphereList = GenerateCore(10, 40, "Striped") # Generate the NP core 
PandasNPDataframe = AttachLigandsMartini([MIMILigandPath, NITLLigandPath], ['N1', 'N1'] , ['R3', 'R4'], 
                                         SphereList, 'Striped')

# Plotting the 3D plot of the new Martini NP 
fig = px.scatter_3d(PandasNPDataframe, x='X', y='Y', z='Z', color='name')
fig.show()


# In[1]:


GenerateXYZ(PandasNPDataframe, 'MARTININP.xyz')


# In[28]:


set(PandasNPDataframe['name'])


# In[21]:


m = Chem.MolFromMolFile('c70.mol') # Read in the C70 mol 
mol1 = u1.atoms.convert_to("RDKIT") # Convert to rdkit 
u = mda.Universe('c70.pdb') # 
mol1 = u.atoms.convert_to("RDKIT", force=True)


# In[199]:


def itpConstraintParser():
    """
    Args:
        Placeholder
    Returns:
        Placeholder
    Raises:
        Placeholder 
    """
    file1 = open('/home/sang/Desktop/GIT/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.itp', 'r')
    Lines = file1.readlines()
    for i,j  in enumerate(Lines):
        if '[constraints]' in j:
            print(i,j) 


# In[179]:


gromacs.fileformats.__dict__


# In[22]:


Molecule = mda.Universe('c70.pdb')
mol1 = u.atoms.convert_to("RDKIT", force=True)
C70List = [] 
for index, atoms in enumerate(u.atoms):
    C70List.append(list(atoms.position))
C70List


# In[59]:


def ConvertC70Core(coordinate):
    """    
    Description 
    -----------
    
    Creates a MDAnalysis universe from the pdb input of the C70 file, and 
    creates a coarse-grained framework of the C70. 
    
    
    Parameters
    ----------
    coordinate: 
        path and file where the pdb file is.
    filename : str
    copy : bool
    dtype : data-type
    iterable : iterable object
    shape : int or tuple of int
    files : list of str
    
    Returns
    -------
    int
        Description of anonymous integer return value.
    """
    u = mda.Universe(coordinate) # Create Universe 
    CorePositions, CoreIndex = [], [] # Get positions and indices of the core 
    C70CGBeadPositions = [] # list to append the new 'coarse-grained' positions of the C70 Martini NP
    DistanceDict = {} 
    
    CorePositions = [[index, atoms.position] for index, atoms in enumerate(u.atoms)]
    
    for index, atoms in enumerate(u.atoms): # loop over the indices and atoms within the universe 
        DistanceDict[index] = []
        # Get the distance between atoms in the core positions  
        for items in CorePositions: # each element here is [index, atom.position]
            Distance = distance.euclidean(items[1], atoms.position) # get the distance between the atoms and the core
            Entry = [index, items[0], Distance] # index of position A and index of position B, and the distance between 
            DistanceDict[index].append(Entry) # Append this list 
    ClosestAtomsDistance = [] # List to store the atoms with the closest distances
    
    # Sort each entry by distance 
    for key in DistanceDict.keys():
        DistanceDict[key] = sorted(DistanceDict[key], key=itemgetter(2)) # sort entries by closest distance 
        DistanceDict[key][1][0:2] = sorted(DistanceDict[key][1][0:2]) # Get first two index 
        ClosestAtomsDistance.append(DistanceDict[key][1]) 
    
    # Remove duplicate entries 
    UniqueSetsCoordinates = [list(x) for x in set(tuple(x) for x in ClosestAtomsDistance)]       
    
    # Take the indices in the unique data and compute the averge coordinates
    for entry in UniqueSetsCoordinates:
        data = [list(CorePositions[entry[0]][1]), list(CorePositions[entry[1]][1])]
        averaged = np.average(data, axis=0)
        C70CGBeadPositions.append(averaged)
        
        
    return C70CGBeadPositions

def ConstructNPNetwork():
    """
    """
    pass


# In[60]:


NPPositions = ConvertC70Core('c70.pdb')
C70 = [[i[0], i[1], i[2]] for i in NPPositions]
C70


# In[61]:


"""
Plotting the new C70 core and plotting the 3D scatter plot to see what it looks like
"""
df = pd.DataFrame(NPPositions, columns = ['X','Y','Z'])
fig = px.scatter_3d(df, x='X', y='Y', z='Z')
#fig = px.scatter_3d(dfCore, x='X', y='Y', z='Z')
fig.show()


# ![a-Cross-sectional-SEM-image-of-the-device-architecture-b-chemical-structures-of.png](attachment:a-Cross-sectional-SEM-image-of-the-device-architecture-b-chemical-structures-of.png)

# In[40]:


"""
The smiles for the ligand 

Links on information about the C60/C70 ligand bonded NP
-------------------------------------------------------

poly(3-hexyl-thiophene) (P3HT)

phenyl-C-61-butyric acid methyl ester (PCBM).

https://pubchem.ncbi.nlm.nih.gov/compound/3-Hexylthiophene

https://pubchem.ncbi.nlm.nih.gov/compound/6_6_-Phenyl-C61-butyric-acid-methyl-ester

https://chemrxiv.org/engage/chemrxiv/article-details/60c74b9abb8c1a44713db259

https://research.rug.nl/en/publications/resolving-donor-acceptor-interfaces-and-charge-carrier-energy-lev

https://www.americanelements.com/60-pcbm-6-6-phenyl-c61-butyric-acid-methyl-ester-160848-22-6

http://www-jmg.ch.cam.ac.uk/data/molecules/misc/c70.html
 
https://en.wikipedia.org/wiki/Phenyl-C61-butyric_acid_methyl_ester - information on its application
as potential solar cells 

"""


# Which ligand do I need to attach by? 


P3HT = 'CCCCCCC1=CSC=C1'
m_P3HT = Chem.MolFromSmiles(P3HT)
# Universe 
u_P3HT = mda.Universe.from_smiles(P3HT)
# Need to add index 
u_P3HT_object = u_P3HT.select_atoms('all')
#mol_with_atom_index(m_aromatic)

"""
<Atom 1: C0 of type CA resid 1 and segid SYSTEM>
<Atom 2: C1 of type C resid 1 and segid SYSTEM>
<Atom 3: C2 of type C resid 1 and segid SYSTEM>
<Atom 4: C3 of type C resid 1 and segid SYSTEM>
<Atom 5: C4 of type C resid 1 and segid SYSTEM>
<Atom 6: C5 of type C resid 1 and segid SYSTEM>
<Atom 7: C6 of type C resid 1 and segid SYSTEM>
<Atom 8: C7 of type C resid 1 and segid SYSTEM>
<Atom 9: S8 of type S resid 1 and segid SYSTEM>
<Atom 10: C9 of type C resid 1 and segid SYSTEM>
<Atom 11: C10 of type C resid 1 and segid SYSTEM>
<Atom 12: H11 of type H resid 1 and segid SYSTEM>
<Atom 13: H12 of type H resid 1 and segid SYSTEM>
<Atom 14: H13 of type H resid 1 and segid SYSTEM>
<Atom 15: H14 of type H resid 1 and segid SYSTEM>
<Atom 16: H15 of type H resid 1 and segid SYSTEM>
<Atom 17: H16 of type H resid 1 and segid SYSTEM>
<Atom 18: H17 of type H resid 1 and segid SYSTEM>
<Atom 19: H18 of type H resid 1 and segid SYSTEM>
<Atom 20: H19 of type H resid 1 and segid SYSTEM>
<Atom 21: H20 of type H resid 1 and segid SYSTEM>
<Atom 22: H21 of type H resid 1 and segid SYSTEM>
<Atom 23: H22 of type H resid 1 and segid SYSTEM>
<Atom 24: H23 of type H resid 1 and segid SYSTEM>
<Atom 25: H24 of type H resid 1 and segid SYSTEM>
<Atom 26: H25 of type H resid 1 and segid SYSTEM>
<Atom 27: H26 of type H resid 1 and segid SYSTEM>


Testing attaching ligands onto C70 
"""
PCBM = 'C1=CC=CC(=C1)CCCCC(=O)OC'
m_PCBM = Chem.MolFromSmiles(PCBM)
u_PCBM = mda.Universe.from_smiles(PCBM)


LigandList = [] 
Sphere = []
Xplot, Yplot, Zplot = [], [], []
XplotSphere, YplotSphere, ZplotSphere = [], [], []
# Select Atom attached to the core 
#FirstAtom = u_P3HT.select_atoms('name {}'.format('C0')) # Pick out the atoms attached
# Select end atom on the ligand 
#LastAtom = u_P3HT.select_atoms('name {}'.format('S8')) # Pick out the last atom at the tip of the ligand 
#FirstAtom = u_P3HT.select_atoms('name {}'.format('S8')) # Pick out the atoms attached
# Select end atom on the ligand 
#LastAtom = u_P3HT.select_atoms('name {}'.format('C0')) # Pick out the last atom at the tip of the ligand 
FirstAtom = u_PCBM.select_atoms('name {}'.format('C6')) # Pick out the atoms attached
# Select end atom on the ligand 
LastAtom = u_PCBM.select_atoms('name {}'.format('C10')) # Pick out the last atom at the tip of the ligand 
TransformationList, NameList = [], []
# Find the vector representing the direction from the sulfur to the tip of the ligand. 
LigandAlignmentVector = (FirstAtom.positions- LastAtom.positions)[0]
Molecule = u_PCBM.select_atoms('all')
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
    
unitVector = np.linalg.norm(LigandAlignmentVector)
vec1 = LigandAlignmentVector.tolist()
for index in range(0, len(C70)):
    vec2 = C70[index] # get cartesian coordinates of the core atom 
    # Find the rotation matrix that aligns ligand vector representation to the NP surface vector point representation 
    TransformationVector = rotation_matrix_from_vectors(vec1, vec2)  
    vec1_rot = TransformationVector.dot(vec1) # Rotate the vector to match the surface point on the sphere 
    # TODO 
    unitVectorabs = np.linalg.norm(LigandAlignmentVector)  
    vecMultiplier = vec1_rot/unitVectorabs * 4 # Controls how far we want the ligands to be placed away from
                                               # the NP surface
    # TODO
    Sphere.append(vec1_rot.tolist())

#LigandList.append(vec1_rot.tolist())
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
LigandConstituent = [atom.name for atom in Molecule]
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
dfLigand['name'] = 'Ligand'
dfCore = pd.DataFrame(list(zip(XplotSphere, YplotSphere, ZplotSphere, SphereName)), columns =['X', 'Y', 'Z', 'NAME'])
dfCore['name'] = 'Core'
Total = dfLigand.append(dfCore)


# In[77]:


Total


# In[78]:


fig = px.scatter_3d(Total, x='X', y='Y', z='Z', color='name')
#fig = px.scatter_3d(dfCore, x='X', y='Y', z='Z')
fig.show()


# In[51]:


GenerateXYZ(Total, 'new.xyz')


# In[34]:


vector


# In[23]:


for i in u_P3HT.atoms:
    print(i)


# In[36]:


SphereList = GenerateCore(10, 20)
SphereList


# In[38]:


# C70
Phenyl_C61_Butyric_Acid_Methyl_Ester = 'C1=CC=CC(=C1)CCCCC(=O)OC'
C70BAME = Chem.MolFromSmiles(Phenyl_C61_Butyric_Acid_Methyl_Ester)
u_C70BAME = mda.Universe.from_smiles(Phenyl_C61_Butyric_Acid_Methyl_Ester)
# Need to add index 
u_C70BAME_object = u_C70BAME.select_atoms('all')
#u_C70BAME_object.write("u_C70BAME.xyz")
mol_with_atom_index(C70BAME) # C6 is the attaching ligand 


# In[41]:


PCBM = 'C1=CC=CC(=C1)CCCCC(=O)OC'
m_PCBM = Chem.MolFromSmiles(PCBM)
mol_with_atom_index(m_PCBM)

#u_PCBM = mda.Universe.from_smiles(PCBM)
#for i in u_PCBM.select_atoms('all'):
#    print(i)
    


# In[17]:


m_PCBM


# In[ ]:


# https://stackoverflow.com/questions/53075481/how-do-i-cluster-a-list-of-geographic-points-by-distance


import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import math

points = np.array([[33.    , 41.    ],
       [33.9693, 41.3923],
       [33.6074, 41.277 ],
       [34.4823, 41.919 ],
       [34.3702, 41.1424],
       [34.3931, 41.078 ],
       [34.2377, 41.0576],
       [34.2395, 41.0211],
       [34.4443, 41.3499],
       [34.3812, 40.9793]])


def distance(origin, destination): #found here https://gist.github.com/rochacbruno/2883505
    lat1, lon1 = origin[0],origin[1]
    lat2, lon2 = destination[0],destination[1]
    radius = 6371 # km
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1))         * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d

def create_clusters(number_of_clusters,points):
    kmeans = KMeans(n_clusters=number_of_clusters, random_state=0).fit(points)
    l_array = np.array([[label] for label in kmeans.labels_])
    clusters = np.append(points,l_array,axis=1)
    return clusters

def validate_solution(max_dist,clusters):
    _, __, n_clust = clusters.max(axis=0)
    n_clust = int(n_clust)
    for i in range(n_clust):
        two_d_cluster=clusters[clusters[:,2] == i][:,np.array([True, True, False])]
        if not validate_cluster(max_dist,two_d_cluster):
            return False
        else:
            continue
    return True

def validate_cluster(max_dist,cluster):
    distances = cdist(cluster,cluster, lambda ori,des: int(round(distance(ori,des))))
    print(distances)
    print(30*'-')
    for item in distances.flatten():
        if item > max_dist:
            return False
    return True

if __name__ == '__main__':
    for i in range(2,len(points)):
        print(i)
        print(validate_solution(20,create_clusters(i,points)))


# In[96]:


DistIndex1Index2
UniqueIndex1Index2Distances = [list(i) for i in list(set(map(tuple, DistIndex1Index2)))]


# In[9]:


import os                                                                                                                                                                                                  
import numpy as np                                                                                                                                                                                         
import itertools                                                                                                                                                                                           
import requests                                                                                                                                                                                            
import collections                                                                                                                                                                                         
import random                                                                                                                                                                                              
                                                                                                                                                                                                           
# RDKit libaries                                                                                                                                                                                           
from rdkit import Chem                                                                                                                                                                                     
from rdkit.Chem import AllChem                                                                                                                                                                             
from rdkit.Chem import ChemicalFeatures                                                                                                                                                                   
from rdkit.Chem import rdchem                                                                                                                                                                              
from rdkit.Chem import rdMolDescriptors                                                                                                                                                                    
from rdkit import RDConfig                                                                                                                                                                                 
                                                                                                                                                                                                           
# Boilerplate libraries                                                                                                                                                                                    
import sys                                                                                                                                                                                                 
import re                                                                                                                                                                                                  
import math                                                                                                                                                                                                
import scipy                                                                                                                                                                                               
                                                                                                                                                                                                           
# scipy libaries                                                                                                                                                                                           
from scipy.sparse import csr_matrix                                                                                                                                                                        
from scipy.sparse.csgraph import floyd_warshall                                                                                                                                                            
from scipy.spatial import ConvexHull, convex_hull_plot_2d           

def read_DG_data(DGfile):                                                                                                                                                                                  
    # Reads Delta G_OW for fragments into dictionary                                                                                                                                                       
    DG_data = {}                                                                                                                                                                                           
    with open(DGfile) as f:                                                                                                                                                                                
        for line in f:                                                                                                                                                                                     
            (key,val) = line.rstrip().split()                                                                                                                                                              
            DG_data[key] = float(val)                                                                                                                                                                      
                                                                                                                                                                                                           
    return DG_data                                                                                                                                                                                         
                                                                                                                                                                                                           
def include_weights(A,w):                                                                                                                                                                                  
    # Weights atoms by setting diagonal components                                                                                                                                                         
    A_weighted = np.copy(A)                                                                                                                                                                                
    for i,weight in enumerate(w):                                                                                                                                                                          
        A_weighted[i,i] = weight     
        
def get_smarts_matches(mol):                                                                                                                                                                               
    #Get matches to SMARTS strings                                                                                                                                                                         
    smarts_strings = {                                                                                                                                                                                     
    'S([O-])(=O)(=O)O'  :    'Qa',                                                                                                                                                                         
    'S([O-])(=O)(=O)[C;!$(*F)]'   :    'Q0'}                                                                                                                                                               
    matched_maps = []                                                                                                                                                                                      
    matched_beads = []                                                                                                                                                                                     
    for smarts in smarts_strings:                                                                                                                                                                          
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))                                                                                                                                      
        for match in matches:                                                                                                                                                                              
            matched_maps.append(list(match))                                                                                                                                                               
            matched_beads.append(smarts_strings[smarts])                                                                                                                                                   
                                                                                                                                                                                                           
    return matched_maps,matched_beads           


# In[130]:


DataPath = '/home/sang/Desktop/GIT/MDNPPackage/temp/cgparam/fragment_DGs.dat'


# In[131]:


read_DG_data(DataPath)


# In[ ]:


def write_gro(mol_name,bead_types,coords0,gro_name):                                                                                                                                                       
    #write gro file                                                                                                                                                                                        
    conf = mol.GetConformer(0)                                                                                                                                                                             
    with open(gro_name,'w') as gro:                                                                                                                                                                        
        gro.write('single molecule of {}\n'.format(mol_name))                                                                                                                                              
        gro.write('{}\n'.format(len(bead_types)))                                                                                                                                                          
        i = 1                                                                                                                                                                                              
        for bead,xyz in zip(bead_types,coords0):                                                                                                                                                           
            gro.write('{:5d}{:5}{:>5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(1,mol_name,bead,i,xyz[0],xyz[1],xyz[2]))                                                                                         
            i += 1                                                                                                                                                                                         
        gro.write('5.0 5.0 5.0')         

