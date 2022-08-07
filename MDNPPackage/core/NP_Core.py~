"""
Author: Sang Young Noh
-----------------------
Last updated: 02/07/2022 
------------------------
"""

import re
import sys
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import math 
from operator import itemgetter
import itertools                                                                
import requests                                                                 
import collections                                                              
import random

# scipy libraries 
import scipy
from scipy.sparse import csr_matrix                                            
from scipy.sparse.csgraph import floyd_warshall                                
from scipy.spatial import ConvexHull, convex_hull_plot_2d 
from scipy.linalg import solve
from scipy.spatial import distance

# rdkit libraries 
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdchem          
from rdkit.Chem import rdMolDescriptors
from rdkit import RDConfig  

# alignment libraries in MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

# parmed functionality 
import parmed as pmd
from parmed.gromacs.gromacstop import GromacsTopologyFile
sys.path.append("..")

from MDNPPackage.connect.NP_Connect import NPConnect
from MDNPPackage.utils.NP_UTILS import generate_core

class CentralCoreGenerator(NPConnect):
    """
    """
    def __init__(self, R, points, gros, firstatoms, lastatoms, top1, top2, CG = 'CG',  option = 'Plain'):
        self.R = R
        self.points = points
        self.gros = gros
        self.firstatoms = firstatoms
        self.lastatoms = lastatoms
        self.top1 = GromacsTopologyFile(top1)
        self.top2 = GromacsTopologyFile(top2)
        self.option = option
        self.spherelist = generate_core(R, points, option)
        super().__init__(gros, firstatoms, lastatoms, self.spherelist, CG, option)
        
    def _fibanocci_sphere(self):
        """ Return a Fibanocci sphere with N number of points on the surface. 
        This will act as the template for the nanoparticle core. 
        """
        points = []
        phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

        for i in range(self.points):
            y = 1 - (i / float(self.points - 1)) * 2  # y goes from 1 to -1
            radius = math.sqrt(1 - y * y)  # radius at y
            theta = phi * i  # golden angle increment
            x = math.cos(theta) * radius
            z = math.sin(theta) * radius
            points.append((x, y, z))    

        return points

    def _rotation_matrix_from_vectors(vec1, vec2):
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

    def _label_np(self, core, nptype = 'Janus'):
        """
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
        xcoordinates = [i[0] for i in core] # Find x coordinates
        ycoordinates = [i[1] for i in core] # Find y coordinates
        zcoordinates = [i[2] for i in core] # Find z coordinates 
        length = 2 * abs(max(zcoordinates)) # From 2 * the radius, we know the total length of the NP 
    
        if nptype == 'Striped': 
            # As we have a spherical structure, we just need to find the minimum/maximum in 
            # one of the axes to find that for the rest 
            # define the threshold for how you wish to generate the NP with striped pattern 
            threshold = length / 3 
            # Find the central band of the sphere where you wish to put 
            # different ligands 
            stripedvalues = [i for i in core if i[2] > (min(zcoordinates) + threshold)
                             and i[2] < (max(zcoordinates) - threshold)]
            
            ceilingvalues = [i for i in core if i not in stripedvalues] 
            return [stripedvalues, ceilingvalues]
            
        elif nptype == 'Janus':
            # Same logic as with the striped example, but with the Janus pattern 
            threshold = length / 2 
            topvalues = [i for i in core if i[2] > (min(zCoordinates) + threshold)] 
            botvalues = [i for i in core if i not in topvalues] # Return bottom hemisphere 
            return [topvalues, botvalues]
           
    def _generate_core(self):
        """ Creates a Fibanocci sphere that represents the NP core 
        and allocates the radius. 

        The core is scaled down/up to the size that one wishes to have. 
        We can generate arrays corresponding  to a plain core, or a tuple with 
        two entries with different parts of the NP core that corresponds to positions 
        with striped or janus type positions.
        """
        spherelist = [] 
        sphere = self._fibanocci_sphere() # Create the fibanocci sphere representing the NP core 
        xsphere, ysphere, zsphere  = [], [], []

        for entry in sphere:
            xsphere.append(entry[0])
            ysphere.append(entry[1])
            zsphere.append(entry[2])
            
        # Append as 2d list
        for index in range(0, len(xsphere)):
            spherelist.append([xsphere[index], ysphere[index], zsphere[index]])
        # Take the radius value, and then multiply the unit vector in each 
        # Direction by that radius value to increase the total volume of the 
        # NP core.
        for index in range(0, len(spherelist)-1):
            spherelist[index][0] = spherelist[index][0] * self.R
            spherelist[index][1] = spherelist[index][1] * self.R
            spherelist[index][2] = spherelist[index][2] * self.R
        # Return just the whole list without any further modifications
        if self.option == 'Plain':
            return [spherelist[1:-1]]
        # Separate out the anisotropy for the Striped variant 
        elif self.option == 'Striped':
            stripedvalues, ceilingvalues = self._label_np(spherelist[1:-1], self.option)[0], self._label_np(spherelist[1:-1], self.option)[1]
            return stripedvalues, ceilingvalues
        # Separate out the anisotropy for the Janus variant 
        elif self.option == 'Janus':
            topvalues, bottomvalues = self._label_np(spherelist[1:-1], self.option)[0], self._label_np(spherelist[1:-1], self.option)[1]  
            return topvalues, bottomvalues

    def pandas_np(ligandstring, firstatom, lastatom, spherelist, 
                   ligandname, corename, length = 1.0):

        """Placeholder
        """
        transformationlist, namelist = [], [] # 
        ligandlist = [] 
        sphere = []
        xplot, yplot, zplot = [], [], []
        xplotsphere, yplotsphere, zplotsphere = [], [], []
        u = mda.Universe.from_smiles(LigandString)
        ligand = u.select_atoms('all')
        logging.info(f"The length of the ligand is {len(Ligand)}")
        firstatomgroup = u.select_atoms('name {}'.format(firstAtom))
        lastatomgroup = u.select_atoms('name {}'.format(lastAtom))
        ligandalignmentvector = (firstatomgroup.positions- lastatomgroup.positions)[0]
        
        for i,j in enumerate(ligand.positions):
            vector = (j - firstatomgroup.positions)[0]
            vector[0] = ligandalignmentvector[0] - vector[0]
            vector[1] = ligandalignmentvector[1] - vector[1]    
            vector[2] = ligandalignmentvector[2] - vector[2]
            if vector[0] == -math.inf:
                pass
            if vector[0] == 0.0:
                pass
            else:
                TransformationList.append([vector, Ligand.atoms[i].type])        
        
        vecligand = ligandalignmentvector.tolist()

        # Loop over the sphere and find the 
        for index in range(0, len(spherelist)):
            vec2 = spherelist[index]
            # Find the transformationvector for the ligand vector to vec2, which is the position of the point on sphere
            transformationvector = self._rotation_matrix_from_vectors(vecligand, vec2)
            # Rotate the vector 
            vec1rot = transformationvector.dot(vecligand) # Rotate the vector to match the surface point on the sphere 
            # Get the absolute length of the unit vector 
            unitvectorabs = np.linalg.norm(ligandalignmentvector)
            # Change the rotation vector in unit vector, then multiply by the absolute 
            # length of the sphere 
            vecmultiplier = vec1rot/unitvectorabs * (np.linalg.norm(np.array(vec2))) + (vec1rot/unitvectorabs * length)
            # Find the difference in length 
            sphere.append(vec2)
            # Translate the vector further out 
            for trans in transformationlist:
                ligandatomcoordinate = transformationvector.dot(trans[0])
                ligandatomcoordinate[0] = ligandatomcoordinate[0] + vecmultiplier[0]
                ligandatomcoordinate[1] = ligandatomcoordinate[1] + vecmultiplier[1]
                ligandatomcoordinate[2] = ligandatomcoordinate[2] + vecmultiplier[2]
                ligandlist.append(ligandatomcoordinate.tolist()) # Append coordinates of the 
                namelist.append(trans[1]) # Append the names of the atoms
        # Append the coordinates of the ligands 
        for index, entry in enumerate(ligandlist):
            xplot.append(entry[0])
            yplot.append(entry[1])
            zplot.append(entry[2])  
    
        ligandconstituent = [atom.name for atom in ligand]
        ligands = []
        for index in range(0, len(sphere)): 
            ligands = ligands + ligandconstituent
        spherename = [] 
        # Append the coordinates of the sphere 
        for entry in sphere:
            xplotsphere.append(entry[0])
            yplotsphere.append(entry[1])
            zplotsphere.append(entry[2])
            spherename.append('P5')
    
        dfligand = pd.DataFrame(list(zip(xplot, yplot, zplot, ligands)),
                                columns =['X', 'Y', 'Z', 'NAME'])
        
        dfcore = pd.DataFrame(list(zip(xplotsphere, yplotsphere, zplotsphere, spherename)),
                              columns =['X', 'Y', 'Z', 'NAME'])
        
        dfligand['RESNAME'] = ligandname
        dfcore['RESNAME'] = corename
        return dfligand, dfcore
    
    def pandas_np_martini(molecule, ligandalignmentvector, transformationlist, spherelist, ligandname, 
                              corename, length = 1.0):
        """ Function to read Martini molecule information and orientate on NP surface"""

        ligandList, namelist = [], []
        sphere = []
        xplot, yplot, zplot = [], [], []
        xplotsphere, yplotsphere, zplotsphere = [], [], []

        # Sulfur/ligand vector 

        vec1 = ligandalignmentvector.tolist()

        for index in range(0, len(spherelist)):
            vec2 = spherelist[index] 
            transformationvector = self._rotation_matrix_from_vectors(vec1, vec2)  
            vec1rot = transformationvector.dot(vec1) # Rotate the vector to match the surface point on the sphere 
            unitvectorabs = np.linalg.norm(liganalignmentvector)  
            vecmultiplier = vec1rot/unitvectorabs * (np.linalg.norm(np.array(vec2))) + (vec1_rot/unitvectorabs * length)
            sphere.append(vec2)
            
            # Get the factors to translate the vector 
            for trans in transformationlist:
                
                ligandatomcoordinate = transformationvector.dot(trans[0])
                ligandatomcoordinate[0] = ligandatomcoordinate[0] + vecmultiplier[0]
                ligandatomcoordinate[1] = ligandatomcoordinate[1] + vecmultiplier[1]
                ligandatomcoordinate[2] = ligandatomcoordinate[2] + vecmultiplier[2]
                ligandlist.append(ligandatomcoordinate.tolist()) 
                namelist.append(trans[1]) # Append the names of the atoms 

            # Append the coordinates of the ligands 
            for index, entry in enumerate(ligandlist):
                xplot.append(entry[0])
                yplot.append(entry[1])
                zplot.append(entry[2])
        
            # Add in the ligand index 
            ligandconstituent = [atom.name for atom in molecule] # Molecule is utilized here 
            for index in range(0, len(sphere)): 
                ligands = ligands + ligandconstituent
            spherename = [] 
            # Append the coordinates of the sphere 
            for entry in sphere:
                xplotsphere.append(entry[0])
                yplotsphere.append(entry[1])
                zplotsphere.append(entry[2])
                sphereName.append('P5')
            dfligand = pd.DataFrame(list(zip(xplot, yplot, zplot, ligands)),
                                    columns =['X', 'Y', 'Z', 'NAME'])
            dfligand['RESNAME'] = ligandname
            dfcore = pd.DataFrame(list(zip(xplotsphere, yplotsphere, zplotsphere, spherename)),
                                  columns =['X', 'Y', 'Z', 'NAME'])
            dfcore['RESNAME'] = corename
            return dfligand, dfcore
    
    def _core_network(self, corename, pandasfile = None, pickfile = None):
        """ Bond restraint allocator for the central core atoms of the NP 
        
        We add a very high spring constant constraint on the central atoms 
        to ensure that the central core can be considered as 'static'
        """
        npdataframe = self.return_ordered_coordinates()
        npcore = npdataframe[npdataframe['RESNAME'] == corename]
        npcore = npcore[['X', 'Y', 'Z']] 
        npcorearray = npcore.to_numpy() 
        duplicatearray = []
        returnstring = []
        
        for index, entry in enumerate(npcorearray):
            posgoal = np.array([entry])
            distmatrix = np.linalg.norm(npcorearray - posgoal, axis=1)
            for index2, entry in enumerate(distmatrix):
                if index == index2: # If we are looking at the same index, then pass
                    pass
                else:
                    if entry/10 >= 0.7: # If the length of the bonds is more than 0.7 nm, then pass - dont use that bond   
                        pass
                    # sorted out the indentation but the rest needs to be fixed 
                    elif index == 1:
                        entrystring = f"{index+1} {index2+1} 1 {entry/10} 5000"
                        sortedinput = [index, index2]
                        sortedinput.sort()
                        duplicatearray.append(sortedinput)
                        returnstring.append(entrystring)
                    elif index > 1:
                        sortedinput = [index, index2]
                        sortedinput.sort()
                        if sortedinput in duplicatearray:
                            pass
                        else:
                            duplicatearray.append(sortedinput)
                            entrystring = f"{index+1} {index2+1} 1 {entry/10} 5000"
                            returnstring.append(entrystring)           
        return returnstring

    def generate_np_itp(self, residuename = 'RES'):
        """ """
        atoms = []
        npdataframe = self.return_ordered_coordinates()
        # Create the initial bond network from the core_network function we have already created
        ligandbonds = self._core_network('Core')
        corelen = len(self.return_ordered_coordinates()[self.return_ordered_coordinates()['RESNAME'] == 'Core'])
        lig1len = len(self.return_ordered_coordinates()[self.return_ordered_coordinates()['RESNAME'] == 'Lig1']) 
        lig2len = len(self.return_ordered_coordinates()[self.return_ordered_coordinates()['RESNAME'] == 'Lig2'])

        indices1 = [i + corelen for i in self.return_ordered_coordinates()[self.return_ordered_coordinates()['RESNAME'] == 'Lig1']['index'].iloc[::len(self.top1.atoms)]]
        indices2 = [i + corelen + lig1len for i in self.return_ordered_coordinates()[self.return_ordered_coordinates()['RESNAME'] == 'Lig2']['index'].iloc[::len(self.top2.atoms)]]    
    
        atoms.append('[atoms]')
        atoms.append('; nr  type  resnr residue atom cgnr charge  mass')
        # Append core atom information 
    
        for index in range(0, corelen): 
            index = index + 1
            atomstring = f"{index} P5 1 {residuename} P5 {index} 0 100"
            atoms.append(atomstring)
        
        ligandstringimpropers = [] 
        ligandstringimpropers.append('[dihedrals]')
        ligandstringimpropers.append('; i j k l  funct  ref.angle   force_k')
        ligandstringimpropers.append('; Ligand 1 improper data')
        ligandbonds.append('; Ligand 1 bond data')
    
        for index in indices1: 
            index = index + 1
            # get bond parameters
            for bond in self.top1.bonds:
                bondstring = f"{bond.atom1.idx + (index)} {bond.atom2.idx + (index)} {bond.funct} {bond.type.req / 10} {bond.type.k}"
                ligandbonds.append(bondstring)
        
        # get improper dihedral parameters 
        for improper in self.top1.impropers:
            dihedralstring = f"{improper.atom1.idx + (index)} {improper.atom2.idx + (index)} {improper.atom3.idx + (index)} {improper.atom4.idx + (index)} {improper.funct} {improper.type.psi_eq} {improper.type.psi_k}"
            ligandstringimpropers.append(dihedralstring)
        
        # get atomic parameters
        for atom in self.top1.atoms:
            atomstring = f"{atom.idx + (index)} {atom.type} 1 {residuename} {atom.name} {atom.idx + (index)} {atom.charge} {atom.mass}"
            atoms.append(atomstring)
    
        ligandbonds.append('; Ligand 2 data')
        ligandstringimpropers.append('; Ligand 2 improper data')
    
        for index in indices2:
            index = index + 1 
            # ditto for the second lot 
            for bond in self.top2.bonds:
                bondstring = f"{bond.atom1.idx + (index)} {bond.atom2.idx + (index)} {bond.funct} {bond.type.req / 10} {bond.type.k}"
                ligandbonds.append(bondstring)
            for improper in self.top2.impropers:
                dihedralstring = f"{improper.atom1.idx + (index)} {improper.atom2.idx + (index)} {improper.atom3.idx + (index)} {improper.atom4.idx + (index)} {improper.funct} {improper.type.psi_eq} {improper.type.psi_k}"
                ligandstringimpropers.append(dihedralstring)            
            for atom in self.top2.atoms:
                atomstring = f"{atom.idx + (index)} {atom.type} 1 {residuename} {atom.name} {atom.idx + (index)} {atom.charge} {atom.mass}"
                atoms.append(atomstring)
    
        return ligandbonds, ligandstringimpropers, atoms

    def generate_coordinates(self, groname):
        """ Generate a gro file from the pandas dataframe 
            we generated from the previous functions 
        """
        coordinatesframe = self.return_ordered_coordinates()[['X','Y','Z']]# NPDataframe[['X','Y','Z']]
        coordinates = coordinatesframe.to_numpy()

        emptyuniverse = mda.Universe.empty(len(self.return_ordered_coordinates()), 
                                   1, atom_resindex=[0] * len(self.return_ordered_coordinates()), 
                                   trajectory=True)

        emptyuniverse.add_TopologyAttr('names')
        emptyuniverse.atoms.names = self.return_ordered_coordinates()['NAME'].to_numpy()
        emptyuniverse.load_new(coordinates)
        emptyuniverse.select_atoms('all').write(f'{groname}.gro')

    def generate_itp(self, itpname):
        """
        generate itp for the NP
        """
        bonds, improper, atoms = self.generate_np_itp()
        attachments = self.attach_ligands_martini()

        data = atoms
        data.extend(bonds)
        data.extend(attachments)
        data.extend(improper)

        with open(f'{itpname}.itp', 'w') as f:
            for item in data:
                f.write("%s\n" % item)
        
        
