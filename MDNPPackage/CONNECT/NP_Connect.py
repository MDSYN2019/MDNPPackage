"""
Author: Sang Young Noh
----------------------

Last Updated: 31/05/2021
------------------------

Code to generate the bond restraints where PyCGtools do not seem to be adequate. Not a 
cricism of PyCGTools, but that it is perhaps not designed to optimize NP structures 

For the CA CA bonds

CA CA         1    0.14000   392459.2 ; 7,(1986),230; BENZENE,PHE,TRP,TYR ; This is the amber aa - between aromatic AA
CA S          1    0.17500   189953.6 ; Au_cluster_ff ; This is the parameters for the amber aa - between aromatic AA and Sulfer

"""
## Module imports

import numpy as np 
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array

universe = MDAnalysis.Universe("out.gro") ## Reading in the .gro file5
# Ligand attached P5 atoms
P5Atoms = universe.select_atoms("name P5")
P5ID = P5Atoms.atoms.ids
# ST atoms 
STAtoms = universe.select_atoms("name SB0 PB0")
STID = STAtoms.atoms.ids
# BEN ligands - Aromatic ligands 
PETAtoms = universe.select_atoms("resname PET")
PETID = PETAtoms.atoms.ids
# BEN ligands - Polar ligands
SAtoms = universe.select_atoms("resname S")
SID = SAtoms.atoms.ids

P5AtomsPositionArray = P5Atoms.positions # Core atoms
STAtomsPositionArray = STAtoms.positions # Sulfur atoms
PETAtomsPositionArray = PETAtoms.positions # Hydrophobic ligand atoms 
SAtomsPositionArray = SAtoms.positions # Polar ligands atoms 

# Print out whole index
# We follow the format of the NP_template of the previous MARTINI NP I made
# which has the following format:

"""
How to print out the atoms part of the itp file5

[ moleculetype ]
NP             1

[ atoms ]
; nr  type  resnr residue atom cgnr charge  mass
     1       P5        1   NP   P5        1    0.0000    30.9738
     2       C1        1   NP   C1        2    0.0000    12.0110
     3       C2        1   NP   C2        3    0.0000    12.0110
     4       Qa        1   NP   Qa        4    0.0000    12.0000
     5       P5        1   NP   P5        5    0.0000    30.9738
     6       C1        1   NP   C1        6    0.0000    12.0110
     7       C2        1   NP   C2        7    0.0000    12.0110
"""

"""
Commentary on the masses of the beads 

In the case of 

"""

#BeadDict = {'P5': 'TC6', 'SB0' : 'TC6', 'SB1' : 'P5', 'SB2' : 'SC4', 'SB3' : 'TP1', 'SB4' : 'TP1', 'SB5': 'TP1', 'SB6' : 'TC3', 'PB0' : 'TC6', 'PB1' : 'TC2', 'PB2' : 'TC2', 'PB3' : 'TC5', 'PB4' : 'TC5', 'PB5' : 'TC5' } # Dict to store name and type of the MARTINI beads used for making

#BeadMass = {'P5' : 33.0, 'SB0' : 33.0, 'SB1' : 43.0, 'SB2' : 25.0, 'SB3' : 28.0, 'SB4' : 28.0, 'SB5' : 28.0, 'SB6': 14.0 , 'PB0' : 30.0, 'PB1' : 14.0, 'PB2': 14, 'PB3' : 15, 'PB4': 26, 'PB5': 26} # Dict to store name and mass of the MARTINI beads used for making the CG NP 

BeadDict = {'P5': 'C6', 'SB0' : 'C6', 'SB1' : 'P5', 'SB2' : 'SC3', 'SB3' : 'TP1', 'SB4' : 'TP1', 'SB5': 'TP1', 'SB6' : 'TC3', 'PB0' : 'C6', 'PB1' : 'SC4', 'PB2' : 'TC2', 'PB3' : 'TC5', 'PB4' : 'TC5', 'PB5' : 'TC5' } # Dict to store name and type of the MARTINI beads used for making
BeadMass = {'P5' : 196.96, 'SB0' : 33.0, 'SB1' : 43.0, 'SB2' : 25.0, 'SB3' : 28.0, 'SB4' : 28.0, 'SB5' : 28.0, 'SB6': 14.0 , 'PB0' : 30.0, 'PB1' : 28.0, 'PB2': 14, 'PB3' : 15, 'PB4': 26, 'PB5': 26} # Dict to store name and mass of the MARTINI beads used for making the CG NP


print ('[ moleculetype ]')
print ('NP	       1')
allAtoms = universe.select_atoms('all')
print ("[ atoms ]")
print ("; nr  type  resnr residue atom cgnr charge  mass")
for row in allAtoms:
	#print (row.id)
	print ("     {}       {}        {}   {}   {}        {}    {}    {}".format(row.id, BeadDict[row.name], 1, 'NP', row.name, row.id, 0.000, BeadMass[row.name]))
	
"""
 We want to print out bond restraints in the following format:

index, index2, bond type (?), distance, spring constant. For example: 

[ bonds ]
1 785  1 0.20        5000			
1 397  1 0.20        5000

# ----------------------------------------------------------------
##  As for angular restraints, we want them to be printed out as:
# -------------------------------------------------------- #
# index, index2, index3, bond type (?), angle, spring constant  #
# -------------------------------------------------------- #
#[ angles ]
#2 3 4 2 180 25
#6 7 8 2 180 25
#10 11 12 2 180 25
#14 15 16 2 180 25

"""
# Some sanity checks - Check that the number of ids selected actually match the number of atoms within the array 
assert (len(P5ID) == len(P5AtomsPositionArray))
assert (len(STID) == len(STAtomsPositionArray))
assert (len(PETID) == len(PETAtomsPositionArray))

#print (STID, P5ID)

## Find the closest ST group to the Gold Core  
print ('[ bonds ]')
print('; i  j  func')
print ('; P5 - P5 ')
DuplicateArray = []
for index, atom in enumerate(P5AtomsPositionArray): # Loop over the P5 position array                                                                                                                      
        distarray = (distance_array(atom, P5AtomsPositionArray)) # Find the distance between the P5 atoms and the selected index                                                                            
        for index2, entry in enumerate(distarray[0]):
                if index == index2: # If we are looking at the same index (same P5 atom), then pass                                                                                                         
                       pass
                else:
                        if distarray[0][index2]/10 >= 0.7: # If the length of the bonds is more than 0.7, then pass - dont use that bond                                                                   
	                        pass
                        elif index == 1:
                                print (P5ID[index], P5ID[index2], 1, distarray[0][index2]/10, 5000) # Else, allocate a strong contraint between the P5 files to hold the core together                     
                                sortedInput = [P5ID[index], P5ID[index2]]
                                sortedInput.sort()
                                DuplicateArray.append(sortedInput)
                        elif index > 1:
                                sortedInput = [P5ID[index], P5ID[index2]]
                                sortedInput.sort()
                                if sortedInput in DuplicateArray: # if already contained in array, then dont do anything                                                                                   
                                        pass
                                else:
                                        print (P5ID[index], P5ID[index2], 1, distarray[0][index2]/10, 5000) # The ones printed out here has gone through the check of the duplicates so print              
                                        DuplicateArray.append(sortedInput) # Store the combination to make sure it is checked for duplicates as it goes down                                                
#for index, atom in enumerate(P5AtomsPositionArray): # Loop over the P5 position arrays
#	distarray = (distance_array(atom, P5AtomsPositionArray)) # Find the distance between the P5 atoms and the selected index 
#	for index2, entry in enumerate(distarray[0]):
#		if index == index2: # If we are looking at the same index (same P5 atom), then pass. 
#			pass
#		else:
#			print (P5ID[index], P5ID[index2], 1, distarray[0][index2]/10, 1000) # Else, allocate a strong contraint between the P5 files to hold the core together 

# Now to identify the gold core atoms that are the closest to the sulfer atoms  
print ('; ST - P5 ')  
for index, atom in enumerate(STAtomsPositionArray): # Loop over the sulfer atoms                     
        distarray = (distance_array(atom, P5AtomsPositionArray)) # Find the shortest distance between the ST and the appropriate gold core and allocate the print                                         
        for NPindex, entry in enumerate(distarray[0]):
                if entry/10 <= 0.40:
                        print (STID[index], P5ID[NPindex], 1, entry/10, 5000)


# Allocate the bond parameters for the hydrophobic ligands 
print ('; PET ligands ')
"""
<AtomGroup [<Atom 181: PB0 of type P of resname PET, resid 85 and segid SYSTEM>, <Atom 182: PB1 of type P of resname PET, resid 85 and segid SYSTEM>, <Atom 183: PB2 of type P of resname PET, resid 85 and segid SYSTEM>, <Atom 184: PB3 of type P of resname PET, resid 85 and segid SYSTEM>, <Atom 185: PB4 of type P of resname PET, resid 85 and segid SYSTEM>, <Atom 186: PB5 of type P of resname PET, resid 85 and segid SYSTEM>]>

[PET]
PB0 PB1
PB1 PB2
PB2 PB3
PB3 PB4
PB4 PB5
PB5 PB3
"""
for res in PETAtoms.residues:
        # First index is ST, second SC1, third SC2 and fourth SC4
        ligandResIDS = res.atoms.ids 
        #	print (ligandsResIDS)
        print (ligandResIDS[0], ligandResIDS[1], 1, 0.216, 5000) #  ST - C1  
        print (ligandResIDS[1], ligandResIDS[2], 1, 0.225, 5000) #  C1 - C2
#        print (ligandResIDS[2], ligandResIDS[3], 1, 0.2179, 5000.2) #  C2 - C3 
#        print (ligandResIDS[3], ligandResIDS[4], 1, 0.2220, 5000.2)
#        print (ligandResIDS[4], ligandResIDS[2], 1, 0.2220, 5000.2)
        
# Allocate the bond parameters for the polar ligands 
print ('; S ligands ')
"""
<AtomGroup [<Atom 69: SB0 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 70: SB1 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 71: SB2 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 72: SB3 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 73: SB4 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 74: SB5 of type S of resname S, resid 69 and segid SYSTEM>, <Atom 75: SB6 of type S of resname S, resid 69 and segid SYSTEM>]>

[S]
SB0 SB1 
SB1 SB6  
SB1 SB3 
SB2 SB3
SB3 SB4
SB4 SB5
SB5 SB6

"""
for res in SAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
        ligandResIDS = res.atoms.ids
        #	print (ligandsResIDS)
        print (ligandResIDS[0], ligandResIDS[1], 1, 0.2400, 5000) # ST - RA 
        print (ligandResIDS[1], ligandResIDS[2], 1, 0.2300, 5000) # SB2 -SB3
        print (ligandResIDS[2], ligandResIDS[3], 1, 0.3040, 5000) # SB5 - SB6
        print (ligandResIDS[2], ligandResIDS[4], 1, 0.2710, 5000) # SB4 - SB5 
        print (ligandResIDS[2], ligandResIDS[5], 1, 0.2850, 5000) # SB5 - SB1
        #print (ligandResIDS[1], ligandResIDS[3], 1, 0.3000, 5000) # SB2 - SB6 
        #print (ligandResIDS[2], ligandResIDS[4], 1, 0.2580, 5000) # SB3 - SB4 
        #print (ligandResIDS[4], ligandResIDS[5], 1, 0.2580, 5000) # SB4 - SB5 
        #print (ligandResIDS[5], ligandResIDS[6], 1, 0.2140, 5000) # SB5 - SB1
        #print (ligandResIDS[6], ligandResIDS[3], 1, 0.2310, 5000) # SB5 - SB6

print ('#ifndef FLEXIBLE')
print ('[ constraints ]')
print ('#endif')
for res in SAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
        ligandResIDS = res.atoms.ids 
        #	print (ligandsResIDS)
        #print (ligandResIDS[0], ligandResIDS[1], 1, 0.2400, 5000.2) # ST - RA 
        #print (ligandResIDS[1], ligandResIDS[2], 1, 0.2160, 5000.2) # SB2 -SB3 
        #print (ligandResIDS[1], ligandResIDS[3], 1, 0.3000, 1000000.2) # SB2 - SB6 
        #print (ligandResIDS[2], ligandResIDS[4], 1, 0.2510, 1000000.2) # SB3 - SB4 
        #print (ligandResIDS[4], ligandResIDS[5], 1, 0.2580, 1000000.2) # SB4 - SB5 
        #print (ligandResIDS[5], ligandResIDS[6], 1, 0.2140, 1000000.2) # SB5 - SB1
        #print (ligandResIDS[2], ligandResIDS[3], 1, 0.3040, 1000000.2) # SB5 - SB6
        #print (ligandResIDS[2], ligandResIDS[4], 1, 0.2710, 1000000.2) # SB4 - SB5 
        #print (ligandResIDS[2], ligandResIDS[5], 1, 0.2850, 1000000.2) # SB5 - SB1
for res in PETAtoms.residues:
        # First index is ST, second SC1, third SC2 and fourth SC4
        ligandResIDS = res.atoms.ids 
        print (ligandResIDS[2], ligandResIDS[3], 1, 0.2179, 1000000.2) #  C2 - C3 
        print (ligandResIDS[3], ligandResIDS[4], 1, 0.2220, 1000000.2)
        print (ligandResIDS[4], ligandResIDS[2], 1, 0.2220, 1000000.2)

print ('[ angles ]')
for res in SAtoms.residues:
#	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	#
	ligandResIDS = res.atoms.ids 	
	print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], 1, 165, 200) # ST - RA - R1  
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], 1, 115, 200) # R1 - R2 - R3
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[4], 1, 146, 200) # R2 - R3 - R1
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[5], 1, 160, 200) # R2 - R3 - R1
#	ligandRes = res.atoms

"""
[angle1]
  1  2  3
[angle2]
  2  3  4
[angle3]
  3  4  5
[angle4]
  4  5  6
[angle5]
  6  4  3
"""
for res in PETAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
	ligandResIDS = res.atoms.ids
	print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], 1, 166, 200) # ST - C1 - C2 
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], 1, 60, 200) # C1 - C2 - C3 
	print (ligandResIDS[3], ligandResIDS[4], ligandResIDS[2], 1, 58, 200) # C2 -  - R3
	ligandRes = res.atoms
        
#print ('; angles ')
#for index, atom in enumerate(STAtomsPositionArray):
#	distarray = (distance_array(atom, P5AtomsPositionArray))
#	print (P5ID[np.argmin(distarray)], STID[index], STID[index]+1, 1, 120, 25) # Explanation required 		

print('[ dihedrals ]')
for res in PETAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
	ligandResIDS = res.atoms.ids
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4], 2, 0, 50) # ST - RA - R1  
#for res in SAtoms.residues:
#	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
#	ligandResIDS = res.atoms.ids 	
#	print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], ligandResIDS[4], 2, 0, 50) # ST - RA - R1  

print('[ exclusions ]')
for res in PETAtoms.residues:
        ligandResIDS = res.atoms.ids
        print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4]) # ST - RA - R1  
        print (ligandResIDS[1], ligandResIDS[0], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4]) # ST - RA - R1  
        print (ligandResIDS[2], ligandResIDS[0], ligandResIDS[1], ligandResIDS[3], ligandResIDS[4]) # ST - RA - R1  
        print (ligandResIDS[3], ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], ligandResIDS[4]) # ST - RA - R1  
        print (ligandResIDS[4], ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], ligandResIDS[3]) # ST - RA - R1
        
for res in SAtoms.residues:
        ligandResIDS = res.atoms.ids 	
        print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4], ligandResIDS[5]) # ST - RA - R1  
        print (ligandResIDS[1], ligandResIDS[0], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4], ligandResIDS[5]) # ST - RA - R1  
        print (ligandResIDS[2], ligandResIDS[1], ligandResIDS[0], ligandResIDS[3], ligandResIDS[4], ligandResIDS[5]) # ST - RA - R1  
        print (ligandResIDS[3], ligandResIDS[1], ligandResIDS[2], ligandResIDS[0], ligandResIDS[4], ligandResIDS[5]) # ST - RA - R1  
        print (ligandResIDS[4], ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], ligandResIDS[0], ligandResIDS[5]) # ST - RA - R1  
        print (ligandResIDS[5], ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], ligandResIDS[4], ligandResIDS[0]) # ST - RA - R1  
        
for index, atom in enumerate(P5AtomsPositionArray): # Loop over the P5 position arrays
        P5IDTemp = np.delete(P5ID, index)
        print(P5ID[index], *P5IDTemp)
