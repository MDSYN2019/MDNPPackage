import math
import itertools
import requests
import collections
import random
from typing import Tuple
import re
import sys

import pandas as pd
import numpy as np
import plotly.graph_objs as go
from operator import itemgetter
import textwrap

# scipy libraries
import scipy
from scipy.spatial import distance

# alignment libraries in MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

# replacing parmed with vermouth..
import vermouth.forcefield
import vermouth.molecule
import vermouth.gmx.itp_read

# parmed functionality - may still need it for all-atomic functionality
import parmed as pmd
from parmed.gromacs.gromacstop import GromacsTopologyFile

sys.path.append("..")

from MDNPPackage.connect.NP_Connect import NPConnect
from MDNPPackage.utils.NP_Utils import generate_core


class CentralCoreGenerator(NPConnect):
    """
    the main generator class for the nanoparticle 
    """

    def __init__(
        self,
        R: float,
        points: list[float],
        gros: list[str],
        first_atoms: list[str],
        last_atoms: list[str],
        top_1: str,
        top_2: str,
        CG: str = "CG",
        option: str = "Plain",
    ):
        self.R = R
        self.points = points
        self.gros = gros
        self.first_atoms = first_atoms
        self.last_atoms = last_atoms
        self.top_1 = top_1
        self.top_2 = top_2
        self.dummy_ff_1 = vermouth.forcefield.ForceField(name="test_ff")
        self.dummy_ff_2 = vermouth.forcefield.ForceField(name="test_ff")
        self.option = option
        self.sphere_list = generate_core(R, points, option)
        super().__init__(gros, first_atoms, last_atoms, self.sphere_list, CG, option)

    def _rotation_matrix_from_vectors(vec_1: np.array, vec_2: np.array):
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
        a, b = (
            (vec_1 / np.linalg.norm(vec_1)).reshape(3),
            (vec_2 / np.linalg.norm(vec_2)).reshape(3),
        )
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        k_mat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + k_mat + k_mat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix

    def _label_np(self, core: list[list[float]], np_type: str = "Janus"):
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
        x_coordinates = [i[0] for i in core]  # Find x coordinates
        y_coordinates = [i[1] for i in core]  # Find y coordinates
        z_coordinates = [i[2] for i in core]  # Find z coordinates
        length = 2 * abs(
            max(z_coordinates)
        )  # From 2 * the radius, we know the total length of the NP

        if np_type == "Striped":
            # As we have a spherical structure, we just need to find the minimum/maximum in
            # one of the axes to find that for the rest
            # define the threshold for how you wish to generate the NP with striped pattern
            threshold = length / 3
            # Find the central band of the sphere where you wish to put
            # different ligands
            striped_values = [
                i
                for i in core
                if i[2] > (min(z_coordinates) + threshold)
                and i[2] < (max(z_coordinates) - threshold)
            ]

            ceiling_values = [i for i in core if i not in striped_values]
            return [striped_values, ceiling_values]

        elif np_type == "Janus":
            # Same logic as with the striped example, but with the Janus pattern
            threshold = length / 2
            top_values = [i for i in core if i[2] > (min(z_coordinates) + threshold)]
            bot_values = [i for i in core if i not in top_values]  # Return bottom hemisphere
            return [top_values, bot_values]

    def _core_network(self, core_name: str) -> list[str]:
        """ Bond restraint allocator for the central core atoms of the NP 
        
        We add a very high spring constant constraint on the central atoms 
        to ensure that the central core can be considered as 'static'
        """
        np_data_frame = self.return_ordered_coordinates()
        np_core = np_data_frame[np_data_frame["RESNAME"] == core_name]
        np_core = np_core[["X", "Y", "Z"]]
        np_core_array = np_core.to_numpy()
        duplicate_array = []
        return_string = []
        return_string.append("[ bonds ]")
        for index, entry in enumerate(np_core_array):
            pos_goal = np.array([entry])
            dist_matrix = np.linalg.norm(np_core_array - pos_goal, axis=1)
            for index_2, entry in enumerate(dist_matrix):
                if index == index_2:  # If we are looking at the same index, then pass
                    pass
                else:
                    if (
                        entry / 10 >= 0.7
                    ):  # If the length of the bonds is more than 0.7 nm, then pass - dont use that bond
                        pass
                    # sorted out the indentation but the rest needs to be fixed
                    elif index == 1:
                        entry_string = f"{index+1} {index_2+1} 1 {entry/10} 5000"
                        sorted_input = [index, index_2]
                        sorted_input.sort()
                        duplicate_array.append(sorted_input)
                        return_string.append(entry_string)
                    elif index > 1:
                        sorted_input = [index, index_2]
                        sorted_input.sort()
                        if sorted_input in duplicate_array:
                            pass
                        else:
                            duplicate_array.append(sorted_input)
                            entry_string = f"{index+1} {index_2+1} 1 {entry/10} 5000"
                            return_string.append(entry_string)
        return return_string

    def _generate_ligand_parameters(self):
        """
        use vermouth to get ligand information to add the molecular mechanics 
        forcefield information of the ligands when attached to the NP 
        """
        with open(self.top_1, "r") as d:
            data = d.read()
            top_1_lines = textwrap.dedent(data).splitlines()
        with open(self.top_2, "r") as d:
            data = d.read()
            top_2_lines = textwrap.dedent(data).splitlines()

        # fill itp_read with the ligand information
        vermouth.gmx.itp_read.read_itp(top_1_lines, self.dummy_ff_1)
        vermouth.gmx.itp_read.read_itp(top_2_lines, self.dummy_ff_2)

        self.ff_1_block = self.dummy_ff_1.blocks[list(self.dummy_ff_1.blocks.keys())[0]]
        self.ff_2_block = self.dummy_ff_2.blocks[list(self.dummy_ff_2.blocks.keys())[0]]
        self.ff_1_length = len(list(self.ff_1_block.find_atoms()))
        self.ff_2_length = len(list(self.ff_2_block.find_atoms()))
        self.ff_1_block_atoms = list(self.ff_1_block.atoms)
        self.ff_2_block_atoms = list(self.ff_2_block.atoms)

    def generate_np_itp(self, residue_name: str = "RES") -> tuple[list[str], list[str], list[str]]:
        """ 
        """
        atoms = []
        self._generate_ligand_parameters()
        np_data_frame = self.return_ordered_coordinates()
        # Create the initial bond network from the core_network function we have already created
        ligand_bonds = self._core_network("Core")
        core_len = len(
            self.return_ordered_coordinates()[
                self.return_ordered_coordinates()["RESNAME"] == "Core"
            ]
        )
        lig_1_len = len(
            self.return_ordered_coordinates()[
                self.return_ordered_coordinates()["RESNAME"] == "Lig1"
            ]
        )
        lig_2_len = len(
            self.return_ordered_coordinates()[
                self.return_ordered_coordinates()["RESNAME"] == "Lig2"
            ]
        )

        indices_1 = [
            i + core_len
            for i in self.return_ordered_coordinates()[
                self.return_ordered_coordinates()["RESNAME"] == "Lig1"
            ]["index"].iloc[
                :: self.ff_1_length
            ]  # len(self.top1.atoms)]
        ]
        indices_2 = [
            i + core_len + lig_1_len
            for i in self.return_ordered_coordinates()[
                self.return_ordered_coordinates()["RESNAME"] == "Lig2"
            ]["index"].iloc[
                :: self.ff_2_length
            ]  # len(self.top1.atoms)]
        ]

        atoms.append("[ atoms ]")
        atoms.append("; nr  type  resnr residue atom cgnr charge  mass")
        # Append core atom information

        for index in range(0, core_len):
            index = index + 1
            atom_string = f"{index} P5 1 {residue_name} P5 {index} 0 100"
            atoms.append(atom_string)

        ligand_string_impropers = []
        ligand_string_impropers.append("[ dihedrals ]")
        ligand_string_impropers.append("; i j k l  funct  ref.angle   force_k")
        ligand_string_impropers.append("; Ligand 1 improper data")
        ligand_bonds.append("; Ligand 1 bond data")

        for index in indices_1:
            index = index + 1
            # get bond parameters
            for information in self.ff_1_block.interactions["bonds"]:
                # for bond in self.top1.bonds:
                bond_string = f"{information[0][0] + (index)} {information[0][1] + (index)} {information[1][0]} {information[1][1]}"
                ligand_bonds.append(bond_string)

                # get improper dihedral parameters
            for dihedrals in self.ff_1_block.interactions["dihedrals"]:
                dihedral_string = f"{dihedrals[0][0] + (index)} {dihedrals[0][1] + (index)} {dihedrals[0][2] + (index)} {dihedrals[0][3] + (index)} {dihedrals[1][0]} {dihedrals[1][1]} {dihedrals[1][2]}"
                ligand_string_impropers.append(dihedral_string)

            for improper in self.ff_1_block.interactions["impropers"]:
                dihedral_string = f"{improper[0][0] + (index)} {improper[0][1] + (index)} {improper[0][2] + (index)} {improper[0][3] + (index)} {improper.funct} {improper.type.psi_eq} {improper.type.psi_k}"
                ligand_string_impropers.append(dihedral_string)
            # get atomic parameters
            # for atom in self.top1.atoms:
            for atom in self.ff_1_block_atoms:
                atomic_index = atom["index"] - 1
                atom_string = f"{atomic_index + (index)} {atom['atype']} 1 {residue_name} {atom['atomname']} {atomic_index + (index)} {atom['charge']} {0}"
                atoms.append(atom_string)

        ligand_bonds.append("; Ligand 2 bond data")
        ligand_string_impropers.append("; Ligand 2 improper data")

        for index in indices_2:
            index = index + 1
            # ditto for the second lot
            # for bond in self.top2.bonds:
            for information in self.ff_2_block.interactions["bonds"]:
                # bond_string = f"{bond.atom1.idx + (index)} {bond.atom2.idx + (index)} {bond.funct} {bond.type.req / 10} {bond.type.k}"
                bond_string = f"{information[0][0] + (index)} {information[0][1] + (index)} {information[1][0]} {information[1][1]}"
                ligand_bonds.append(bond_string)

            for dihedrals in self.ff_2_block.interactions["dihedrals"]:
                dihedral_string = f"{dihedrals[0][0] + (index)} {dihedrals[0][1] + (index)} {dihedrals[0][2] + (index)} {dihedrals[0][3] + (index)} {dihedrals[1][0]} {dihedrals[1][1]} {dihedrals[1][2]}"
                ligand_string_impropers.append(dihedral_string)

            # for atom in self.top2.atoms:
            for atom in self.ff_2_block_atoms:
                atomic_index = atom["index"] - 1
                atom_string = f"{atomic_index + (index)} {atom['atype']} 1 {residue_name} {atom['atomname']} {atomic_index + (index)} {atom['charge']} {0}"
                atoms.append(atom_string)

        return ligand_bonds, ligand_string_impropers, atoms

    def _tile_universe(self, universe: mda.Universe, n_x: int, n_y: int, n_z: int):
        box = universe.dimensions[:3]
        copied = []
        for x in range(n_x):
            for y in range(n_y):
                for z in range(n_z):
                    u_ = universe.copy()
                    move_by = box * (x, y, z)
                    u_.atoms.translate(move_by)
                    copied.append(u_.atoms)

        new_universe = mda.Merge(*copied)
        new_box = box * (n_x, n_y, n_z)
        new_universe.dimensions = list(new_box) + [90] * 3
        return new_universe

    def generate_coordinates(self, gro_name: str):
        """ Generate a gro file from the pandas dataframe 
            we generated from the previous functions 
        """
        coordinates_frame = self.return_ordered_coordinates()[
            ["X", "Y", "Z"]
        ]  # NPDataframe[['X','Y','Z']]
        coordinates = coordinates_frame.to_numpy()

        empty_universe = mda.Universe.empty(
            len(self.return_ordered_coordinates()),
            1,
            atom_resindex=[0] * len(self.return_ordered_coordinates()),
            trajectory=True,
        )
        empty_universe.add_TopologyAttr("names")
        empty_universe.atoms.names = self.return_ordered_coordinates()["NAME"].to_numpy()
        empty_universe.load_new(coordinates)
        empty_universe.dimensions = [73, 73, 73, 90, 90, 90]
        empty_universe.select_atoms("all").write(f"{gro_name}.gro")

    def generate_itp(self, itp_name: str):
        """
        generate itp for the NP
        """
        bonds, improper, atoms = self.generate_np_itp()
        attachments = self.attach_ligands_martini()
        data = atoms
        data.extend(bonds)
        data.extend(attachments)
        data.extend(improper)
        with open(f"{itp_name}.itp", "w") as f:
            for item in data:
                f.write("%s\n" % item)
