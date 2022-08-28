from typing import Tuple
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
import textwrap

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

# replacing parmed with vermouth..
import vermouth.forcefield
import vermouth.molecule
import vermouth.gmx.itp_read

# parmed functionality
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

    def _fibanocci_sphere(self) -> list[list[float]]:
        """ Return a Fibanocci sphere with N number of points on the surface. 
        This will act as the template for the nanoparticle core. 
        """
        points = []
        phi = math.pi * (3.0 - math.sqrt(5.0))  # golden angle in radians

        for i in range(self.points):
            y = 1 - (i / float(self.points - 1)) * 2  # y goes from 1 to -1
            radius = math.sqrt(1 - y * y)  # radius at y
            theta = phi * i  # golden angle increment
            x = math.cos(theta) * radius
            z = math.sin(theta) * radius
            points.append((x, y, z))

        return points

    def _rotation_matrix_from_vectors(vec_1, vec_2):
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

    def _label_np(self, core, np_type="Janus"):
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

    def _generate_core(self) -> list[list[float]]:
        """ Creates a Fibanocci sphere that represents the NP core 
        and allocates the radius. 

        The core is scaled down/up to the size that one wishes to have. 
        We can generate arrays corresponding  to a plain core, or a tuple with 
        two entries with different parts of the NP core that corresponds to positions 
        with striped or janus type positions.
        """
        sphere_list = []
        sphere = self._fibanocci_sphere()  # Create the fibanocci sphere representing the NP core
        x_sphere, y_sphere, z_sphere = [], [], []

        for entry in sphere:
            x_sphere.append(entry[0])
            y_sphere.append(entry[1])
            z_sphere.append(entry[2])
        # Append as 2d list
        for index in range(0, len(x_sphere)):
            sphere_list.append([x_sphere[index], y_sphere[index], z_sphere[index]])
        # Take the radius value, and then multiply the unit vector in each
        # Direction by that radius value to increase the total volume of the
        # NP core.
        for index in range(0, len(sphere_list) - 1):
            sphere_list[index][0] = sphere_list[index][0] * self.R
            sphere_list[index][1] = sphere_list[index][1] * self.R
            sphere_list[index][2] = sphere_list[index][2] * self.R
        # Return just the whole list without any further modifications
        if self.option == "Plain":
            return [sphere_list[1:-1]]
        # Separate out the anisotropy for the Striped variant
        elif self.option == "Striped":
            striped_values, ceiling_values = (
                self._label_np(sphere_list[1:-1], self.option)[0],
                self._label_np(sphere_list[1:-1], self.option)[1],
            )
            return striped_values, ceiling_values
        # Separate out the anisotropy for the Janus variant
        elif self.option == "Janus":
            top_values, bottom_values = (
                self._label_np(sphere_list[1:-1], self.option)[0],
                self._label_np(sphere_list[1:-1], self.option)[1],
            )
            return top_values, bottom_values

    def pandas_np(
        ligand_string, first_atom, last_atom, sphere_list, ligand_name, core_name, length=1.0
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Placeholder
        """
        transformation_list, name_list = [], []  #
        ligand_list = []
        sphere = []
        x_plot, y_plot, z_plot = [], [], []
        x_plot_sphere, y_plot_sphere, z_plot_sphere = [], [], []
        u = mda.Universe.from_smiles(ligand_string)
        ligand = u.select_atoms("all")
        logging.info(f"The length of the ligand is {len(Ligand)}")
        first_atom_group = u.select_atoms("name {}".format(first_atom))
        last_atom_group = u.select_atoms("name {}".format(last_atom))
        ligand_alignment_vector = (first_atom_group.positions - last_atom_group.positions)[0]
        for i, j in enumerate(ligand.positions):
            vector = (j - first_atom_group.positions)[0]
            vector[0] = ligand_alignment_vector[0] - vector[0]
            vector[1] = ligand_alignment_vector[1] - vector[1]
            vector[2] = ligand_alignment_vector[2] - vector[2]
            if vector[0] == -math.inf:
                pass
            if vector[0] == 0.0:
                pass
            else:
                transformation_list.append([vector, ligand.atoms[i].type])
        vec_ligand = ligand_alignment_vector.tolist()
        # Loop over the sphere and find the
        for index in range(0, len(sphere_list)):
            vec_2 = sphere_list[index]
            # Find the transformationvector for the ligand vector to vec2, which is the position of the point on sphere
            transformation_vector = self._rotation_matrix_from_vectors(vec_ligand, vec_2)
            # Rotate the vector
            vec_1_rot = transformation_vector.dot(
                vec_ligand
            )  # Rotate the vector to match the surface point on the sphere
            # Get the absolute length of the unit vector
            unit_vector_abs = np.linalg.norm(ligand_alignment_vector)
            # Change the rotation vector in unit vector, then multiply by the absolute
            # length of the sphere
            vec_multiplier = vec_1_rot / unit_vector_abs * (np.linalg.norm(np.array(vec_2))) + (
                vec_1_rot / unit_vector_abs * length
            )
            # Find the difference in length
            sphere.append(vec_2)
            # Translate the vector further out
            for trans in transformation_list:
                ligand_atom_coordinate = transformation_vector.dot(trans[0])
                ligand_atom_coordinate[0] = ligand_atom_coordinate[0] + vec_multiplier[0]
                ligand_atom_coordinate[1] = ligand_atom_coordinate[1] + vec_multiplier[1]
                ligand_atom_coordinate[2] = ligand_atom_coordinate[2] + vec_multiplier[2]
                ligand_list.append(ligand_atom_coordinate.tolist())  # Append coordinates of the
                name_list.append(trans[1])  # Append the names of the atoms
        # Append the coordinates of the ligands
        for index, entry in enumerate(ligand_list):
            x_plot.append(entry[0])
            y_plot.append(entry[1])
            z_plot.append(entry[2])

        ligand_constituent = [atom.name for atom in ligand]
        ligands = []
        for index in range(0, len(sphere)):
            ligands = ligands + ligand_constituent
        sphere_name = []
        # Append the coordinates of the sphere
        for entry in sphere:
            x_plot_sphere.append(entry[0])
            y_plot_sphere.append(entry[1])
            z_plot_sphere.append(entry[2])
            sphere_name.append("P5")

        df_ligand = pd.DataFrame(
            list(zip(x_plot, y_plot, z_plot, ligands)), columns=["X", "Y", "Z", "NAME"]
        )
        df_core = pd.DataFrame(
            list(zip(x_plot_sphere, y_plot_sphere, z_plot_sphere, sphere_name)),
            columns=["X", "Y", "Z", "NAME"],
        )
        df_ligand["RESNAME"] = ligand_name
        df_core["RESNAME"] = core_name
        return df_ligand, df_core

    def pandas_np_martini(
        molecule,
        ligand_alignment_vector,
        transformation_list,
        sphere_list,
        ligand_name,
        core_name,
        length=1.0,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """ Function to read Martini molecule information and orientate on NP surface"""

        ligand_list, name_list = [], []
        sphere = []
        x_plot, y_plot, z_plot = [], [], []
        x_plot_sphere, y_plot_sphere, z_plot_sphere = [], [], []
        # Sulfur/ligand vector
        vec_1 = ligand_alignment_vector.tolist()
        for index in range(0, len(sphere_list)):
            vec_2 = sphere_list[index]
            transformation_vector = self._rotation_matrix_from_vectors(vec1, vec2)
            vec_1_rot = transformation_vector.dot(
                vec1
            )  # Rotate the vector to match the surface point on the sphere
            unit_vector_abs = np.linalg.norm(ligand_alignment_vector)
            vec_multiplier = vec_1_rot / unit_vector_abs * (np.linalg.norm(np.array(vec_2))) + (
                vec_1_rot / unit_vector_abs * length
            )
            sphere.append(vec_2)
            # Get the factors to translate the vector
            for trans in transformation_list:
                ligand_atom_coordinate = transformation_vector.dot(trans[0])
                ligand_atom_coordinate[0] = ligand_atom_coordinate[0] + vec_multiplier[0]
                ligand_atom_coordinate[1] = ligand_atom_coordinate[1] + vec_multiplier[1]
                ligand_atom_coordinate[2] = ligand_atom_coordinate[2] + vec_multiplier[2]
                ligand_list.append(ligand_atom_coordinate.tolist())
                name_list.append(trans[1])  # Append the names of the atoms

            # Append the coordinates of the ligands
            for index, entry in enumerate(ligand_list):
                x_plot.append(entry[0])
                y_plot.append(entry[1])
                z_plot.append(entry[2])

            # Add in the ligand index
            ligand_constituent = [atom.name for atom in molecule]  # Molecule is utilized here
            for index in range(0, len(sphere)):
                ligands = ligands + ligand_constituent

            sphere_name = []
            # Append the coordinates of the sphere
            for entry in sphere:
                x_plot_sphere.append(entry[0])
                y_plot_sphere.append(entry[1])
                z_plot_sphere.append(entry[2])
                sphere_name.append("P5")

            df_ligand = pd.DataFrame(
                list(zip(x_plot, y_plot, z_plot, ligands)), columns=["X", "Y", "Z", "NAME"]
            )
            df_ligand["RESNAME"] = ligand_name
            df_core = pd.DataFrame(
                list(zip(x_plot_sphere, y_plot_sphere, z_plot_sphere, sphere_name)),
                columns=["X", "Y", "Z", "NAME"],
            )
            dfcore["RESNAME"] = core_name
            return df_ligand, df_core

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
        self.ff_1_block_atoms = self.ff_1_block.atoms
        self.ff_2_block_atoms = self.ff_2_block.atoms

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

        atoms.append("[atoms]")
        atoms.append("; nr  type  resnr residue atom cgnr charge  mass")
        # Append core atom information

        for index in range(0, core_len):
            index = index + 1
            atom_string = f"{index} P5 1 {residue_name} P5 {index} 0 100"
            atoms.append(atom_string)

        ligand_string_impropers = []
        ligand_string_impropers.append("[dihedrals]")
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

        ligand_bonds.append("; Ligand 2 data")
        ligand_string_impropers.append("; Ligand 2 improper data")

        for index in indices_2:
            index = index + 1
            # ditto for the second lot
            # for bond in self.top2.bonds:
            for information in self.ff_2_block.interactions["bonds"]:
                # bond_string = f"{bond.atom1.idx + (index)} {bond.atom2.idx + (index)} {bond.funct} {bond.type.req / 10} {bond.type.k}"
                bond_string = f"{information[0][0] + (index)} {information[0][1] + (index)} {information[1][0]} {information[1][1]} {bond.type.k}"
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

    def _tile_universe(self, universe, n_x, n_y, n_z):
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
