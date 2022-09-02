from typing import Tuple
import MDAnalysis as mda
import numpy as np
import pandas as pd
import math


def label_np(core: list[list[float]], np_type: str = "Janus") -> list[list[float]]:
    """
    Depending on the type of NP we want in the input, we can try to generate 
    different patterns on the surface of the spehre, which will help us generate the 
    lists correponding to the anisotropic nature of the NP. 
    """
    x_coordinates = [i[0] for i in core]  # Find x coordinates
    y_coordinates = [i[1] for i in core]  # Find y coordinates
    z_coordinates = [i[2] for i in core]  # Find z coordinates
    length = 2 * abs(max(z_coordinates))  # From 2 * the radius, we know the total length of the NP

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
            if i[2] > (min(z_coordinates) + threshold) and i[2] < (max(z_coordinates) - threshold)
        ]

        ceiling_values = [i for i in core if i not in striped_values]
        return [striped_values, ceiling_values]

    elif np_type == "Janus":
        # Same logic as with the striped example, but with the Janus pattern
        threshold = length / 2
        top_values = [i for i in core if i[2] > (min(z_coordinates) + threshold)]
        bot_values = [i for i in core if i not in top_values]  # Return bottom hemisphere
        return [top_values, bot_values]


def fibanocci_sphere(sample_points : int) -> list[Tuple[float]]:
    """ Return a Fibanocci sphere with N number of points on the surface. 
        This will act as the template for the nanoparticle core. 
    """
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))  # golden angle in radians

    for i in range(sample_points):
        y = 1 - (i / float(sample_points - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append((x, y, z))
    return points

def generate_core(radius: float, n : int , option="Plain") -> Tuple[list[list[float]], list[list[float]]]:
    """ Creates a Fibanocci sphere that represents the NP core 
        and allocates the radius. 

        The core is scaled down/up to the size that one wishes to have. 
        We can generate arrays corresponding  to a plain core, or a tuple with 
        two entries with different parts of the NP core that corresponds to positions 
        with striped or janus type positions.
    """
    spherelist = []
    sphere = fibanocci_sphere(n)  # Create the fibanocci sphere representing the NP core
    xsphere, ysphere, zsphere = [], [], []

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
    for index in range(0, len(spherelist) - 1):
        spherelist[index][0] = spherelist[index][0] * radius
        spherelist[index][1] = spherelist[index][1] * radius
        spherelist[index][2] = spherelist[index][2] * radius
    # Return just the whole list without any further modifications
    if option == "Plain":
        return [spherelist[1:-1]]
    # Separate out the anisotropy for the Striped variant
    elif option == "Striped":
        stripedvalues, ceilingvalues = (
            label_np(spherelist[1:-1], option)[0],
            label_np(spherelist[1:-1], option)[1],
        )
        return stripedvalues, ceilingvalues
    # Separate out the anisotropy for the Janus variant
    elif option == "Janus":
        topvalues, bottomvalues = (
            label_np(spherelist[1:-1], option)[0],
            label_np(spherelist[1:-1], option)[1],
        )
        return topvalues, bottomvalues


def rotation_matrix_from_vectors(vec_1: list[float],
                                 vec_2: list[float]):
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
    a, b = (vec_1 / np.linalg.norm(vec_1)).reshape(3), (vec_2 / np.linalg.norm(vec_2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    k_mat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + k_mat + k_mat.dot(k_mat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def pandas_np(ligand_string : str,
              first_atom : str,
              last_atom : str,
              sphere_list : list[list[float]],
              ligand_name : str,
              core_name : str,
              length : float = 1.0):
    """Placeholder
    """
    transformation_list, name_list = [], []
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

    unit_vector = np.linalg.norm(ligand_alignment_vector)
    vec_ligand = ligand_alignment_vector.tolist()

    # Loop over the sphere and find the
    for index in range(0, len(sphere_list)):
        vec_2 = sphere_list[index]
        # Find the transformationvector for the ligand vector to vec2, which is the position of the point on sphere
        transformation_vector = rotation_matrix_from_vectors(vecligand, vec2)
        # Rotate the vector
        vec_1_rot = transformationvector.dot(
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
        x_plotsphere.append(entry[0])
        y_plotsphere.append(entry[1])
        z_plotsphere.append(entry[2])
        sphere_name.append("P5")

    df_ligand = pd.DataFrame(
        list(zip(x_plot, y_plot, z_plot, ligands)), columns=["X", "Y", "Z", "NAME"]
    )

    df_core = pd.DataFrame(
        list(zip(x_plot_sphere, y_plot_sphere, z_plot_sphere, sphere_name)),
        columns=["X", "Y", "Z", "NAME"],
    )

    df_ligand["RESNAME"] = ligand_name
    df_core["RESNAME"] = cor_ename
    return df_ligand, df_core


def pandas_np_martini(
    molecule,
    ligand_alignment_vector,
    transformation_list,
    sphere_list,
    ligand_name,
    core_name,
    length=1.0,
):
    """ Function to read Martini molecule information and orientate on NP surface"""

    ligand_list, name_list = [], []
    sphere = []
    x_plot, y_plot, z_plot = [], [], []
    x_plot_sphere, y_plot_sphere, z_plot_sphere = [], [], []
    # Sulfur/ligand vector
    unit_vector = np.linalg.norm(ligand_alignment_vector)
    vec_1 = ligand_alignment_vector.tolist()
    for index in range(0, len(sphere_list)):
        vec_2 = sphere_list[index]
        transformation_vector = rotation_matrix_from_vectors(vec_1, vec_2)
        vec_1_rot = transformation_vector.dot(
            vec_1
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

    df_ligand["RESNAME"] = ligand_name

    df_core = pd.DataFrame(
        list(zip(x_plot_sphere, y_plot_sphere, z_plot_sphere, sphere_name)),
        columns=["X", "Y", "Z", "NAME"],
    )

    df_core["RESNAME"] = core_name
    return df_ligand, df_core


def read_martini_molecules(grofile: str, first: str, last: str):
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
    transformation_list = []
    martini_universe = mda.Universe(grofile)  # Load the Martini gro file in as a universe
    ids = [i.name for i in martini_universe.atoms]
    molecule = martini_universe.select_atoms("all")
    # In this case, the atoms will be N1 and R3
    first_atom = molecule.select_atoms("name {}".format(first))
    last_atom = molecule.select_atoms("name {}".format(last))
    ligand_alignment_vector = (first_atom.positions - last_atom.positions)[
        0
    ]  # Get the alignment vector created from the first and COM

    # Loop over the positions
    for i, j in enumerate(molecule.positions):
        vector = (j - first_atom.positions)[0]
        vector[0] = ligand_alignment_vector[0] - vector[0]
        vector[1] = ligand_alignment_vector[1] - vector[1]
        vector[2] = ligand_alignment_vector[2] - vector[2]

        if vector[0] == -math.inf:
            pass
        if vector[0] == 0.0:
            pass
        else:
            transformation_list.append([vector, molecule.atoms[i].type])

    # Return the universe, the transformed (normalized) coordinate list of the ligand molecule, and the
    # alignment vector that shows the arrow of direction of the vector, which we will be able to reorientate
    return molecule, transformation_list, ligand_alignment_vector
