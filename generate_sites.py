#! /usr/bin/env python3

import ase
from ase.io import *
from ase import Atoms
import sys
import os
import math
import numpy as np


def recursive_neighbours():  

    '''
    Returns the minimum distance of all Si/Al atoms in the structure
    relative to a given Si/Al atom.
    '''

    old_neighbours, positions, symbols = find_neighbours(
        input_structure, Al_of_the_first_acid_site, silicon_expected_neighbours
    )  # find the first 4 oxygen neighbours

    found_already = list(old_neighbours)
    neighbours_by_Si_to_Si_distance = []
    expected = oxygen_expected_neighbours  # search for the next neighbours which should be Si/Al atoms

    while len(found_already) < number_of_atoms_in_structure - 1:

        current_neighbours = []

        for stuff in old_neighbours:

            neighbours, positions, symbols = find_neighbours(
                input_structure, stuff, expected
            )

            for n in neighbours:

                if (
                    n not in found_already
                    and n != hydrogen_acid_site
                    and n != oxygen_of_the_first_acid_site
                    and n != Al_of_the_first_acid_site
                ):  # only search through atoms whichs have not been labelled with neighbour information yet

                    found_already.append(n)
                    current_neighbours.append(n)

        old_neighbours = current_neighbours
        expected = expected_neighbours(
            symbols[0]
        )  # keep switching between 2 and 4 expected neighbours

        if (
            expected == silicon_expected_neighbours
        ):  # only keep info on the Silcon tetrahedra
            neighbours_by_Si_to_Si_distance.append(current_neighbours)

    return neighbours_by_Si_to_Si_distance


def find_neighbours(infile, atom, expected):  

    '''
    Returns the nearest neighbours of a given atom.
    '''

    atoms = read(infile)
    pos = atoms.get_positions()
    sym = atoms.get_chemical_symbols()

    distances = []
    neighbours = []
    positions = []
    symbols = []

    for stuff in atoms:

        if stuff.index != atom:

            dist = atoms.get_distance(
                atom, stuff.index, mic=True
            )  # mic stands for minimum images convention where the unit cell is reflected in all 6 faces to find to minimum distance between two given atoms

            symbol = sym[stuff.index]

            if dist < cut_off_distance:

                distances.append(dist)
                neighbours.append(stuff.index)
                positions.append(pos[stuff.index])
                symbols.append(sym[stuff.index])

    if not (len(neighbours) == expected):

        print(
            "warning "
            + str(expected)
            + " neighbours expected but "
            + str(len(neighbours))
            + " neighbours found "
            + " for atom with index "
            + str(atom)
            
        )

    return neighbours, positions, symbols


def expected_neighbours(symbol):

    '''
    Returns expected number of neighbours for a given atom type.
    '''

    if symbol == "Si" or symbol == "Al":

        x = silicon_expected_neighbours

    elif symbol == "O":

        x = oxygen_expected_neighbours

    elif symbol == "H":

        x = hydrogen_expected_neighbours

    return x


def replace_Si_with_Al(infile, Si):

    '''
    Replaces a Si atom with Al keeping original the indices intact.
    '''

    atoms = read(infile)
    pos = atoms.get_positions()
    sym = atoms.get_chemical_symbols()
    cell = atoms.get_cell()

    replaceed_atom_pos = pos[Si]
    sym.pop(Si)
    sym.insert(Si, "Al")

    atoms = Atoms(sym, positions=pos)
    atoms.set_cell(cell)
    atoms.set_pbc([True, True, True])

    atoms.write("Al_" + str(Si) + ".traj")


def add_hydrogen(infile, atom):

    '''
    Adds the hydrogen on the outer edge of the triangle made by two Si/Al stoms and an oxygen.
    '''

    atoms = read(infile)
    pos = atoms.get_positions()

    neighbours, positions, symbols = find_neighbours(
        infile, atom, oxygen_expected_neighbours
    )

    midpoint = (positions[0] + positions[1]) / 2  # midpoint between the Silica

    shift = (
        midpoint - pos[atom]
    )  # Oxygen atom minus the midpoint, shift 1A in the opposite direction so that the hydrogen ends up on the outer edge of the triagle made by two Si and an O
    shift /= np.linalg.norm(shift)

    new_atoms = Atoms("H", positions=[pos[atom] - shift])
    atoms.extend(new_atoms)

    atoms.write("in.traj")


def make_files_by_Si_to_Si_distance(neighbours_by_Si_to_Si_distance):

    '''
    Makes the top layer of folders according to neighbour number relative to a given site.
    '''

    folder_name = output_folder
    change_folder(folder_name)

    for i, sites in enumerate(neighbours_by_Si_to_Si_distance):

        folder_name = "neighbours_" + str(i)
        change_folder(folder_name)
        make_subfolders_for_neighbours(sites)
        os.chdir("..")


def make_subfolders_for_neighbours(sites):

    '''
    Make folders for each silicon and another subfolder for its oxygen neighbours. 
    Folders are named with the indices from the input file.
    '''

    for silicons in sites:

        folder_name = "Si_" + str(silicons)
        change_folder(folder_name)
        replace_Si_with_Al(input_structure, silicons)
        neighbours, positions, symbols = find_neighbours(
            input_structure, silicons, silicon_expected_neighbours
        )

        for oxygens in neighbours:

            folder_name = "O_" + str(oxygens)
            change_folder(folder_name)
            add_hydrogen("../Al_" + str(silicons) + ".traj", oxygens)
            os.chdir("..")

        os.chdir("..")


def change_folder(
    folder_name,
):  

    '''
    Change the folder and create a folder if the path doesnt exist already.
    '''

    if not os.path.isdir(folder_name):

        os.mkdir(folder_name)

    os.chdir(folder_name)


def check_for_fails():  

    '''
    Returns the indices for which an incorrect number of neighbours are found.
    '''

    atoms = read(input_structure)

    indicies = [
        [atom.index, atom.symbol]
        for atom in atoms
        if (atom.symbol == "H")
        or (atom.symbol == "Si")
        or (atom.symbol == "O")
        or (atom.symbol == "Al")
    ]

    exceptions = []

    for stuff in indicies:

        try:

            neighbours, positions, symbols = find_neighbours(
                input_structure, stuff[0], expected_neighbours(stuff[1])
            )

        except:

            exceptions.append(stuff)

        if len(neighbours) != expected_neighbours(stuff[1]):

            exceptions.append(stuff)

    return exceptions


#####################################
# Input parameters and required files
#####################################

# filepaths

input_structure = os.getcwd() + "/Input_structures/CHA_no_proton.traj"
output_folder = "Generated_structures"

# acid site information

number_of_atoms_in_structure = 108
oxygen_of_the_first_acid_site = 107  # information about the first acid site that we find the neighbours of
Al_of_the_first_acid_site = 7
hydrogen_acid_site = 108  # extra parameter for when there is a proton on the first site

# geometry and distance criteria

cut_off_distance = 2
silicon_expected_neighbours = 4
oxygen_expected_neighbours = 2
hydrogen_expected_neighbours = 1


##########
# Run it
##########

neighbours_by_Si_to_Si_distance = recursive_neighbours()
make_files_by_Si_to_Si_distance(neighbours_by_Si_to_Si_distance)
