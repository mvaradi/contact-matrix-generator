import logging
import math
from pathlib import Path

import numpy as np
from gemmi import Position, cif


class ContactMap(object):
    """
    The methods in this class parse an mmCIF file, calculate pairwise C-alpha
    distances per chains and saves the distances as numpy .npy matrices
    """

    def __init__(self, input_path, output_path):
        """

        :param input_path: String; path to a single mmCIF file
        :param output_path: String; path to directory to save the .npy files
        """

        self.input_path = input_path
        self.output_path = output_path
        self.distances = {}

    def check_atom(self, atom):
        """

        Checks if an atom is a C-alpha atom

        :param atom:
        :return: Boolean
        """
        return atom["group_PDB"] == "ATOM" and atom["label_atom_id"] == "CA"

    def get_atoms(self, mmcif):
        """
        Returns a row from the "_atom_site" loop of the mmcif data

        :param mmcif: Gemmi mmCIF object
        :return: Gemmi mmCIF row
        """

        return mmcif.find(
            "_atom_site.",
            [
                "group_PDB",
                "label_atom_id",
                "label_asym_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
            ],
        )

    def get_position(self, atom):
        """
        Returns the x, y, z coordinates of an atom

        :param atom: Gemmi mmCIF atom
        :return: Gemmi Position object
        """

        x = float(atom["Cartn_x"])
        y = float(atom["Cartn_y"])
        z = float(atom["Cartn_z"])
        return Position(x, y, z)

    def get_file_name(self, chain_id):
        """
        Returns a file name

        :param chain_id: String, a PDB chain identifier
        :return: String, a file name to save the matrix data
        """
        file_name = Path(self.input_path).parts[-1].split(".")[0]
        file_name_with_chain = f"{file_name}_{chain_id}_matrix.npy"
        return str(Path(self.output_path, file_name_with_chain))

    def save_matrices(self):
        """
        Save all the matrices per chain to individual .npy files

        :return: numpy matrix
        """
        for chain_id in self.distances.keys():
            dimension = int(math.sqrt(len(self.distances[chain_id])))
            distances_matrix = np.array(self.distances[chain_id]).reshape(
                dimension, dimension
            )
            self.save_matrix(self.get_file_name(chain_id), distances_matrix)

    def save_matrix(self, path, matrix):
        """
        Save a single matrix to an .npy file

        :param matrix: numpy matrix
        :param path: String, path to where the file should be saved
        :return: None
        """
        with open(path, "wb") as outfile:
            np.save(outfile, matrix)

    def get_distances(self, mmcif):
        """
        Calculate the pairwise distances between two atoms across an mmCIF file

        :param mmcif: Gemmi mmCIF object
        :return: None
        """
        for atom_1 in self.get_atoms(mmcif):
            chain_1, pos_1 = self.get_chain_and_position(atom_1)
            for atom_2 in self.get_atoms(mmcif):
                chain_2, pos_2 = self.get_chain_and_position(atom_2)
                if chain_1 and chain_2 and chain_1 == chain_2:
                    distance = round(pos_1.dist(pos_2), 2)
                    self.distances[chain_1].append(distance)
        print(self.distances)

    def get_chain_and_position(self, atom):
        """
        Get the chain identifier and the x, y, z coordinates of an atom

        :return: String, GEMMI Position object
        """
        if self.check_atom(atom):
            chain = atom["label_asym_id"]
            if chain not in self.distances.keys():
                self.distances[chain] = []
            position = self.get_position(atom)
            return chain, position
        return None, None

    def run(self):
        """
        Run the process:
        1.) Read a single mmCIF file
        2.) Calculate the pairwise distances between C-alpha atoms
        3.) Save the resulting matrices per chains in .npy format

        :return: None
        """
        try:
            mmcif = cif.read_file(self.input_path).sole_block()
            self.get_distances(mmcif)
            self.save_matrices()
            return True
        except Exception as e:
            logging.error("Error: %s" % e)
            return
