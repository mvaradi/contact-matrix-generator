from pathlib import Path
from unittest import TestCase

from gemmi import cif

from contact_matrix_generator.generate_contact_matrix import ContactMap


class TestContactMap(TestCase):
    def test_get_file_name(self):
        """
        Should create a file name using the input path,
        the output path and the chain identifier
        :return: None
        """
        cm = ContactMap("/input/path/file.cif", "/output/path")
        file_name = cm.get_file_name("A")
        expected = str(Path("/", "output", "path", "file_A_matrix.npy"))
        self.assertEqual(file_name, expected)

        cm = ContactMap("/input/path/file.foo", "/output/path")
        file_name = cm.get_file_name("B")
        expected = str(Path("/", "output", "path", "file_B_matrix.npy"))
        self.assertEqual(file_name, expected)

    def test_check_atom(self):
        """
        Should return True for CA ATOMs and False for everything else
        :return:
        """
        cm = ContactMap("/input/path/file.cif", "/output/path")
        atom = {"group_PDB": "ATOM", "label_atom_id": "CA"}
        self.assertTrue(cm.check_atom(atom))

        atom = {"group_PDB": "HETATOM", "label_atom_id": "CA"}
        self.assertFalse(cm.check_atom(atom))

        atom = {"group_PDB": "ATOM", "label_atom_id": "CB"}
        self.assertFalse(cm.check_atom(atom))

    def test_get_position(self):
        test_file = str(Path("tests", "test.cif"))
        cm = ContactMap(test_file, "/output/path")
        mmcif = cif.read_file(cm.input_path).sole_block()
        row = mmcif.find(
            "_atom_site.",
            [
                "group_PDB",
                "label_atom_id",
                "label_asym_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
            ],
        )[0]
        position = cm.get_position(row)
        self.assertEqual(position[0], 17.023)
        self.assertEqual(position[1], -10.577)
        self.assertEqual(position[2], 32.291)

    def test_get_atoms(self):
        test_file = str(Path("tests", "test.cif"))
        cm = ContactMap(test_file, "/output/path")
        mmcif = cif.read_file(cm.input_path).sole_block()
        row = mmcif.find(
            "_atom_site.",
            [
                "group_PDB",
                "label_atom_id",
                "label_asym_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
            ],
        )[0]
        rows = cm.get_atoms(mmcif)
        self.assertEqual(rows[0][0], row[0])
        self.assertEqual(rows[0][1], row[1])
        self.assertEqual(rows[0][2], row[2])
        self.assertEqual(rows[0][3], row[3])
        self.assertEqual(rows[0][4], row[4])
        self.assertEqual(rows[0][5], row[5])

    def test_get_distances(self):
        test_file = str(Path("tests", "test.cif"))
        cm = ContactMap(test_file, "/output/path")
        mmcif = cif.read_file(cm.input_path).sole_block()
        cm.get_distances(mmcif)
        expected = {
            "A": [
                0.0,
                3.82,
                6.89,
                8.91,
                12.25,
                3.82,
                0.0,
                3.82,
                6.33,
                9.05,
                6.89,
                3.82,
                0.0,
                3.8,
                6.02,
                8.91,
                6.33,
                3.8,
                0.0,
                3.79,
                12.25,
                9.05,
                6.02,
                3.79,
                0.0,
            ]
        }
        self.assertEqual(cm.distances, expected)
