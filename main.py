import argparse

from contact_matrix_generator.generate_contact_matrix import ContactMap


def main():
    """
    Usage: python generate_contact_matrix.py -i [FILE] -o [DIR]
    Where [FILE] is an mmCIF file and [DIR] is the output directory

    :return: None
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_path",
        required=True,
        help="Path to a single mmCIF file",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        required=True,
        help="The path to where the output .npy files will be saved",
    )
    args = parser.parse_args()

    cm = ContactMap(args.input_path, args.output_path)
    cm.run()


if __name__ == "__main__":
    main()
