import os
from Bio.PDB import PDBParser, DSSP
import confix
import pandas as pd


def save_dict_to_tsv(data_dict, output_file):
    # Create a DataFrame from the dictionary
    df = pd.DataFrame.from_dict(data_dict, orient="index")
    df.rename(
        columns={
            "0": "DSSP index",
            "1": "Amino acid",
            "2": "Secondary structure",
            "3": "Relative ASA",
            "4": "Phi",
            "5": "Psi",
            "6": "NH–>O_1_relidx",
            "7": "NH–>O_1_energy",
            "8": "O–>NH_1_relidx",
            "9": "O–>NH_1_energy",
            "10": "NH–>O_2_relidx",
            "11": "NH–>O_2_energy",
            "12": "O–>NH_2_relidx",
            "13": "O–>NH_2_energy",
        }
    )
    # Save the DataFrame as a TSV file
    df.to_csv(output_file, sep="\t")


def run_dssp_on_folder(pdb_folder, output_folder):
    # Ensure the output_folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Get a list of .pdb files in the pdb_folder
    pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        pdb_path = os.path.join(pdb_folder, pdb_file)
        output_path = os.path.join(
            output_folder, f"{os.path.splitext(pdb_file)[0]}.dssp"
        )

        # Parse the PDB file using PDBParser
        p = PDBParser()
        structure = p.get_structure("pdb", pdb_path)
        model = structure[0]

        # Run DSSP on the model and save the output to a .dssp file
        dssp = dict(DSSP(model, pdb_path, dssp="mkdssp"))

        save_dict_to_tsv(dssp, output_path)


if __name__ == "__main__":
    pdb_folder = (
        confix.PATH_OUT_PDB_FILE
    )  # Replace with the path to the folder containing the .pdb files
    output_folder = (
        confix.PATH_OUT_DSSP_TO_TSV
    )  # Replace with the path to the folder where DSSP output will be saved
    run_dssp_on_folder(pdb_folder, output_folder)
