import os
import subprocess
import confix


def run_mkdssp(pdb_folder, dssp_output_folder):
    # Ensure the dssp_output_folder exists
    os.makedirs(dssp_output_folder, exist_ok=True)

    # Get a list of downloaded PDB files in the pdb_folder
    pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith(".ent")]

    # Run mkdssp for each PDB file
    count = 1
    for pdb_file in pdb_files:
        pdb_path = os.path.join(pdb_folder, pdb_file)
        pdb_code = os.path.splitext(pdb_file)[0]
        dssp_file = f"{pdb_code}.dssp"
        dssp_path = os.path.join(dssp_output_folder, dssp_file)

        # Run the mkdssp command
        try:
            subprocess.run(["mkdssp", "-i", pdb_path, "-o", dssp_path], check=True)
            print(f"Generated DSSP file for {pdb_code}. [{count},{len(pdb_files)}]")
        except subprocess.CalledProcessError as e:
            print(f"Error generating DSSP file for {pdb_code}: {e}")
        count += 1


if __name__ == "__main__":
    pdb_folder = (
        confix.PATH_OUT_PDB_FILE
    )  # Replace with the path to the folder containing the downloaded PDB files
    dssp_output_folder = (
        confix.PATH_OUT_DSSP_FILE
    )  # Replace with the path to the folder where DSSP files will be saved
    run_mkdssp(pdb_folder, dssp_output_folder)
