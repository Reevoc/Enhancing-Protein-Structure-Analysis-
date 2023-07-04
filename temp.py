import os
import confix


def eliminate_first_3_chars_in_folder(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        return

    for root, dirs, files in os.walk(folder_path):
        for filename in files:
            # Get the new file name eliminating the first 3 characters
            new_filename = filename[3:]

            # Generate the old and new file paths
            old_file_path = os.path.join(root, filename)
            new_file_path = os.path.join(root, new_filename)

            try:
                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f"File '{old_file_path}' has been renamed to '{new_file_path}'.")
            except OSError as e:
                print(f"Error occurred while renaming the file: {e}")


if __name__ == "__main__":
    eliminate_first_3_chars_in_folder(confix.PATH_OUT_DSSP_FILE)
