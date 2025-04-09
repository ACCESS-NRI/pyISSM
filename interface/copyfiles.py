"""
This script copies over relevant Pyton and MATLAB files from the existing ISSM repository.

Author: Lawrence Bird
Email: lawrence.bird@anu.edu.au
"""
import os
import shutil

def copy_files(src_dir, dest_dir, file_extension):
    """
    Recursively copies files with a specific extension from the source directory to
    the destination directory, preserving the original directory structure.

    Parameters:
    - src_dir (str): Path to the source directory to copy files from.
    - dest_dir (str): Path to the destination directory where files will be copied.
    - file_extension (str): File extension to filter by (e.g., '.py' for Python files).

    Notes:
    - Only files matching the given extension will be copied.
    - Directory structure from the source will be replicated in the destination.
    - File metadata (such as modification times) is preserved during copying.
    - Creates any necessary subdirectories in the destination path.

    Example:
        copy_files('/path/to/source', '/path/to/destination', '.py')
    """

    for root, dirs, files in os.walk(src_dir):

        # Determine relative path from the source directory
        rel_path = os.path.relpath(root, src_dir)
        dest_path = os.path.join(dest_dir, rel_path)

        # Ensure the corresponding directory exists in destination
        os.makedirs(dest_path, exist_ok=True)

        # Copy only .py files
        for file in files:
            if file.endswith(file_extension):
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path, file)
                shutil.copy2(src_file, dest_file)