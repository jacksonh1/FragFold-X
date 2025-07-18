from pathlib import Path


import pymol
from pymol import cmd

def align_pdbs(input_pdb_files, output_dir):
    """
    Parameters
    ----------
    input_pdb_files : list of str or pathlib.Path
        List of paths to input PDB files. The first file in the list is used as the reference structure.
    output_dir : str or pathlib.Path
        Directory where the aligned PDB files will be saved.

    Raises
    ------
    ValueError
        If no input PDB files are provided.

    Notes
    -----
    - Requires the PyMOL Python library to be installed and accessible.
    - The reference PDB file (the first in the list) is saved as-is to the output directory.
    - All other PDB files are aligned to the reference and saved to the output directory.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pymol.finish_launching(['pymol', '-cq'])  # Quiet, no GUI

    ref_name = "ref"
    input_pdb_files = [Path(f) for f in input_pdb_files]
    cmd.load(str(input_pdb_files[0]), ref_name)

    # Save the reference as is
    ref_output = output_dir / Path(input_pdb_files[0]).name
    cmd.save(str(ref_output), ref_name)

    for pdb_file in input_pdb_files[1:]:
        obj_name = pdb_file.stem
        cmd.load(str(pdb_file), obj_name)
        cmd.align(obj_name, ref_name)
        output_path = output_dir / pdb_file.name
        cmd.save(str(output_path), obj_name)
        cmd.delete(obj_name)

    cmd.delete("all")
    pymol.cmd.quit()

def align_pdbs_in_dir_and_overwrite(input_dir):
    """
    Aligns all PDB files in the specified directory to the first PDB file found in that directory.
    Overwrites the original PDB files with their aligned versions.

    Parameters
    ----------
    input_dir : str or pathlib.Path
        Directory containing PDB files to be aligned.

    Notes
    -----
    - The first PDB file found in the directory is used as the reference structure.
    """
    input_dir = Path(input_dir)
    pdb_files = list(input_dir.glob("*.pdb"))
    align_pdbs(pdb_files, input_dir)  # Overwrites original files
