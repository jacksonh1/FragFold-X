from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextlib
import os
import sys

import pymol
from pymol import cmd


@contextlib.contextmanager
def suppress_pymol_output():
    """Temporarily silence PyMOL's C-level stdout/stderr chatter."""
    sys.stdout.flush()
    sys.stderr.flush()
    with open(os.devnull, "w") as devnull:
        old_stdout_fd = os.dup(sys.stdout.fileno())
        old_stderr_fd = os.dup(sys.stderr.fileno())
        try:
            os.dup2(devnull.fileno(), sys.stdout.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                yield
        finally:
            os.dup2(old_stdout_fd, sys.stdout.fileno())
            os.dup2(old_stderr_fd, sys.stderr.fileno())
            os.close(old_stdout_fd)
            os.close(old_stderr_fd)


with suppress_pymol_output():
    pymol.finish_launching(['pymol', '-cq'])  # Quiet, no GUI
    cmd.feedback("disable", "all", "actions results details warnings errors")



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

    cmd.reinitialize()

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
    # cmd.quit()

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


def color_residues_pymol(df, pdb_file, chain='B', colormap='magma', session_name=None, vmax=None, vmin=None):
    """
    Color residues in PyMOL based on residue-value mapping for a specific chain
    df must have columns 'resi' and 'value' for residue number and value to color by
    """
    cmd.reinitialize()
    cmd.load(pdb_file)
    cmap = plt.get_cmap(colormap)
    if vmax is None:
        vmax = df['value'].max()
    if vmin is None:
        vmin = df['value'].min()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    for idx, row in df.iterrows():
        resi = row['resi']
        value = row['value']
        color_rgba = cmap(norm(value))
        color_rgb = color_rgba[:3]  # Remove alpha channel
        color_name = f"res_{resi}_color"
        cmd.set_color(color_name, color_rgb)
        selection = f"resi {resi} and chain {chain}"
        cmd.color(color_name, selection)
    cmd.show("cartoon", f"chain {chain}")
    cmd.hide("lines")
    cmd.color("gray80", f"not chain {chain}")
    if session_name is None:
        session_filename = Path(pdb_file).stem + f"_{chain}_colored.pse"
        session_name = Path(pdb_file).parent / session_filename
    cmd.save(session_name)
    print(f"PyMOL session saved as: {session_name}")


def align_and_get_chain_rmsd(pdb_file1, pdb_file2, chain1='A', chain2='A'):
    """
    Aligns two PDB files on specified chains and returns the RMSD between those chains.

    Parameters
    ----------
    pdb_file1 : str or pathlib.Path
        Path to the first PDB file (reference).
    pdb_file2 : str or pathlib.Path
        Path to the second PDB file (to be aligned).
    chain1 : str
        Chain identifier in the first PDB file.
    chain2 : str
        Chain identifier in the second PDB file.

    Returns
    -------
    float
        The RMSD between the aligned chains.
    """
    cmd.reinitialize()
    obj1 = "obj1"
    obj2 = "obj2"
    cmd.load(str(pdb_file1), obj1)
    cmd.load(str(pdb_file2), obj2)
    sel1 = f"{obj1} and chain {chain1}"
    sel2 = f"{obj2} and chain {chain2}"
    result = cmd.align(sel2, sel1)
    rmsd = result[0]  # RMSD is the first element
    cmd.delete("all")
    return rmsd