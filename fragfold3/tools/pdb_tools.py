import os
from Bio.PDB import PDBIO # type: ignore
from IPython.display import display
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder # type: ignore
from Bio.SeqIO.PdbIO import CifSeqresIterator
from pathlib import Path
import numpy as np
from Bio.SeqUtils import seq1
from Bio.Data import IUPACData
from Bio import PDB
import tempfile
import py3Dmol
from io import StringIO

def contacts_to_strings(contacts):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    return [_get_contact_name(rA, rB) for rA, rB in contacts]


def _get_residue_name(res):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    return f"{res.get_parent().id}_{res.id[1]}"


def _get_contact_name(resi, resj):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    # not sure why the order of the residues is important here.
    if resi.get_parent().id <= resj.get_parent().id:
        return f"{_get_residue_name(resi)}-{_get_residue_name(resj)}"
    else:
        return f"{_get_residue_name(resj)}-{_get_residue_name(resi)}"


def get_chains_from_structure(path):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    path = Path(path)
    parser = PDBParser()
    structure = parser.get_structure(Path(path).stem, path)
    if structure is None:
        raise ValueError(f"Could not parse structure from {path}")
    chain_id_list = [c.id for c in structure.get_chains()]
    return chain_id_list


def extract_chain_sequence_from_resnames(input_file, chain_id='B'):
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.pdb':
        parser = PDBParser(QUIET=True)
    elif file_extension == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please use .pdb or .cif files.")
    structure = parser.get_structure("structure", input_file)
    for chain in structure.get_chains(): # type: ignore
        if chain.id == chain_id:
            # Extract residues and convert to one-letter code
            residues = [residue for residue in chain if residue.get_id()[0] == " "]
            sequence = "".join(seq1(residue.get_resname()) for residue in residues)
            return sequence
    raise ValueError(f"Chain {chain_id} not found in PDB file {input_file}")


def extract_sequences_from_pdb(input_file):
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.pdb':
        parser = PDBParser(QUIET=True)
    elif file_extension == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please use .pdb or .cif files.")
    structure = parser.get_structure("structure", input_file)
    ppb = PPBuilder()
    sequences = {}
    for pp in ppb.build_peptides(structure[0]): # type: ignore
        chain_id = pp[0].get_parent().id
        sequences[chain_id] = str(pp.get_sequence())
    return sequences


def extract_chain_sequence_from_coords(input_file, chain_id='B'):
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.pdb':
        parser = PDBParser(QUIET=True)
    elif file_extension == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please use .pdb or .cif files.")

    structure = parser.get_structure(os.path.basename(input_file), input_file)
    ppb = PPBuilder()
    
    for model in structure: # type: ignore
        for chain in model:
            if chain.id == chain_id:
                full_sequence = ""
                for pp in ppb.build_peptides(chain):
                    full_sequence += str(pp.get_sequence())
                return full_sequence
    return None


def _extract_chain_sequence_from_seqres_pdb(pdb_file, chain_id='A'):
    sequence = ""
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("SEQRES") and line[11] == chain_id:
                sequence += ''.join(aa_three_to_one(aa) for aa in line[19:].split())
    if sequence == "":
        raise ValueError(f"Chain {chain_id} not found in {pdb_file}")
    return sequence


def _extract_chain_sequence_from_seqres_cif(cif_file, chain_id='A'):
    cif_sequences = {}
    with open(cif_file, "r") as handle:
        for record in CifSeqresIterator(handle):
            cif_sequences[record.annotations['chain']] = str(record.seq)
    if chain_id not in cif_sequences:
        raise ValueError(f"Chain {chain_id} not found in {cif_file}")
    return cif_sequences[chain_id]


def extract_chain_sequence_from_seqres(input_file, chain_id='A'):
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.pdb':
        return _extract_chain_sequence_from_seqres_pdb(input_file, chain_id)
    elif file_extension == '.cif':
        return _extract_chain_sequence_from_seqres_cif(input_file, chain_id)
    else:
        raise ValueError("Unsupported file format. Please use .pdb or .cif files.")


def get_chain_center_of_mass(pdb_file, chain_id):
    """
    CAUTION: This was generated by an LLM and hasn't been tested.
    Calculate the center of mass of a specific chain in a PDB file.

    Parameters:
    - pdb_file: Path to the PDB file.
    - chain_id: ID of the chain for which to calculate the center of mass.

    Returns:
    - center_of_mass: A numpy array representing the center of mass coordinates.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Structure", pdb_file)

    # Extract atoms from the specified chain
    chain_atoms = []
    for model in structure: # type: ignore
        for chain in model:
            if chain.get_id() == chain_id:
                for residue in chain:
                    for atom in residue:
                        chain_atoms.append(atom)

    # Calculate the center of mass
    coords = np.array([atom.get_coord() for atom in chain_atoms])
    center_of_mass = np.mean(coords, axis=0)
    return [float(coord) for coord in center_of_mass]


def aa_three_to_one(aa):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(aa, 'X')  # 'X' for unknown amino acids
 

def check_residue_numbering_with_pdb(pdb_file, residue, number, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    if res.get_id()[1] == number:
                        resname = res.get_resname().capitalize()
                        if resname in IUPACData.protein_letters_3to1:
                            one_letter = IUPACData.protein_letters_3to1[resname]
                        else:
                            raise ValueError(f"Unknown residue name: {resname}")
                        if one_letter == residue:
                            return True
                        else:
                            raise ValueError(f"Residue mismatch: expected {residue} ({number}), found {one_letter} ({res.get_id()[1]})")
    raise ValueError(f"Residue {residue} with number {number} not found in chain {chain_id} of the PDB file.")


def biopython_Structure_to_string(structure):
    """Convert a Biopython Structure object to a PDB string."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    io.save(temp_file.name)
    temp_file.close()
    with open(temp_file.name, "r") as f:
        pdb_string = f.read()
    os.remove(temp_file.name)  # Clean up the temporary file
    return pdb_string


def show_structure_with_py3Dmol(structure, style="cartoon"):
    """
    Displays a Biopython Structure object using py3Dmol in a Jupyter environment.
    """
    # Write structure to a PDB string
    pdb_buf = StringIO()
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_buf)
    pdb_str = pdb_buf.getvalue()

    # Visualize with py3Dmol
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_str, "pdb")
    view.setStyle({}, {style: {"colorscheme": "chain"}})
    view.zoomTo()
    display(view)


def show_pdb_file_with_py3Dmol(pdb_file, style="cartoon"):
    """
    Displays a PDB file using py3Dmol in a Jupyter environment.
    
    Parameters:
    - pdb_file: Path to the PDB file
    - style: Visualization style (default: "cartoon")
    """
    with open(pdb_file, 'r') as f:
        pdb_str = f.read()
    
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_str, "pdb")
    view.setStyle({}, {style: {"colorscheme": "chain"}})
    view.zoomTo()
    display(view)