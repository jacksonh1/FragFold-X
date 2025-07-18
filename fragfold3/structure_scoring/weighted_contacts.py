from Bio.PDB import PDBParser, NeighborSearch, Superimposer, Select # type: ignore
from pathlib import Path
import json
import re
import fragfold3.tools.colabfold_tools as colabfold_tools
import fragfold3.tools.pdb_tools as pdb_tools



def is_interchain_contact(
    res_1,
    res_2,
    chain_group_a: None | set | list = None,
    chain_group_b: None | set | list = None,
):
    """
    returns True if the residues are in different chains or in different chain groups (if defined)
    This function adapted from original FragFold https://github.com/swanss/FragFold
    """
    res1_chain = res_1.get_parent().id
    res2_chain = res_2.get_parent().id
    if chain_group_a is None and chain_group_b is None:
        return res1_chain != res2_chain
    assert (
        chain_group_a is not None and chain_group_b is not None
    ), f"if 1 chain group is defined, you must define the other. chain groups: {chain_group_a=}, {chain_group_b=}"
    if res1_chain in chain_group_a and res2_chain in chain_group_b:
        return True
    if res1_chain in chain_group_b and res2_chain in chain_group_a:
        return True
    return False


def get_interchain_contacts(
    structure,
    contact_distance=4.0,
    chain_group_a=None,
    chain_group_b=None,
):
    '''
    This function adapted from original FragFold https://github.com/swanss/FragFold
    '''
    ns = NeighborSearch([x for x in structure.get_atoms()])
    nearby_res = ns.search_all(contact_distance, "R")
    contacts = [
        (x, y)
        for x, y in nearby_res # type: ignore
        if is_interchain_contact(x, y, chain_group_a, chain_group_b)
    ]
    return contacts


def get_interchain_contacts_from_pdb(
    pdb_file: str | Path,
    distance_cutoff: float | int = 4.0,
    chain_group_a: list[str] | None = None,
    chain_group_b: list[str] | None = None,
):
    '''
    This function adapted from original FragFold https://github.com/swanss/FragFold
    '''
    pdb_file = Path(pdb_file)
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("s", pdb_file)
    contacts = pdb_tools.contacts_to_strings(
        get_interchain_contacts(
            s,
            contact_distance=distance_cutoff,
            chain_group_a=chain_group_a,
            chain_group_b=chain_group_b,
        )
    )
    return contacts


def calculate_weighted_contacts(
    pdb_file: str | Path,
    score_file: str | Path | None = None,
    distance_cutoff: float | int = 4.0,
    chain_groups: list[list[str]] | None = None,
):
    """get the interchain contacts, number of contacts, iptm, and iptm
    weighted number of contacts from a pdb file

    This function adapted from original FragFold https://github.com/swanss/FragFold
    though it has been modified quite a bit

    Parameters
    ----------
    pdb_file : str | Path
        pdb file of predicted structure
    score_file : str | Path | None, optional
        a json file with scores corresponding to the `pdb_file`, by default None.
        If None, the score file will be inferred from the pdb file name (using
        the colabfold naming convention) and assumed to be in the same
        directory as the `pdb_file`.
    distance_cutoff : float | int, optional
        The distance in angstroms between 2 residues to be considered a contact,
        by default 4.0
    chain_groups : list[list[str]] | None, optional
        The groups of chain ids to be considered "intermolecular", by default None. If None,
        the first chain will be considered group A and the rest group B. For a
        contact to be considered "intermolecular", the residues have to belong to
        different groups. If not None, the groups should be a list of 2 lists,
        where each inner list contains the chain ids. For example, to find
        contacts between chains A and B, you would pass chain_groups=[["A"], ["B"]].

    Returns
    -------
    dict
        dictionary with the interchain contacts ("contacts"), number of
        contacts ("n_contacts"), iptm ("iptm"), and iptm weighted number of
        contacts ("weighted_contacts")
    """
    pdb_file = Path(pdb_file)
    if score_file is None:
        score_file = pdb_file.parent / colabfold_tools.colabfold_pdb_filename_2_score_filename(pdb_file)
    with open(score_file) as f:
        score_data = json.load(f)
    iptm = score_data["iptm"]
    if chain_groups is None:
        chains = pdb_tools.get_chains_from_structure(pdb_file)
        chain_group_a = chains[:-1]
        chain_group_b = [chains[-1]]
    else:
        chain_group_a, chain_group_b = chain_groups
    res = get_interchain_contacts_from_pdb(
        pdb_file,
        distance_cutoff=distance_cutoff,
        chain_group_a=chain_group_a,
        chain_group_b=chain_group_b,
    )
    res_dict = {
        "contacts": res,
        "n_contacts": len(res),
        "iptm": iptm,
        "weighted_contacts": len(res) * iptm,
    }
    return res_dict
