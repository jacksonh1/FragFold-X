from Bio.PDB import PDBParser, NeighborSearch, Superimposer, Select # type: ignore
from pathlib import Path
import json
import re
import fragfoldx.tools.colabfold_tools as colabfold_tools
import fragfoldx.tools.pdb_tools as pdb_tools
import fragfoldx.structure_scoring.contacts as contacts


def calculate_weighted_contacts(
    pdb_file: str | Path,
    score_file: str | Path | None = None,
    distance_cutoff: float | int = 4.0,
    chain_group_a: list[str] | None = None,
    chain_group_b: list[str] | None = None,
    # chain_groups: list[list[str]] | None = None,
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
    chain_group_a : list[str] | None, optional
        The chain ids to be considered group A, by default None. If None, all chains
        except the last chain will be considered group A and the last chain will
        be group B. If not None, this should be a list of chain ids and 
        chain_group_b must also be provided. For example, to find contacts 
        between chains A and B, you would pass chain_group_a=["A"] and 
        chain_group_b=["B"]. If the structure only has 2 chains, this will be the
        default behavior.
    chain_group_b : list[str] | None, optional
        The chain ids to be considered group B, by default None. If None, the last
        chain will be considered group B and all other chains will be group A. If
        not None, this should be a list of chain ids and chain_group_a must also
        be provided.

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
    if any([chain_group_a, chain_group_b]) and not all([chain_group_a, chain_group_b]):
        raise ValueError(f"If one chain_group is defined, both must be defined: {chain_group_a=}, {chain_group_b=}")
    if chain_group_a is None:
        chains = pdb_tools.get_chains_from_structure(pdb_file)
        chain_group_a = chains[:-1]
        chain_group_b = [chains[-1]]
    res = contacts.get_interchain_contacts_from_pdb(
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
        "contact_distance_cutoff": distance_cutoff,
        "chain_group_a": chain_group_a,
        "chain_group_b": chain_group_b,
    }
    return res_dict
