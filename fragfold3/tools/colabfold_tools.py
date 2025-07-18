import re
import json
from pathlib import Path
from fragfold3 import config


def get_colabfold_scores(pdb_file: str | Path):
    pdb_file = Path(pdb_file)
    score_file = pdb_file.stem.replace("_unrelaxed", "_scores") + ".json"
    score_file = pdb_file.parent / score_file
    with open(score_file) as f:
        scores = json.load(f)
    return score_file, scores


def colabfold_pdb_filename_2_score_filename(pdb_filename):
    return pdb_filename.stem.replace("_unrelaxed", "_scores") + ".json"


def parse_rank_from_pdb_filename(
    pdb_filename, pattern=r".+_unrelaxed_rank_(\d\d\d)_.+"
):
    """
    extracts the rank from the name of the `pdb_file` using the provided regex pattern.

    Parameters
    ----------
    pdb_filename : str
        name of the pdb file
    pattern : regexp, optional
        regex pattern to extract structure rank from the `pdb_filename`,
        by default r".+_unrelaxed_rank_(\\d\\d\\d)_.+"
    """
    filename_pattern = re.compile(pattern)
    match = filename_pattern.match(pdb_filename)
    assert (
        match is not None
    ), f"File name {pdb_filename} does not match pattern {pattern}"
    rank = match.group(1)
    return int(rank)


def parse_pdb_filename_general(pdb_filename, pattern=config.COLABFOLD_PDB_PREDICTION_FILENAME_REGEX):
    """
    extracts the rank, filename and weights from the name of the `pdb_file` using the provided regex pattern.

    Parameters
    ----------
    pdb_filename : str
        name of the pdb file
    pattern : regexp, optional
        regex pattern to extract structure rank from the `pdb_filename`,
        by default r"(?P<name>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"
    """
    pdb_filename = Path(pdb_filename)
    p = re.compile(pattern)
    m = p.match(pdb_filename.name)
    if m is None:
        raise ValueError(f"Filename {pdb_filename.name} does not match pattern {pattern}")
    return m.groupdict()
