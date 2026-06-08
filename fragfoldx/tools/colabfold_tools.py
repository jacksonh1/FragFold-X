import re
import os
from typing import Literal
import json
from pathlib import Path
import subprocess
import fragfoldx.config as env


def colabfold_batch_wrapper(
    input_file_or_directory: str | Path,
    output_dir: str | Path,
    weights: Literal[
        "alphafold2",
        "alphafold2_ptm",
        "alphafold2_multimer_v1",
        "alphafold2_multimer_v2",
        "alphafold2_multimer_v3",
        "deepfold_v1",
    ] = "alphafold2_ptm",
    colabfold_executable: str | Path = env.COLABFOLD_BATCH,
    colabfold_data: str | Path = env.COLABFOLD_DATA,
    extra_args: str = "",
    run = True,
    # gpu_number: int = 0,
) -> str:
    # pair-mode is fixed to "unpaired": fragfoldx builds its inputs with a3mcat, which only
    # produces unpaired MSAs, so paired mode would silently produce wrong predictions.
    # `extra_args` is a raw string of additional colabfold_batch flags (e.g. "--num-models 3").
    colab_command = f"{colabfold_executable} --model-type {weights} --data '{colabfold_data}' --pair-mode unpaired {extra_args} '{input_file_or_directory}' '{output_dir}'"
    if run:
        subprocess.run(colab_command, shell=True, check=True)
        return colab_command
    else:
        return colab_command


def colabfold_batch_MSA_wrapper(
    input_file: str | Path,
    output_dir: str | Path,
    colabfold_executable: str = env.COLABFOLD_BATCH,
    colabfold_data: str = env.COLABFOLD_DATA,
):
    # subprocess.run("export MPLBACKEND=Agg", shell=True, check=True)
    os.environ["JAX_PLATFORMS"] = "cpu"
    colab_command = f'{colabfold_executable} --msa-only --data "{colabfold_data}" --msa-only "{input_file}" "{output_dir}"'
    subprocess.run(colab_command, shell=True, check=True)



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


def parse_pdb_filename_general(pdb_filename, pattern=env.COLABFOLD_PDB_PREDICTION_FILENAME_REGEX):
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
