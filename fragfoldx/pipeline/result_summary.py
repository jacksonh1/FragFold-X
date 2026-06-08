"""
input:
- input directory
- reference pdb
- output file (defaults to "structure_scores.csv" in the input directory)


find all pdb files in the input directory
calculate a variety of scores for each pdb file
- which scores to calculate and their parameters should be adjustable

should this be a separate pipeline?


output:
- structure_scores.csv (in the input directory)
    - pdb information
    - fragment position
    - name for fragment source and receptor
    - parameter file
    - dockQ score (vs reference)
        - this will require defining the chains in the reference pdb
        - for model, the peptide chain is the last chain and the receptor chains are everything before it
    - rank (from pdb filename)
    - ipTM score
"""
# x = "Q96F46-33to62_vs_Q16552-Q16552_unrelaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb"
# regex = r"(?P<fragment_protein>.+)-(?P<fragment_start>\d+)to(?P<fragment_end>\d+)_vs_(?P<receptor_proteins>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d\.pdb"
# p = re.compile(regex)
# m = p.match(x)
# print(
#     m.group("fragment_protein"),
#     int(m.group("fragment_start")),
#     int(m.group("fragment_end")),
#     m.group("receptor_proteins"),
#     int(m.group("rank")),
# )

import os
import pandas as pd
from pathlib import Path
import fragfoldx.config as config
import fragfoldx.structure_scoring.weighted_contacts as weighted_contacts
import re
import tqdm
from functools import partial
import multiprocessing
from loguru import logger


def parse_pdb_filename(
    filename: str | Path,
    pdb_filename_regex=config.PDB_FILENAME_REGEX,
):
    filename = Path(filename)
    p = re.compile(pdb_filename_regex)
    m = p.match(filename.name)
    if m is None:
        raise ValueError(f"Filename {filename.name} does not match regex {pdb_filename_regex}")
    d = m.groupdict()
    d["fragment_start"] = int(d["fragment_start"])
    d["fragment_end"] = int(d["fragment_end"])
    d["rank"] = int(d["rank"])
    # Split on _vs_ separator, then split receptor proteins on '-' only between protein names (not within isoforms)
    # This assumes receptor_proteins is a string like "P15692-4-P15692-4" or "Q5VWK5-Q5VWK5"
    d["receptor_proteins"] = re.findall(r'[a-zA-Z0-9]+(?:-\d+)?', d["receptor_proteins"])
    return d


def score_pdb_files(
    input_directory: Path | str,
    output_file: Path | str | None = None,
    regex=config.PDB_FILENAME_REGEX,
):
    input_directory = Path(input_directory)
    if output_file is None:
        output_file = input_directory / "structure_scores.csv"
    res = []
    for p in tqdm.tqdm(list(input_directory.rglob("*.pdb"))):
        temp = parse_pdb_filename(p, regex)
        temp["fragment_center"] = ((temp["fragment_start"] + temp["fragment_end"]) / 2)# + 1  # type: ignore
        d = weighted_contacts.calculate_weighted_contacts(p)
        temp.update(d)
        temp["pdb_file"] = str(p.resolve())
        res.append(temp)
    structure_scores = pd.DataFrame(res)
    structure_scores = structure_scores.sort_values(["fragment_start", "rank"])


# parse_pdb_filename(
#     "Q16552-118to147_vs_Q96F46-Q8NAC3_unrelaxed_alphafold2_multimer_v3_model_2_seed_000.pdb"
# )

# %%


def score_pdb(
    pdb_file: Path | str,
    regex=config.PDB_FILENAME_REGEX,
    root: Path | str | None = None,
    distance_cutoff=3.5,
    **kwargs,
):
    pdb_file = Path(pdb_file)
    temp = parse_pdb_filename(pdb_file, regex)
    temp["fragment_center"] = ((temp["fragment_start"] + temp["fragment_end"]) / 2)# + 1 # type: ignore
    d = weighted_contacts.calculate_weighted_contacts(pdb_file, distance_cutoff=distance_cutoff, **kwargs)
    temp.update(d)
    temp["pdb_file"] = str(pdb_file.resolve())
    fragfold_processing_params = (
        pdb_file.resolve().parent.parent / "fragfold_params.yaml"
    )
    if fragfold_processing_params.exists():
        if root is not None:
            temp["fragfold_processing_params"] = str(
                fragfold_processing_params.resolve().relative_to(root)
            )
        else:
            temp["fragfold_processing_params"] = str(
                fragfold_processing_params.resolve()
            )
    else:
        logger.warning(
            f"{fragfold_processing_params} does not exist. No fragfold processing params found for {pdb_file.name}."
        )
    if root is not None:
        temp["pdb_file_relative"] = pdb_file.resolve().relative_to(root)
    return temp



def get_cpu_count():
    # Try to respect SLURM allocations
    if "SLURM_NTASKS" in os.environ:
        return int(os.environ["SLURM_NTASKS"])
    else:
        # Fallback to all CPUs on the node
        return multiprocessing.cpu_count()


def score_pdb_files_multiprocessing(
    input_directory: Path | str,
    output_file: Path | str | None = None,
    regex=config.PDB_FILENAME_REGEX,
    pdb_file_pat="*_rank_*.pdb",
    n_processes=None,
    distance_cutoff=4.0,
    **kwargs,
):
    input_directory = Path(input_directory)
    if output_file is None:
        output_file = input_directory / "structure_scores.csv"
    res = []
    pdb_files = [
        i for i in input_directory.rglob(pdb_file_pat) if "-aligned" not in i.name
    ]  # aligned?
    if len(pdb_files) == 0:
        logger.warning(f"No pdb files found in {input_directory}")
        return
    if n_processes is None:
        n_processes = get_cpu_count()
    logger.info(f"Using {n_processes} processes")
    # Use the "forkserver" start method rather than the default "fork": the parent process is
    # multi-threaded (tqdm's monitor thread, loguru), and forking a multi-threaded process is
    # unsafe (it can deadlock) and is deprecated as of Python 3.12. forkserver forks each worker
    # from a fresh, single-threaded server process, which sidesteps that and — unlike "spawn" —
    # imports the package only once (in the server) instead of in every worker.
    ctx = multiprocessing.get_context("forkserver")
    with ctx.Pool(n_processes) as p:
        results_iterator = p.imap_unordered(
            partial(
                score_pdb,
                regex=regex,
                distance_cutoff=distance_cutoff,
                **kwargs,
            ),
            pdb_files,
            chunksize=1,
        )
        for result in results_iterator:
            res.append(result)
    structure_scores = pd.DataFrame(res)
    # print(structure_scores.head())
    # breakpoint()
    structure_scores = structure_scores.sort_values(by=["fragment_start", "rank"])
    structure_scores.to_csv(output_file, index=False)


