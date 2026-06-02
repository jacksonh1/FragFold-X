# %%
"""
inputs:
- receptor MSAs - <list or path> (if isinstance, convert to list)
- receptor slice coords (1 based)
- fragment source - <path>
- fragment slice coords (1 based)
- stride
- fragment length
- msa_cache_dir
- model weights
- output_directory
- reference pdb (optional)

"""

from pathlib import Path
import a3mcat
# import os

# import fragfold3
import fragfold3.tools.plotting as plotting
import matplotlib.pyplot as plt
import fragfold3.pipeline.parameters as ffparams
import fragfold3.tools.colabfold_tools as colabfold_tools
import fragfold3.tools.sequence_utils as seq_utils
import fragfold3.job_schedulers.slurm_job_submitter as slurm_job_submitter
import fragfold3.tools.pymol_utils as pymol_utils
import fragfold3.pipeline.result_summary as result_summary
import shutil
import time
from loguru import logger


class EmptyMSA:

    def __init__(self, fasta: str | Path, domain_start: int, domain_end: int):
        self.fasta = Path(fasta)
        protein = a3mcat.import_fasta(self.fasta)
        if len(protein) != 1:
            raise ValueError(
                f"Expected exactly one sequence in {self.fasta} to build empty a3m, found {len(protein)}."
            )
        protein = protein[0]
        self.name = protein.header
        self.domain_start = domain_start
        self.domain_end = domain_end
        self.msa = a3mcat.MSAa3m.empty_MSA(sequence=protein.seq_str)
        self.sliced_msa = self.msa[domain_start : domain_end + 1]

    def __repr__(self):
        return f"EmptyMSA(\nfasta={self.fasta}, \ndomain_start={self.domain_start}, \ndomain_end={self.domain_end}\nname={self.name})"



class DomainMSA:

    def __init__(self, msa_file: str | Path, domain_start: int, domain_end: int):
        self.msa_file = Path(msa_file)
        self.name = self.msa_file.stem
        self.domain_start = domain_start
        self.domain_end = domain_end
        self.msa = self._import_a3m()
        # if domain_start <0:
        #     domain_start = 0
        # if domain_end <0:
        #     domain_end = len(self.msa.query.seq_str) - 1
        self.sliced_msa = self.msa[domain_start : domain_end + 1]

    def _import_a3m(self) -> a3mcat.MSAa3m:
        """
        Import the receptor MSA using a3mcat.
        """
        if not self.msa_file.exists():
            raise FileNotFoundError(f"MSA file {self.msa_file} does not exist.")
        a3m = a3mcat.MSAa3m.from_a3m_file(self.msa_file)
        return a3m

    def __repr__(self):
        return f"DomainMSA(\nmsa_file={self.msa_file}, \ndomain_start={self.domain_start}, \ndomain_end={self.domain_end})"

# ==============================================================================
# // colabfold related functions
# ==============================================================================

def get_colabfold_msa(
    fasta_file: str | Path, msa_cache_dir: str | Path, **kwargs
) -> Path:
    fasta_file = Path(fasta_file)
    seqs = seq_utils.import_fasta(fasta_file)
    if len(seqs) != 1:
        raise ValueError(
            f"Input fasta file {fasta_file} must contain exactly one sequence."
        )
    seqid = str(seqs[0].id)  # type: ignore
    # filestem = fasta_file.stem
    msa_file = Path(msa_cache_dir) / f"{seqid}.a3m"
    if not msa_file.exists():
        # Generate MSA using colabfold or any other method
        # This is a placeholder for the actual MSA generation logic
        colabfold_tools.colabfold_batch_MSA_wrapper(
            input_file=fasta_file,
            output_dir=msa_cache_dir,
            **kwargs,  # pass any additional arguments needed for MSA generation
        )
        if not msa_file.exists():
            raise FileNotFoundError(
                f"MSA file {msa_file} was not found despite attempting download. \n check {msa_cache_dir} and {fasta_file} sequence id"
            )
    else:
        print(
            f"MSA file ({msa_file}) with same name found in {msa_cache_dir}, skipping download."
        )
    return msa_file


def run_colabfold_batch(
    a3m_input_dir: str | Path, output_dir: str | Path, weights: str, **kwargs
):
    a3m_input_dir = Path(a3m_input_dir)
    colabfold_tools.colabfold_batch_wrapper(
        input_file_or_directory=a3m_input_dir,
        output_dir=output_dir,
        weights=weights,  # type: ignore
        **kwargs,  # pass any additional arguments needed for colabfold
    )
    # for input_file in a3m_input_dir.glob("*.a3m"):
    #     done_file = Path(output_dir) / f"{input_file.stem}.done.txt"
    #     if done_file.exists():
    #         print(
    #             f"Skipping {input_file} as {done_file} already exists. Remove it to rerun."
    #         )
    #         continue
    #     colabfold_tools.colabfold_batch_wrapper(
    #         input_file_or_directory=input_file,
    #         output_dir=output_dir,
    #         weights=weights,  # type: ignore
    #         **kwargs,  # pass any additional arguments needed for colabfold
    #     )


def get_colabfold_batch_commands(
    a3m_files_or_directories: list[Path],
    output_dir: str | Path,
    parameters: ffparams.Fragfold3Params,
):
    # param_list = []
    commands = []
    for input_path in a3m_files_or_directories:
        param_dict = {
            "input_file_or_directory": input_path,
            "output_dir": output_dir,
            "weights": parameters.model_weights,  # type: ignore
            "colabfold_executable": parameters.colabfold_batch,
            "colabfold_data": parameters.colabfold_data,
        }
        param_dict["extra_args"] = parameters.extra_colabfold_params
        command = colabfold_tools.colabfold_batch_wrapper(
            **param_dict,
            run=False,
        )
        commands.append(command)
        # param_list.append(param_dict)
    return commands

# ==============================================================================
# // generate fragment indices
# ==============================================================================

def gen_fragment_indices_by_overlap(
    length: int, overlap: int, fragment_length: int
) -> list[tuple[int, int]]:
    """
    Generate overlapping window indices for a given length, overlap, and fragment length.
    Returns a list of tuples representing the start and end indices (0-based).

    At the end of the sequence, if the last fragment does not fit, the start index
    will be adjusted to ensure the fragment length is maintained. So the last fragment
    may have more overlap with the previous fragment than the specified overlap.

    parameters:
    ----------
    length : int
        The total length of the sequence.
    overlap : int
        The number of residues to overlap between fragments.
    fragment_length : int
        The length of the fragment to extract.
    """
    if overlap < 0 or fragment_length <= 0:
        raise ValueError(
            "Overlap must be non-negative and fragment length must be positive."
        )
    if length < fragment_length:
        raise ValueError("Length must be greater than or equal to fragment length.")

    indices = []
    i = 0
    b = False
    while not b:
        if i + fragment_length >= length:
            i = length - fragment_length
            b = True
        indices.append((i, i + fragment_length - 1))
        i += fragment_length - overlap
    return indices


def gen_fragment_indices_by_sliding(
    length: int, stride: int, fragment_length: int
) -> list[tuple[int, int]]:
    """
    Generate sliding window indices for a given length, stride, and fragment length.
    Returns a list of tuples representing the start and end indices (0-based).

    If the window hits the end of the sequence, the start index will be
    adjusted to ensure the fragment length is maintained. So the last fragment
    may have more overlap with the previous fragment than the specified stride.

    parameters:
    ----------
    length : int
        The total length of the sequence.
    stride : int
        The number of residues to slide the window each time.
    fragment_length : int
        The length of the fragment to extract.
    """
    if stride <= 0 or fragment_length <= 0:
        raise ValueError("Stride and fragment length must be positive integers.")
    if length < fragment_length:
        raise ValueError(
            f"Length ({length}) must be greater than or equal to fragment length ({fragment_length})."
        )
    indices = []
    # for start in range(0, length - fragment_length + 1, stride):
    #     end = start + fragment_length - 1
    #     indices.append((start, end))
    i = 0
    b = False
    while not b:
        if i + fragment_length >= length:
            i = length - fragment_length
            b = True
        indices.append((i, i + fragment_length - 1))
        i += stride
    return indices

# ==============================================================================
# // a3m preparation functions
# ==============================================================================

def prepare_receptor_msa(
    receptor_fastas: list[str | Path],
    receptor_slice_coords: list[tuple[int, int]],
    msa_cache_dir: str | Path,
    use_msas: bool = True,
    **kwargs,
):
    _receptor_msas = []
    for fasta_file, domain_coords in zip(receptor_fastas, receptor_slice_coords):
        if use_msas:
            msa_file = get_colabfold_msa(
                fasta_file=fasta_file,
                msa_cache_dir=msa_cache_dir,
                **kwargs,
            )
            domain_msa = DomainMSA(
                msa_file=msa_file,
                domain_start=domain_coords[0],
                domain_end=domain_coords[1],
            )
        else:
            domain_msa = EmptyMSA(
                fasta=fasta_file,
                domain_start=domain_coords[0],
                domain_end=domain_coords[1],
            )
        _receptor_msas.append(domain_msa)
    receptor_msa = _receptor_msas[0].sliced_msa
    receptor_name = _receptor_msas[0].name
    for domain_msa in _receptor_msas[1:]:
        receptor_msa += domain_msa.sliced_msa
        receptor_name += f"-{domain_msa.name}"
    return receptor_msa, receptor_name


def prepare_fragment_source_msa(
    fragment_source_fasta: str | Path,
    fragment_slice_coords: tuple[int, int],
    msa_cache_dir: str | Path,
    use_msa: bool = True,
    **kwargs,
):
    if use_msa:
        fragment_source_msa_file = get_colabfold_msa(
            fasta_file=fragment_source_fasta,
            msa_cache_dir=msa_cache_dir,
            **kwargs,
        )
        fragment_source_domainMSA = DomainMSA(
            msa_file=fragment_source_msa_file,
            domain_start=fragment_slice_coords[0],
            domain_end=fragment_slice_coords[1],
        )
    else:
        fragment_source_domainMSA = EmptyMSA(
            fasta=fragment_source_fasta,
            domain_start=fragment_slice_coords[0],
            domain_end=fragment_slice_coords[1],
        )
    fragment_source_msa = fragment_source_domainMSA.sliced_msa
    fragment_source_name = fragment_source_domainMSA.name
    return fragment_source_msa, fragment_source_name


def prepare_input_a3ms(params: ffparams.Fragfold3Params):
    af_input_dir = Path(params.output_directory) / "input_files"
    af_input_dir.mkdir(exist_ok=True, parents=True)
    # clear any existing a3m files in the input directory
    if params.overwrite:
        clear_input_a3ms(af_input_dir)
    # Slice coords are stored in the user's indexing base. Convert to 0-based here (and only
    # here) for slicing the MSAs, since a3mcat uses python slicing. The output filenames
    # below are built from the display-base coords directly, so they honor the input base.
    base = int(params.indexing_base)
    receptor_slice_coords_0 = [
        (start - base, end - base) for start, end in params.receptor_slice_coords
    ]
    fragment_slice_coords_0 = (
        params.fragment_slice_coords[0] - base,
        params.fragment_slice_coords[1] - base,
    )
    receptor_msa, receptor_name = prepare_receptor_msa(
        receptor_fastas=params.receptor_fastas,
        receptor_slice_coords=receptor_slice_coords_0,
        msa_cache_dir=params.msa_cache_dir,
        colabfold_executable=params.colabfold_batch,
        colabfold_data=params.colabfold_data,
        use_msas=params.use_receptor_msas,
    )
    fragment_source_msa, fragment_source_name = prepare_fragment_source_msa(
        fragment_source_fasta=params.fragment_source_fasta,
        fragment_slice_coords=fragment_slice_coords_0,
        msa_cache_dir=params.msa_cache_dir,
        colabfold_executable=params.colabfold_batch,
        colabfold_data=params.colabfold_data,
        use_msa=params.use_fragment_msa,
    )
    # print(f"fragment source msa: {fragment_source_msa}")
    if params.fragmentation_method == "overlap":
        fragment_indices = gen_fragment_indices_by_overlap(
            length=len(fragment_source_msa.query.seq_str),
            overlap=params.overlap_length,
            fragment_length=params.fragment_length,
        )
    elif params.fragmentation_method == "sliding_window":
        fragment_indices = gen_fragment_indices_by_sliding(
            length=len(fragment_source_msa.query.seq_str),
            stride=params.stride,
            fragment_length=params.fragment_length,
        )
    else:
        raise ValueError(
            f"Invalid fragmentation_method: {params.fragmentation_method}. Must be 'overlap' or 'sliding'."
        )
    a3m_files = []
    for fragment_index in fragment_indices:
        start, end = fragment_index
        name = f"{fragment_source_name}-{params.fragment_slice_coords[0] + start}to{params.fragment_slice_coords[0] + end}_vs_{receptor_name}"
        file_name = af_input_dir / f"{name}.a3m"
        a3m_files.append(file_name)
        if file_name.exists():
            logger.info(f"{file_name} exists. Skipping")
            continue
        fragment_msa = fragment_source_msa[start : end + 1]
        prediction_input_a3m = receptor_msa + fragment_msa
        prediction_input_a3m.save(file_name)
    return a3m_files, af_input_dir

# ==============================================================================
# // pre-colabfold processing functions
# ==============================================================================

def setup(
    params: ffparams.Fragfold3Params,
):
    main_output_dir = Path(params.output_directory)
    if main_output_dir.exists() and params.warn_output_exists:
        raise FileExistsError(
            f"Output directory {main_output_dir} already exists. Set `warn_output_exists=False` to suppress. exiting..."
        )
    main_output_dir.mkdir(exist_ok=True, parents=True)
    a3m_files, af_input_dir = prepare_input_a3ms(params)
    predictions_dir = main_output_dir / "predictions"
    predictions_dir.mkdir(exist_ok=True, parents=True)
    return {
        "a3m_files": a3m_files,
        "af_input_dir": af_input_dir,
        "predictions_dir": predictions_dir,
    }

# ==============================================================================
# // post-colabfold processing functions
# ==============================================================================

def clear_input_a3ms(input_a3m_dir: str | Path):
    """
    Clear the input a3m files in the specified directory.
    """
    input_a3m_dir = Path(input_a3m_dir)
    print(
        f"clearing {len(list(input_a3m_dir.glob('*.a3m')))} a3m files from {input_a3m_dir}."
    )
    for input_file in input_a3m_dir.glob("*.a3m"):
        input_file.unlink()  # remove the file


def align_pdbs(
    params: ffparams.Fragfold3Params,
) -> None:
    predictions_dir = Path(params.output_directory) / "predictions"
    logger.info(f"aligning pdb files in {predictions_dir}")
    # get the list of PDB files in the predictions directory
    pdb_files = list(predictions_dir.glob("*.pdb"))
    pymol_utils.align_pdbs(
        input_pdb_files=pdb_files,
        output_dir=predictions_dir,
    )
    logger.info(f"finished aligning pdb files in {predictions_dir}")


def cleanup(params: ffparams.Fragfold3Params):
    """
    Cleanup function to remove temporary files and directories.
    This is a placeholder for any cleanup logic you might want to implement.
    """
    af_input_dir = Path(params.output_directory) / "input_files"
    shutil.rmtree(af_input_dir)
    # clear_input_a3ms(af_input_dir)
    # consider removing the png files in the predictions directory
    predictions_dir = Path(params.output_directory) / "predictions"
    for png_file in predictions_dir.glob("*.png"):
        png_file.unlink()


def create_summary_csv(params: ffparams.Fragfold3Params):
    """
    Post-processing function to handle the results of the ColabFold batch predictions.
    """
    predictions_dir = Path(params.output_directory) / "predictions"
    logger.info(f"summarizing results in {predictions_dir}")
    result_summary.score_pdb_files_multiprocessing(
        input_directory=predictions_dir,
        output_file=Path(params.output_directory) / "structure_scores.csv",
        n_processes=params.structure_score_params.n_processes,
        chain_group_a=params.structure_score_params.chain_group_a,
        chain_group_b=params.structure_score_params.chain_group_b,
        distance_cutoff=params.structure_score_params.contact_distance_cutoff,
    )
    logger.info(f"finished summarizing results in {predictions_dir}")


def plot_results(params: ffparams.Fragfold3Params):
    """
    Generate plots from the summary CSV file.
    """
    summary_csv = Path(params.output_directory) / "structure_scores.csv"
    if not summary_csv.exists():
        raise FileNotFoundError(
            f"Summary CSV file {summary_csv} does not exist. Run create_summary_csv first."
        )
    filename1 = Path(params.output_directory) / "position_plot.html"
    receptors = [i.stem for i in params.receptor_fastas]
    fragment_source = params.fragment_source_fasta.stem
    fig, df = plotting.plotly_fragfold_results(
        summary_csv,
        fragment_source_label=fragment_source,
        receptor_labels=receptors,
        xcol="fragment_start",
    )
    fig.write_html(filename1)
    title = f"fragment source: {fragment_source} - receptor(s): {' + '.join(receptors)}"
    filename2 = Path(params.output_directory) / "position_plot.png"
    fig, ax = plotting.plot_fragfold_results(
        df=df,
        title=title,
        xcol="fragment_start",
    )
    fig.savefig(filename2, dpi=300)
    plt.close(fig)


# ==============================================================================
# // main pipeline functions
# ==============================================================================

def check_all_done(input_a3ms: list, predictions_dir: str | Path):
    """
    Check if all colabfold predictions are done for the given input A3M files.
    """
    for a3m_file in input_a3ms:
        done_flag_file = Path(predictions_dir) / f"{a3m_file.stem}.done.txt"
        if not done_flag_file.exists():
            return False
    return True


def get_unfinished_a3m_files(input_a3ms: list, predictions_dir: str|Path):
    """
    for a list of a3m files, filter out the ones that are already done (presence of done flag in predictions_dir
    """
    unfinished_files = []
    for a3m_file in input_a3ms:
        done_flag_file = Path(predictions_dir) / f"{a3m_file.stem}.done.txt"
        if not done_flag_file.exists():
            unfinished_files.append(a3m_file)
    return unfinished_files


# @app.command()
def fragfold3_pipeline(
    params: ffparams.Fragfold3Params,
    root: str | Path | None = None,
    clean_files: bool = False,
):
    if root is not None:
        root = Path(root)
        params.convert_paths2abs(root=root)
    prepared_data = setup(params)
    a3m_files = prepared_data["a3m_files"]
    af_input_dir = prepared_data["af_input_dir"]
    predictions_dir = prepared_data["predictions_dir"]
    if check_all_done(input_a3ms=a3m_files, predictions_dir=predictions_dir):
        print("All jobs already done, skipping colabfold predictions.")
    else:
        run_colabfold_batch(
            a3m_input_dir=af_input_dir,
            output_dir=predictions_dir,
            weights=params.model_weights,  # type: ignore
            colabfold_executable=params.colabfold_batch,
            colabfold_data=params.colabfold_data,
            extra_args=params.extra_colabfold_params,
        )
    if params.reference_pdb is not None:
        if not Path(params.reference_pdb).exists():
            logger.warning(f"reference pdb {params.reference_pdb} does not exist, skipping copy.")
        else:
            new_filename = f"{params.reference_pdb.stem}-reference_structure.pdb"
            shutil.copy(params.reference_pdb, predictions_dir / new_filename)
    align_pdbs(params)
    create_summary_csv(params)
    if clean_files:
        cleanup(params)
    plot_results(params)
    main_output_dir = Path(params.output_directory)
    if root is not None:
        root = Path(root)
        params.convert_paths2relative(root=root)
    # TODO: consider adding a function to convert everything back to 1-based indexing if the user set 1-based indexing in the params
    # this would require renaming all of the output files
    params.save(main_output_dir / "fragfold_params.yaml")


def fragfold3_pipeline_scheduler(
    params: ffparams.Fragfold3Params,
    job_submitter: slurm_job_submitter.SlurmJobSubmitter = slurm_job_submitter.colabfold_sbatch_submitter,
    root: str | Path | None = None,
    max_jobs_allowed=2,
    clean_files: bool = False,
    **job_submitter_kwargs,
):
    """distribute the individual predictions across any number of nodes using SLURM."""
    if root is not None:
        root = Path(root)
        params.convert_paths2abs(root=root)
    logger.info("running setup")
    prepared_data = setup(params)
    logger.info("finished setup")
    a3m_files = prepared_data["a3m_files"]
    # af_input_dir = prepared_data["af_input_dir"]
    predictions_dir = prepared_data["predictions_dir"]
    logger.info("getting colabfold batch commands")
    logger.info("finished getting colabfold batch commands")
    if check_all_done(input_a3ms=a3m_files, predictions_dir=predictions_dir):
        logger.info("All jobs already done, skipping colabfold job submission.")
        all_done = True
    else:
        unfinished_a3m_files = get_unfinished_a3m_files(input_a3ms=a3m_files, predictions_dir=predictions_dir)
        colab_cmds = get_colabfold_batch_commands(
            a3m_files_or_directories=unfinished_a3m_files,
            output_dir=predictions_dir,
            parameters=params,
        )
        # submit the commands to the job scheduler
        job_submitter.watch_and_submit(colab_cmds, max_jobs_allowed=max_jobs_allowed, **job_submitter_kwargs)
        all_done = False
    timeout_limit = 60 * 60 * 6  # 6 hours
    start_time = time.time()
    while not all_done:
        all_done = check_all_done(input_a3ms=a3m_files, predictions_dir=predictions_dir)
        if not all_done:
            logger.info("Waiting for all jobs to finish...")
            time.sleep(5)
        if time.time() - start_time > timeout_limit:
            logger.info("Timeout reached. Exiting...")
            break
    if params.reference_pdb is not None:
        if not Path(params.reference_pdb).exists():
            logger.warning(f"reference pdb {params.reference_pdb} does not exist, skipping copy.")
        else:
            new_filename = f"{params.reference_pdb.stem}-reference_structure.pdb"
            shutil.copy(params.reference_pdb, predictions_dir / new_filename)
    align_pdbs(params)
    create_summary_csv(params)
    if clean_files:
        logger.info("cleaning up files")
        cleanup(params)
        logger.info("finished cleaning up files")
    plot_results(params)
    main_output_dir = Path(params.output_directory)
    if root is not None:
        root = Path(root)
        params.convert_paths2relative(root=root)
    params.save(main_output_dir / "fragfold_params.yaml")
