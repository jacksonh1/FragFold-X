import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO

from Bio.SeqRecord import SeqRecord
import fragfold3.config as env
import fragfold3.tools.sequence_utils as tools
from typing import Literal


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
    pairmode: Literal[
        "unpaired",
        "paired",
        "unpaired_paired",
    ] = "unpaired",
    colabfold_executable: str | Path = env.COLABFOLD_BATCH,
    colabfold_data: str | Path = env.COLABFOLD_DATA,
    extra_args: str | Path = "",
    run = True,
    # gpu_number: int = 0,
) -> str:
    # subprocess.run(
    #     "AVAILABLE_GPU=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | sort -n -k2 | head -n1 | cut -d, -f1)",
    #     shell=True,
    #     check=True,
    # )
    # subprocess.run("export CUDA_VISIBLE_DEVICES=$AVAILABLE_GPU", shell=True, check=True)
    # subprocess.run(f"export CUDA_VISIBLE_DEVICES={gpu_number}", shell=True, check=True)
    # colab_command = f"{colabfold_executable} --model-type alphafold2_ptm --amber --num-relax 5 --use-gpu-relax --data '{colabfold_data}' --pair-mode unpaired {input_file_or_directory} {output_dir}"
    # colab_command = f"export CUDA_VISIBLE_DEVICES={gpu_number}; {colabfold_executable} --model-type {weights} --data '{colabfold_data}' --pair-mode {pairmode} {extra_args} {input_file_or_directory} {output_dir}"
    colab_command = f"{colabfold_executable} --model-type {weights} --data '{colabfold_data}' --pair-mode {pairmode} {extra_args} '{input_file_or_directory}' '{output_dir}'"
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
    subprocess.run("export MPLBACKEND=Agg", shell=True, check=True)
    colab_command = f'{colabfold_executable} --msa-only --data "{colabfold_data}" --msa-only "{input_file}" "{output_dir}"'
    subprocess.run(colab_command, shell=True, check=True)


def mafft_align_wrapper(
    input_seqrecord_list: list[SeqRecord],
    # mafft_executable: str = env.MAFFT_EXECUTABLE,
    mafft_executable: str,
    extra_args: str = "",
    n_align_threads: int = 8,
    output_format: str = "dict",
) -> tuple[str, dict | list]:
    # example extra_args: "--retree 1"
    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    # run mafft
    alignment_filename = f"{temp_file.name}-mafft.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")
    else:
        mafft_command = f'{mafft_executable} --thread {n_align_threads} --quiet --anysymbol {extra_args} "{temp_file.name}" > "{alignment_filename}"'
    # print(mafft_command)
    subprocess.run(mafft_command, shell=True, check=True)
    mafft_output = tools.import_fasta(alignment_filename, output_format=output_format)
    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return mafft_command, mafft_output  # type: ignore


# def clustal_align_wrapper(
#     input_seqrecord_list,
#     alignment_type="basic",
#     output_type="list",
#     n_align_threads: int = 8,
# ):
#     assert output_type in [
#         "list",
#         "dict",
#         "alignment",
#     ], f'`output_type` must be one of ["list", "dict", "alignment"]'
#     assert alignment_type in [
#         "basic",
#         "full",
#     ], f'`output_type` must be one of ["basic", "full"]'
#     # create temporary file
#     temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
#     # write seqrecords to temporary file
#     SeqIO.write(input_seqrecord_list, temp_file, "fasta")
#     temp_file.close()
#     alignment_filename = f"{temp_file.name}-clustal.fa"
#     # raise an error if the alignment file already exists. (it won't but just in case)
#     if os.path.exists(alignment_filename):
#         raise FileExistsError(f"{alignment_filename} already exists")

#     if alignment_type == "basic":
#         clustal_command = f'clustalo -i "{temp_file.name}" -o "{alignment_filename}" -v --outfmt=fa --threads={n_align_threads}'
#     # elif alignment_type == "full":
#     else:
#         clustal_command = f'clustalo -i "{temp_file.name}" -o "{alignment_filename}" -v --outfmt=fa --full --threads={n_align_threads}'
#     subprocess.run(clustal_command, shell=True, check=True)

#     # read in clustal output
#     if output_type == "list":
#         clustal_output = tools.import_fasta(alignment_filename, output_format="list")
#     elif output_type == "dict":
#         clustal_output = tools.import_fasta(alignment_filename, output_format="dict")
#     # elif output_type == "alignment":
#     else:
#         clustal_output = AlignIO.read(alignment_filename, "fasta")
#     # delete temporary file
#     os.remove(alignment_filename)
#     os.remove(temp_file.name)
#     return clustal_output


# def muscle_align_wrapper(
#     input_seqrecord_list: list[SeqRecord],
#     muscle_binary: str = "/Users/jackson/tools/muscle/muscle-5.1.0/src/Darwin/muscle",
#     output_type: str = "list",
#     n_align_threads: int = 8,
# ) -> list[SeqRecord] | dict[str, SeqRecord] | Align.MultipleSeqAlignment:
#     assert output_type in [
#         "list",
#         "dict",
#         "alignment",
#     ], f'`output_type` must be one of ["list", "dict", "alignment"]'

#     temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
#     SeqIO.write(input_seqrecord_list, temp_file, "fasta")
#     temp_file.close()
#     alignment_filename = f"{temp_file.name}-muscle.fa"
#     # raise an error if the alignment file already exists. (it won't but just in case)
#     if os.path.exists(alignment_filename):
#         raise FileExistsError(f"{alignment_filename} already exists")

#     muscle_command = f'{muscle_binary} -super5 "{temp_file.name}" -output "{alignment_filename}" -threads {n_align_threads}'
#     subprocess.run(muscle_command, shell=True, check=True)

#     if output_type == "list":
#         muscle_output = tools.import_fasta(alignment_filename, output_format="list")
#     elif output_type == "dict":
#         muscle_output = tools.import_fasta(alignment_filename, output_format="dict")
#     # elif output_type == "alignment":
#     else:
#         muscle_output = AlignIO.read(alignment_filename, "fasta")

#     # delete temporary file
#     os.remove(alignment_filename)
#     os.remove(temp_file.name)
#     return muscle_output
