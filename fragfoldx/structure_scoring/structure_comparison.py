"""
tools for comparing a predicted structure to a known structure
"""

import os
import subprocess
import tempfile
import fragfoldx.config as config


# USalign pdb1.pdb pdb2.pdb -mol prot -mm 1 -TMscore 7 -chimerax --script sup.cxc


def run_usalign(
    model_pdb_file,
    native_pdb_file,
    usalign_executable=config.USALIGN_EXECUTABLE,
    extra_args="-mol prot -do -TMscore 7 -outfmt 2",
):
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    command = f"{usalign_executable} {model_pdb_file} {native_pdb_file} {extra_args} > {temp_file.name}"
    subprocess.run(command, shell=True, check=True)
    with open(temp_file.name) as f:
        lines = f.readlines()
    os.remove(temp_file.name)
    keys = lines[0].strip().split("\t")
    values = lines[1].strip().split("\t")
    return dict(zip(keys, values))


def usalign_aln_rmsd(
    model_pdb_file,
    native_pdb_file,
    usalign_executable=config.USALIGN_EXECUTABLE,
    extra_args="-mol prot -do -TMscore 7 -outfmt 2",
):
    results = run_usalign(
        model_pdb_file, native_pdb_file, usalign_executable, extra_args
    )
    return float(results["RMSD"])


def usalign_aln_tm_score_estimate(
    model_pdb_file,
    native_pdb_file,
    usalign_executable=config.USALIGN_EXECUTABLE,
    extra_args="-mol prot -do -TMscore 7 -outfmt 2",
):
    """This is an estimate of the alignment TM-score. It's an estimate because
    I had to derive it from the TMscore returned by USalign, which only provides
    the score normalized to the length of one of the proteins and not the length
    of the alignment. The value used for d0 can't be changed in USalign, or
    rather that feature does not seem to work for me. So I've only recalculated
    the normalization factor using the alignment length, but have not changed 
    the d0 value, which is dependent on length in the way it's usually 
    calculated.
    """
    results = run_usalign(
        model_pdb_file, native_pdb_file, usalign_executable, extra_args
    )
    x = float(results["TM1"])*float(results["L1"])
    return x/float(results["Lali"])
    


# x = run_usalign(
# )
# print(x)
# do = 1.24*((float(x['L2'])-15)**(1/3)) - 1.8
# x = run_usalign(
#         extra_args=f"-d {do} -mol prot -do -TMscore 7 -outfmt 2",
# )
# print(x)
