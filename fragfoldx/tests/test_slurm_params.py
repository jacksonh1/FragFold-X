"""Tests for the colabfold SLURM parameter config and job-submitter throttling."""

import os
from unittest import mock

import pytest
import yaml

from fragfoldx.pipeline import parameters as ffparams
from fragfoldx.pipeline.parameters import SlurmParameters, Fragfold3Params
from fragfoldx.job_schedulers.slurm_job_submitter import SlurmJobSubmitter


# --- SlurmParameters: verbatim sbatch flags + generated job name ---


def test_defaults():
    sp = SlurmParameters()
    assert sp.job_name_prefix == "ff"
    assert sp.max_jobs_allowed == 2
    assert sp.sbatch == {}


def test_sbatch_flags_verbatim_dashes_optional():
    sp = SlurmParameters(sbatch={"partition": "p", "--gres": "gpu:1", "cpus-per-task": 5})
    sp.job_name = "ff-x-1234abcd"
    flags = sp.sbatch_flag_dict()
    # user's sbatch flags passed through under their real names (leading dashes added if absent)
    assert flags["--partition"] == "p"
    assert flags["--gres"] == "gpu:1"
    assert flags["--cpus-per-task"] == "5"
    # fragfoldx-generated flags always present
    assert flags["--job-name"] == "ff-x-1234abcd"
    assert flags["--output"].endswith("ff-x-1234abcd_%A_%a.out")


def test_reserved_sbatch_key_raises():
    # job-name / output / error are set by fragfoldx; supplying them is an error
    for reserved in ("job-name", "--output", "error"):
        with pytest.raises(ValueError):
            SlurmParameters(sbatch={reserved: "x"})


def test_sbatch_flag_dict_requires_job_name():
    sp = SlurmParameters()  # job_name empty until Fragfold3Params sets it
    with pytest.raises(ValueError):
        sp.sbatch_flag_dict()


def test_job_name_generated_hybrid_and_unique():
    p1 = Fragfold3Params(fragment_source_fasta="x.fasta", target_fastas=["TNFa.fasta"])
    p2 = Fragfold3Params(fragment_source_fasta="x.fasta", target_fastas=["TNFa.fasta"])
    for p in (p1, p2):
        jn = p.colabfold_slurm.job_name
        assert jn.startswith("ff-TNFa-")
        assert len(jn.rsplit("-", 1)[1]) == 8  # uuid suffix
    assert p1.colabfold_slurm.job_name != p2.colabfold_slurm.job_name


def test_round_trip_drops_persisted_job_name():
    p = Fragfold3Params(fragment_source_fasta="x.fasta", target_fastas=["TNFa.fasta"])
    d = p.to_writable_dict()
    # job_name is persisted (provenance) but regenerated on reload (not fed back to init)
    assert "job_name" in d["colabfold_slurm"]
    p2 = Fragfold3Params.from_dict(d)
    assert p2.colabfold_slurm.job_name != p.colabfold_slurm.job_name


def test_unknown_colabfold_slurm_key_fails_loudly():
    with pytest.raises(TypeError):
        Fragfold3Params.from_dict(
            {
                "fragment_source_fasta": "x.fasta",
                "colabfold_slurm": {"nonsense_key": 1},
            }
        )


# --- cluster-defaults merge: user wins, extra deep-merges ---


def test_cluster_defaults_merge(tmp_path, monkeypatch):
    frag = tmp_path / "frag.fasta"
    frag.write_text(">frag\nMKKLLAAVGGSSPPQQRRTTYYWWFFEEDDNN\n")
    tgt = tmp_path / "TNFa.fasta"
    tgt.write_text(">TNFa\nVRSSSRTPSDKPVAHVVANPQAEGQLQW\n")
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        yaml.dump(
            {
                "fragment_source_fasta": "frag.fasta",
                "fragment_slice_coords": [1, 20],
                "target_fastas": ["TNFa.fasta"],
                "target_slice_coords": [[1, 20]],
                "output_directory": "out",
                "colabfold_slurm": {
                    "max_jobs_allowed": 8,
                    "sbatch": {"partition": "userpartition", "qos": "high"},
                },
            }
        )
    )
    cluster = tmp_path / "cluster.yaml"
    cluster.write_text(
        yaml.dump(
            {
                "sbatch": {
                    "partition": "clusterpartition",
                    "gres": "gpu:1",
                    "time": "4:00:00",
                    "qos": "low",
                },
            }
        )
    )
    monkeypatch.setenv("FRAGFOLDX_COLABFOLD_SLURM_DEFAULTS", str(cluster))
    p = ffparams.load_config(config_file=cfg)
    cs = p.colabfold_slurm
    assert cs.max_jobs_allowed == 8  # user control knob
    # sbatch flags deep-merged: cluster fills gres/time, user wins on partition/qos
    assert cs.sbatch == {
        "partition": "userpartition",
        "qos": "high",
        "gres": "gpu:1",
        "time": "4:00:00",
    }


# --- count_jobs: exact-match isolation, array expansion ---


def _submitter(job_name="ff-TNFa-aaaa1111"):
    sp = SlurmParameters()
    sp.job_name = job_name
    return SlurmJobSubmitter(sp)


def test_count_jobs_exact_match_isolates_runs():
    js = _submitter()
    fake = (
        " JOBID NAME\n"
        " 101 ff-TNFa-aaaa1111\n"
        " 102 ff-TNFa-bbbb2222\n"  # another run, shares prefix — must NOT count
        " 103 ff-TNFa-aaaa1111-extra\n"  # superstring — must NOT count
    )
    with mock.patch("subprocess.check_output", return_value=fake.encode()):
        assert js.count_jobs() == 1


def test_count_jobs_array_expansion():
    js = _submitter()
    fake = " JOBID NAME\n 200_[1-4] ff-TNFa-aaaa1111\n"
    with mock.patch("subprocess.check_output", return_value=fake.encode()):
        assert js.count_jobs() == 4


def test_count_jobs_empty_job_name_raises():
    js = SlurmJobSubmitter(SlurmParameters())  # job_name still ""
    with pytest.raises(ValueError):
        js.count_jobs()
