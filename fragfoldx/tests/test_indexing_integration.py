"""
Integration checks that residue indexing is correct at every step, for BOTH 0-based and
1-based input. These verify actual SEQUENCE CONTENT (a3m queries, committed PDB chains, and
the scorer CSV) rather than just filenames, which is what catches real off-by-one bugs.

They depend on example fixtures that may be absent in a bare clone (the cached ftsZ MSA lives
under the gitignored ``data/MSAs/``), so each test skips when its fixtures are missing. The
pure-logic guarantees (``resolve_slice_coords``, base preservation, YAML round-trip) live in
``test_fragfoldx.py`` and need no fixtures.
"""

from pathlib import Path

import pytest

import fragfoldx.config as config
import fragfoldx.pipeline.main_pipeline as main_pipeline
import fragfoldx.pipeline.parameters as ffparams
import fragfoldx.pipeline.result_summary as result_summary
import fragfoldx.tools.pdb_tools as pdb_tools

REPO = Path(__file__).resolve().parents[2]
FASTA = REPO / "examples/example1/ftsZ_A0A140NFM6.fasta"
MSA_CACHE = config.MSA_CACHE_DIR / "ftsZ.a3m"
PREDICTIONS = REPO / "examples/example1/output/predictions"

# ftsZ source sequence (0-indexed). example1 is a 0-based config; its 1-based equivalent is
# the same physical region with every coordinate shifted +1.
S = (
    "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))
    if FASTA.exists()
    else ""
)

needs_msa = pytest.mark.skipif(
    not MSA_CACHE.exists(), reason=f"cached MSA {MSA_CACHE} not present"
)
needs_pdbs = pytest.mark.skipif(
    not (PREDICTIONS.exists() and list(PREDICTIONS.glob("*.pdb"))),
    reason="committed example prediction PDBs not present",
)


def _parse_name(stem):
    """ftsZ-162to191_vs_ftsZ -> (162, 191)"""
    a, b = stem.split("_vs_")[0].split("-")[-1].split("to")
    return int(a), int(b)


def _generate(base, frag, rec, out_dir):
    """Generate the input a3ms for one config; return (filenames, windows).

    windows is a list of (start, end, receptor_query, fragment_query) where start/end are the
    display-base residue positions parsed from the filename.
    """
    import a3mcat

    params = ffparams.load_config(
        fragment_source_fasta=str(FASTA),
        fragment_slice_coords=list(frag),
        receptor_fastas=[str(FASTA)],
        receptor_slice_coords=[list(rec)],
        indexing_base=str(base),
        fragment_length=30,
        stride=1,
        output_directory=str(out_dir),
        warn_output_exists=False,
    )
    a3m_files, _ = main_pipeline.prepare_input_a3ms(params)
    rec_len = rec[1] - rec[0] + 1  # inclusive, base-independent
    names, windows = [], []
    for f in sorted(a3m_files):
        names.append(f.name)
        start, end = _parse_name(f.stem)
        query = a3mcat.MSAa3m.from_a3m_file(f).query.seq_str
        windows.append((start, end, query[:rec_len], query[rec_len:]))
    return names, windows


@needs_msa
@pytest.mark.parametrize(
    "base, frag, rec",
    [(0, (162, 193), (9, 316)), (1, (163, 194), (10, 317))],
)
def test_input_a3m_sequences_match_source(tmp_path, base, frag, rec):
    """The sliced fragment/receptor queries in each generated a3m are the right residues."""
    names, windows = _generate(base, frag, rec, tmp_path)
    assert windows, "no input a3ms were generated"
    for start, end, receptor_q, fragment_q in windows:
        # display-base residue N is python index N - base into S
        assert fragment_q == S[start - base : end - base + 1]
        assert receptor_q == S[rec[0] - base : rec[1] - base + 1]


@needs_msa
def test_cross_base_invariants(tmp_path):
    """0-based and 1-based runs over the same region: identical sequences, labels shifted +1."""
    names0, w0 = _generate(0, (162, 193), (9, 316), tmp_path / "b0")
    names1, w1 = _generate(1, (163, 194), (10, 317), tmp_path / "b1")
    # same biology — only the labels differ
    assert sorted(fq for *_, fq in w0) == sorted(fq for *_, fq in w1)
    # 1-based filenames are exactly the 0-based names with start/end shifted +1
    shifted = set()
    for n in names1:
        s, e = _parse_name(Path(n).stem)
        shifted.add(n.replace(f"-{s}to{e}_", f"-{s - 1}to{e - 1}_"))
    assert shifted == set(names0)


@needs_pdbs
def test_committed_pdb_chains_match_source():
    """Real predicted-structure chains match the residues their filename claims (0-based example)."""
    pdbs = sorted(PREDICTIONS.glob("ftsZ-*_unrelaxed_rank_001_*.pdb"))
    assert pdbs, "no committed rank_001 PDBs found"
    for pdb in pdbs:
        start, end = _parse_name(pdb.name.split("_unrelaxed")[0])
        seqs = pdb_tools.extract_sequences_from_pdb(pdb)
        chains = sorted(seqs)
        # a3m is built receptor + fragment, so chain A = receptor, last chain = fragment
        assert seqs[chains[-1]] == S[start : end + 1]
        assert seqs[chains[0]] == S[9:317]


@needs_pdbs
def test_scorer_csv_columns_match_filenames(tmp_path):
    """structure_scores.csv fragment_start/end inherit the display base from the filenames."""
    import pandas as pd

    csv = tmp_path / "scores.csv"
    result_summary.score_pdb_files_multiprocessing(
        input_directory=PREDICTIONS, output_file=csv, n_processes=1
    )
    df = pd.read_csv(csv)
    assert len(df) > 0
    for _, r in df.iterrows():
        stem = Path(r["pdb_file"]).name.split("_unrelaxed")[0]
        fs, fe = _parse_name(stem)
        assert int(r["fragment_start"]) == fs
        assert int(r["fragment_end"]) == fe
        assert float(r["fragment_center"]) == (fs + fe) / 2


@needs_pdbs
def test_scorer_records_params_file(tmp_path):
    """The CSV's `fragfold_processing_params` column is populated when fragfold_params.yaml exists.

    Guards the pipeline ordering fix: params are now saved (right after setup) *before* scoring,
    so score_pdb finds `<output>/fragfold_params.yaml` and records it, instead of warning that the
    file is missing and leaving the column empty.
    """
    import pandas as pd

    params_yaml = PREDICTIONS.parent / "fragfold_params.yaml"
    if not params_yaml.exists():
        pytest.skip(f"{params_yaml} not present")
    csv = tmp_path / "scores.csv"
    result_summary.score_pdb_files_multiprocessing(
        input_directory=PREDICTIONS, output_file=csv, n_processes=1
    )
    df = pd.read_csv(csv)
    assert len(df) > 0
    assert df["fragfold_processing_params"].notna().all()
    for val in df["fragfold_processing_params"]:
        assert Path(val).name == "fragfold_params.yaml"
        assert Path(val).exists()
