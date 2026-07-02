"""
Unit and regression test for the fragfoldx package.
"""

# Import package, test suite, and other packages as needed
import sys
from pathlib import Path

import pytest

import fragfoldx
from fragfoldx.pipeline import parameters as ffparams
from fragfoldx.pipeline import main_pipeline
from fragfoldx.tools import colabfold_tools

# ftsZ test fasta (383 residues), shipped with the repo
FASTA = Path(__file__).resolve().parents[2] / "examples/example1/ftsZ_A0A140NFM6.fasta"


def test_fragfoldx_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "fragfoldx" in sys.modules


# --- resolve_slice_coords: coords stay in the user's base, -1 resolves, bounds enforced ---


def test_resolve_slice_coords_passthrough_both_bases():
    # plain coords are returned unchanged in whatever base they came in
    assert ffparams.resolve_slice_coords((162, 193), FASTA, base=0) == (162, 193)
    assert ffparams.resolve_slice_coords((163, 194), FASTA, base=1) == (163, 194)


def test_resolve_slice_coords_end_sentinel():
    # ftsZ is 383 residues. end == -1 means "last residue" in the given base.
    assert ffparams.resolve_slice_coords((10, -1), FASTA, base=0) == (10, 382)
    assert ffparams.resolve_slice_coords((10, -1), FASTA, base=1) == (10, 383)


def test_resolve_slice_coords_start_sentinel():
    # start == -1 means "first residue" in the given base
    assert ffparams.resolve_slice_coords((-1, 50), FASTA, base=0) == (0, 50)
    assert ffparams.resolve_slice_coords((-1, 50), FASTA, base=1) == (1, 50)


@pytest.mark.parametrize(
    "coords, base",
    [
        ((0, 50), 1),       # 0 is below the first 1-based residue
        ((10, 384), 1),     # 384 is past the last 1-based residue (383)
        ((10, 383), 0),     # 383 is past the last 0-based residue (382)
        ((100, 50), 0),     # start > end
    ],
)
def test_resolve_slice_coords_out_of_range_raises(coords, base):
    with pytest.raises(ValueError):
        ffparams.resolve_slice_coords(coords, FASTA, base=base)


# --- load_config keeps the base and does not mutate / convert coords ---


def _load(base, frag, rec):
    return ffparams.load_config(
        fragment_source_fasta=str(FASTA),
        fragment_slice_coords=list(frag),
        target_fastas=[str(FASTA)],
        target_slice_coords=[list(rec)],
        indexing_base=str(base),
        output_directory="unused",
        warn_output_exists=False,
    )


def test_load_config_does_not_force_zero_based():
    p = _load(1, (163, 194), (10, 317))
    assert p.indexing_base == "1"
    assert tuple(p.fragment_slice_coords) == (163, 194)
    assert tuple(p.target_slice_coords[0]) == (10, 317)


def test_load_config_roundtrips_through_yaml(tmp_path):
    for base, frag, rec in [(0, (162, 193), (9, 316)), (1, (163, 194), (10, 317))]:
        p = _load(base, frag, rec)
        yp = tmp_path / f"params_base{base}.yaml"
        p.save(yp)
        p2 = ffparams.load_config(config_file=yp)
        assert p2.indexing_base == str(base)
        assert tuple(p2.fragment_slice_coords) == frag
        assert tuple(p2.target_slice_coords[0]) == rec


# --- pair-mode is fixed to unpaired; extra_colabfold_params is a raw arg string ---


def test_colabfold_command_is_unpaired_with_extra_args():
    cmd = colabfold_tools.colabfold_batch_wrapper(
        input_file_or_directory="in.a3m",
        output_dir="out",
        extra_args="--num-models 3",
        run=False,
    )
    assert "--pair-mode unpaired" in cmd
    assert "--num-models 3" in cmd


def test_colabfold_wrapper_rejects_pairmode():
    # pairmode is no longer a parameter — passing it must fail loudly
    with pytest.raises(TypeError):
        colabfold_tools.colabfold_batch_wrapper(
            input_file_or_directory="in.a3m", output_dir="out", pairmode="paired", run=False
        )


def test_extra_colabfold_params_accepts_string():
    p = ffparams.load_config(
        fragment_source_fasta=str(FASTA),
        fragment_slice_coords=[1, 30],
        target_fastas=[str(FASTA)],
        target_slice_coords=[[1, 100]],
        extra_colabfold_params="--num-models 3",
        warn_output_exists=False,
    )
    assert p.extra_colabfold_params == "--num-models 3"


def _params_from_dict(**overrides):
    """Build a Fragfold3Params via from_dict (no FASTA reads) with sensible required fields."""
    d = {
        "fragment_source_fasta": str(FASTA),
        "fragment_slice_coords": [1, 30],
        "target_fastas": [str(FASTA)],
        "target_slice_coords": [[1, 100]],
    }
    d.update(overrides)
    return ffparams.Fragfold3Params.from_dict(d)


def test_extra_colabfold_params_migrates_old_dict():
    # backwards-compat: an old dict is converted to the new string (keeps extra_args, drops pairmode)
    p = _params_from_dict(
        extra_colabfold_params={"pairmode": "unpaired", "extra_args": "--num-models 3"}
    )
    assert p.extra_colabfold_params == "--num-models 3"


def test_extra_colabfold_params_empty_dict_migrates_to_empty_string():
    p = _params_from_dict(extra_colabfold_params={})
    assert p.extra_colabfold_params == ""


def test_extra_colabfold_params_rejects_non_string_non_dict():
    # a genuinely wrong type still fails loudly
    with pytest.raises((TypeError, ValueError)):
        _params_from_dict(extra_colabfold_params=42)


def test_deprecated_receptor_keys_migrate_to_target():
    # backwards-compat: old `receptor_*` config keys are remapped to their `target_*` names
    d = {
        "fragment_source_fasta": str(FASTA),
        "fragment_slice_coords": [1, 30],
        "receptor_fastas": [str(FASTA)],
        "receptor_slice_coords": [[1, 100]],
        "use_receptor_msas": False,
    }
    p = ffparams.Fragfold3Params.from_dict(d)
    assert [str(f) for f in p.target_fastas] == [str(FASTA)]
    assert list(p.target_slice_coords[0]) == [1, 100]
    assert p.use_target_msas is False


def test_conflicting_receptor_and_target_keys_fail_loudly():
    # supplying both the deprecated and the new key is ambiguous — it must raise
    with pytest.raises(ValueError):
        _params_from_dict(receptor_fastas=[str(FASTA)])


def test_structure_score_params_ignores_removed_chain_groups_key():
    # old saved params include a `chain_groups` key that is no longer a field; it must be ignored
    p = _params_from_dict(
        structure_score_params={
            "chain_group_a": None,
            "chain_group_b": None,
            "chain_groups": None,
            "contact_distance_cutoff": 4.0,
            "n_processes": None,
        }
    )
    assert p.structure_score_params.contact_distance_cutoff == 4.0
    assert not hasattr(p.structure_score_params, "chain_groups")


# --- MSA cache guards against FASTA-header collisions ---


def _write_cached_a3m(path: Path, header: str, seq: str):
    # ColabFold a3m format: a "#<len>\t<count>" header line, then the query record.
    path.write_text(f"#{len(seq)}\t1\n>{header}\n{seq}\n")


def _write_fasta(path: Path, header: str, seq: str):
    path.write_text(f">{header}\n{seq}\n")


def test_assert_cached_query_matches_ok(tmp_path):
    a3m = tmp_path / "prot.a3m"
    _write_cached_a3m(a3m, "prot", "MKVLAAGGTT")
    # matching query (and case-insensitive) returns without raising
    main_pipeline._assert_cached_query_matches(a3m, "MKVLAAGGTT", "prot")
    main_pipeline._assert_cached_query_matches(a3m, "mkvlaaggtt", "prot")


def test_assert_cached_query_matches_mismatch_raises(tmp_path):
    a3m = tmp_path / "prot.a3m"
    _write_cached_a3m(a3m, "prot", "MKVLAAGGTT")
    with pytest.raises(ValueError) as exc:
        main_pipeline._assert_cached_query_matches(a3m, "QQQQQQQQQQ", "prot")
    # the message must point at the cache file and offer the `rm` fix
    msg = str(exc.value)
    assert str(a3m) in msg
    assert "rm " in msg


def test_get_colabfold_msa_reuses_matching_cache(tmp_path):
    seq = "MKVLAAGGTT"
    fasta = tmp_path / "input.fasta"
    _write_fasta(fasta, "prot", seq)
    cache = tmp_path / "cache"
    cache.mkdir()
    _write_cached_a3m(cache / "prot.a3m", "prot", seq)
    # cache hit + matching sequence: returns the cached path, no download attempted
    assert main_pipeline.get_colabfold_msa(fasta, cache) == cache / "prot.a3m"


def test_get_colabfold_msa_rejects_colliding_header(tmp_path):
    fasta = tmp_path / "input.fasta"
    _write_fasta(fasta, "prot", "MKVLAAGGTT")
    cache = tmp_path / "cache"
    cache.mkdir()
    # cached MSA shares the header but has a DIFFERENT sequence -> collision
    _write_cached_a3m(cache / "prot.a3m", "prot", "QQQQQQQQQQ")
    with pytest.raises(ValueError):
        main_pipeline.get_colabfold_msa(fasta, cache)


# --- root resolves paths for reading only; it never rewrites the stored representation ---


def _make_project(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "p.fasta").write_text(">p\nMKVLAAGGTT\n")   # 10 residues


def test_load_config_with_root_preserves_relative_paths(tmp_path):
    _make_project(tmp_path)
    p = ffparams.load_config(
        fragment_source_fasta="data/p.fasta",
        fragment_slice_coords=[1, 5],
        target_fastas=["data/p.fasta"],
        target_slice_coords=[[1, -1]],   # -1 must resolve via the root-applied read
        msa_cache_dir="msas",
        output_directory="out",
        warn_output_exists=False,
        root=tmp_path,
    )
    # stored paths are untouched (still relative); the -1 sentinel resolved against the 10-mer
    assert str(p.fragment_source_fasta) == "data/p.fasta"
    assert str(p.msa_cache_dir) == "msas"
    assert tuple(p.target_slice_coords[0]) == (1, 10)
    # save echoes them verbatim — the saved yaml stays portable
    out = tmp_path / "saved.yaml"
    p.save(out)
    text = out.read_text()
    assert "data/p.fasta" in text
    assert "msas" in text


def test_load_config_with_root_keeps_absolute_paths_absolute(tmp_path):
    # an absolute msa_cache_dir outside root must NOT be forced relative or crash (the old bug)
    _make_project(tmp_path)
    fasta = tmp_path / "data" / "p.fasta"
    abs_cache = tmp_path / "shared_cache"
    p = ffparams.load_config(
        fragment_source_fasta=str(fasta),     # absolute -> overrides root
        fragment_slice_coords=[1, 5],
        target_fastas=[str(fasta)],
        target_slice_coords=[[1, 5]],
        msa_cache_dir=str(abs_cache),         # absolute, outside root
        output_directory="out",
        warn_output_exists=False,
        root=tmp_path / "elsewhere",          # root the abs paths are NOT under
    )
    assert Path(p.msa_cache_dir) == abs_cache
    assert Path(p.fragment_source_fasta) == fasta


def test_load_config_with_root_dot_does_not_crash(tmp_path, monkeypatch):
    # `--root .` (root unresolved) used to crash convert_paths2relative; now it just works
    _make_project(tmp_path)
    monkeypatch.chdir(tmp_path)
    p = ffparams.load_config(
        fragment_source_fasta="data/p.fasta",
        fragment_slice_coords=[1, 5],
        target_fastas=["data/p.fasta"],
        target_slice_coords=[[1, 5]],
        warn_output_exists=False,
        root=Path("."),
    )
    assert str(p.fragment_source_fasta) == "data/p.fasta"
