"""Programmatically generate fragfoldx input param yamls.

Pattern: keep every path RELATIVE to a project root and pass `root=` to load_config. load_config
uses root only to read the FASTAs (and resolve slice coords); the stored paths stay relative, so
save() writes portable yamls. Run a generated file later with:

    fragfoldx --input_params params/<name>.yaml --root <this dir>

Generating the yamls needs no GPU — only running the predictions does.
"""
from pathlib import Path
import fragfoldx

ROOT = Path(__file__).resolve().parent          # the root all config paths are relative to
PARAMS_DIR = ROOT / "params"
PARAMS_DIR.mkdir(exist_ok=True)

FRAGMENT_FASTA = "ftsZ_A0A140NFM6.fasta"        # relative to ROOT
MODEL_WEIGHTS = ["alphafold2_ptm", "alphafold2_multimer_v3"]


def main():
    for weights in MODEL_WEIGHTS:
        params = fragfoldx.load_config(
            root=ROOT,                          # used only to read the FASTAs
            fragment_source_fasta=FRAGMENT_FASTA,
            fragment_slice_coords=[163, 194],
            target_fastas=[FRAGMENT_FASTA],
            target_slice_coords=[[10, 317]],
            indexing_base="1",
            fragment_length=30,
            stride=1,
            model_weights=weights,
            output_directory=f"output/{weights}",   # relative to ROOT
            warn_output_exists=False,
        )
        out = PARAMS_DIR / f"ftsZ_{weights}.yaml"   # paths stay relative -> portable
        params.save(out)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
