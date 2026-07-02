# FragFold-X

Predict how short peptide **fragments** of a protein bind to a full-length **target** protein, using AlphaFold2 (via ColabFold). fragfoldx slides a window along a protein sequence and, for each fragment, co-folds it against the target and scores the predicted interface. A peak in the score vs. fragment position points to a likely binding fragment.

This is a reimplementation of the original [FragFold](https://github.com/swanss/FragFold) with a rewritten codebase. A few functions are adapted from the original and credited in their docstrings (in `fragfoldx/tools/pdb_tools.py` and `fragfoldx/structure_scoring/weighted_contacts.py`).

Please cite the original FragFold paper if you use this:
A. Savinov, S. Swanson, A. E. Keating, G.-W. Li. High-throughput discovery of inhibitory protein fragments with AlphaFold. Proceedings of the National Academy of Sciences 122(6), e2322412122 (2025). doi: 10.1073/pnas.2322412122

Written by Jackson C. Halpin.

**What's different from the original:**
- multiple domains/chains as the target
- flexible fragmentation (sliding window of any stride, or fixed overlap)
- a clean Python API rather than nextflow — generate inputs, run predictions, and score structures from python (see [Python API](#python-api))

---

## Installation

**Requirements:** a CUDA-capable GPU (for predictions), git, and conda/mamba.

```bash
git clone https://github.com/jacksonh1/FragFold-X.git
cd FragFold-X
make create_environment      # creates the `fragfoldx` conda env and installs the package (editable)
conda activate fragfoldx
```

Then install **localColabFold** (required to run structure predictions) by following its [instructions](https://github.com/YoshitakaMo/localcolabfold), and point fragfoldx at it (next section).

### Point fragfoldx at ColabFold

fragfoldx needs the path to your `colabfold_batch` executable and the ColabFold data directory. Set them either way:

- **Edit `fragfoldx/executables.yaml`** (used by default in an editable install):
  ```yaml
  colabfold_batch: "/path/to/colabfold_batch"
  colabfold_data: "/path/to/colabfold_data"
  ```
- **Or set environment variables** (these override the yaml; add them to your `~/.bashrc` to make them persistent):
  ```bash
  export COLABFOLD_BATCH=/path/to/colabfold_batch
  export COLABFOLD_DATA=/path/to/colabfold_data
  ```

You can also set `colabfold_batch` / `colabfold_data` per-run in the config YAML.

---

## Usage

Run with a YAML config file (see [Configuration](#configuration)):

```bash
fragfoldx --input_params path/to/config.yaml
```

The `fragfoldx` command is installed on your PATH (equivalently, run the script directly: `python fragfoldx/scripts/run_fragfoldx.py --input_params path/to/config.yaml`). See `examples/example1/` for a complete, runnable example.

Or use it from Python:

```python
import fragfoldx
params = fragfoldx.load_config("path/to/config.yaml")
fragfoldx.fragfoldx_pipeline(params)
```

### Paths and the root directory


**TL;DR:** relative paths in the config are resolved against the config file's folder by default. Use `--root DIR` to override that base.


Paths inside the config (`fragment_source_fasta`, `target_fastas`, `msa_cache_dir`, `output_directory`, `reference_pdb`) may be absolute or relative:

- **Absolute paths** are used exactly as written.
- **Relative paths** are resolved by joining them onto a base directory. **By default the base is the directory containing the config file** — so a relative path like `data/x.fasta` in `/projects/ftsZ/config.yaml` points at `/projects/ftsZ/data/x.fasta`, no matter which directory you run `fragfoldx` from.

Config-file-relative resolution is what makes a config portable: write every path relative to the config's own folder, and the config runs the same from any working directory. Move the whole folder — to a scratch dir or another machine — and nothing inside the config needs to change.

**Example.** A config at `/projects/ftsZ/config.yaml` containing:

```yaml
fragment_source_fasta: data/fragment_protein.fasta
target_fastas:
  - data/target_domain1.fasta
output_directory: output
```

resolves to `/projects/ftsZ/data/fragment_protein.fasta`, `/projects/ftsZ/data/target_domain1.fasta`, and outputs to `/projects/ftsZ/output/` — and does so whether you run:

```bash
cd /projects/ftsZ && fragfoldx --input_params config.yaml
# ...or, from anywhere:
fragfoldx --input_params /projects/ftsZ/config.yaml
```

**Overriding the base with `--root`.** If your config's relative paths are written against some *other* directory (e.g. the config lives in a `configs/` folder but its paths are relative to a project root, or you want to redirect a run into a scratch tree), pass `--root DIR` to use `DIR` as the base instead of the config file's directory:

```bash
fragfoldx --input_params configs/run1.yaml --root /projects/ftsZ
```

**Two things the base does *not* apply to:**

- **`colabfold_batch` / `colabfold_data`.** Give these as absolute paths, or set them in the global config / environment variables.
- **The `--input_params` path itself.** You point at the config file relative to your current directory (or absolutely); the base only governs the paths *inside* it.

**Input configs vs. the saved run record.** An input config you write stays portable: `Fragfold3Params.save()` and the programmatic API preserve your path strings verbatim (relative stays relative). The `fragfold_params.yaml` that a run writes into its output directory is different — it records the **resolved, absolute** paths of that run, so it is self-contained: you can reload it later (e.g. to re-score) from anywhere with no `--root`.

From Python, the base defaults to the config file's directory too; pass `root=` to `load_config` only if you need to override it:

```python
params = fragfoldx.load_config("/projects/ftsZ/config.yaml")   # base = /projects/ftsZ
fragfoldx.fragfoldx_pipeline(params)
```

> **Generating configs programmatically?** `Fragfold3Params.save()` writes every path exactly as given, so you can build a params object with relative paths, `save()` it into the folder those paths are relative to, copy that folder elsewhere, and run it there unchanged.

> **Note (behavior change):** relative paths now resolve against the **config file's directory** by default. Previously they resolved against the current working directory (and `--root` was required for portability). Configs that use absolute paths, or that are run from the config's own directory, are unaffected; otherwise pass `--root` to reproduce the old base.

### Running on a SLURM cluster

To run the many fragment predictions in parallel as separate SLURM jobs, add `--colabfold_scheduler slurm`:

```bash
fragfoldx --input_params config.yaml --colabfold_scheduler slurm --colabfold_max_jobs_allowed 2
```

fragfoldx submits each prediction with `sbatch` and keeps the queue topped up to `--colabfold_max_jobs_allowed` until all are done. Because it monitors until completion, you'll usually submit the `fragfoldx` command itself as a job too — see `examples/example2_slurm/run_example2.sh`.

**sbatch settings** come from a JSON file. By default that's `fragfoldx/job_schedulers/colabfold_sbatch_params.json` (edit it in place, since the install is editable). To use a different file, either set the `FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE` environment variable or pass `--colabfold_sbatch_param_file /path/to/file.json`. Only one file is used — they are **not** merged — with precedence: CLI flag > environment variable > default file. Example file:

```json
{
    "--job-name": "colabfold",
    "-c": 5,
    "--partition": "pi_keating",
    "--output": "logs/colabfold_%A_%a.out",
    "--error": "logs/colabfold_%A_%a.err",
    "--gres": "gpu:1",
    "--mem": 20000
}
```

> ColabFold uses essentially a whole GPU per job, so run one job per GPU rather than several on one.

---

## Configuration

A run is fully defined by a YAML file. Example:

```yaml
# fragment: the protein you slice into windows
fragment_source_fasta: data/fragment_protein.fasta
fragment_slice_coords: [5, 34]        # region to fragment (inclusive)
fragment_length: 30                   # window width, in residues
stride: 5                             # step between windows
fragmentation_method: sliding_window  # or "overlap"

# target: what each fragment is folded against (one or more domains/chains)
target_fastas:
  - data/target_domain1.fasta
  - data/target_domain2.fasta
target_slice_coords:
  - [12, 145]
  - [50, 210]

indexing_base: "1"                    # are the coords above 1-based or 0-based?
use_fragment_msa: true
use_target_msas: true
msa_cache_dir: data/MSAs/colabfold_mmseqs
output_directory: examples/example1/output
model_weights: alphafold2_ptm
extra_colabfold_params: "--num-models 3"   # extra colabfold_batch flags, as a raw string
structure_score_params:
  contact_distance_cutoff: 4.5
  chain_group_a: ["A"]
  chain_group_b: ["B"]
  n_processes: 4
```

**Coordinates** are residue positions in the base set by `indexing_base` (`"1"` = 1-based, the default and the UniProt/PDB convention; `"0"` = 0-based). Both ends are inclusive, so `[162, 191]` is a 30-residue fragment, and `-1` means "the last residue." **Outputs use the same base you put in** — 1-based in, 1-based out.

> This applies to the fragment positions in filenames (the `fragment_start` / `fragment_end` / `fragment_center` columns and plot axes). Residue numbers inside the `contacts` column and the predicted PDB files are AlphaFold's per-chain numbering (each chain renumbered from 1), **not** source-protein coordinates.

**Paths** are absolute or relative — see [Paths and the root directory](#paths-and-the-root-directory) for how relative paths are resolved (and the `--root` option).

### Parameters

| Parameter | Type | Purpose | Default |
| --- | --- | --- | --- |
| `fragment_source_fasta` | str | Protein that is sliced into fragments | required |
| `fragment_slice_coords` | [int, int] | Region of the fragment protein to fragment (inclusive) | required |
| `target_fastas` | list[str] | Target sequence(s) each fragment is folded against. Repeat a file for homo-oligomers | required |
| `target_slice_coords` | list[[int, int]] | Slice for each target entry (same length as `target_fastas`) | required |
| `fragment_length` | int | Width of each fragment window, in residues | `30` |
| `fragmentation_method` | `sliding_window` or `overlap` | How fragments are generated (see below) | `sliding_window` |
| `stride` | int | Step between windows — **`sliding_window` mode only** | `1` |
| `overlap_length` | int | Residues each window overlaps the previous — **`overlap` mode only** | `15` |
| `indexing_base` | `"1"` or `"0"` | Whether the coords above are 1- or 0-based. Outputs honor this base | `"1"` |
| `use_fragment_msa` | bool | Use a fragment MSA (`true`) or just the single query sequence (`false`, faster/less accurate) | `true` |
| `use_target_msas` | bool | Use target MSAs (vs. single sequences) | `true` |
| `msa_cache_dir` | str | Directory where downloaded MSAs are cached and reused — see [MSA cache](#msa-cache) | `config.MSA_CACHE_DIR` |
| `output_directory` | str | Output folder | `fragfoldx_output` |
| `model_weights` | str | ColabFold model preset: `alphafold2`, `alphafold2_ptm`, `alphafold2_multimer_v1`/`_v2`/`_v3`, or `deepfold_v1` | `alphafold2_ptm` |
| `extra_colabfold_params` | str | Extra `colabfold_batch` flags appended verbatim to the command, e.g. `"--num-models 3"`. (Pair-mode is always `unpaired` — see note below) | `""` |
| `colabfold_batch` | str | Path to `colabfold_batch` (overrides global config) | `config.COLABFOLD_BATCH` |
| `colabfold_data` | str | Path to ColabFold data dir (overrides global config) | `config.COLABFOLD_DATA` |
| `structure_score_params` | dict | Scoring settings (see below) | defaults below |
| `reference_pdb` | str | Optional reference structure; currently just copied to the output (scoring not yet implemented) | `null` |
| `overwrite` | bool | Regenerate the input a3m files (only). To redo everything, delete the output dir | `true` |
| `warn_output_exists` | bool | Error out if the output dir already exists | `true` |

> ColabFold's pair-mode is always `unpaired` and cannot be changed — fragfoldx builds its inputs with `a3mcat`, which produces unpaired MSAs.

**Fragmentation methods:**
- `sliding_window` — windows of `fragment_length` residues, stepping by `stride` (e.g. `stride: 1` is every position; `stride: 5` is every 5th).
- `overlap` — consecutive windows of `fragment_length` that overlap by `overlap_length` residues.

In both cases the last window is shifted back to stay in-bounds, so it may overlap the previous one more than requested.

**`structure_score_params`:**
- `contact_distance_cutoff` (float, default `4.0`) — max atom-atom distance (Å) counted as a contact.
- `chain_group_a` / `chain_group_b` (list[str] or `null`) — the two sides of the interface to count contacts between; set both or neither. **If left null (default), the last chain (the fragment) is one side and all other chains (the target) are the other** — i.e. fragment-vs-target contacts.
- `n_processes` (int or `null`) — parallel workers for scoring (`null` = use all available / the SLURM allocation).

---

## MSA cache

Before predicting, fragfoldx needs a multiple-sequence alignment (MSA) for the fragment-source protein and for each target (unless you disable them with `use_fragment_msa` / `use_target_msas`). Each MSA is downloaded once from the ColabFold MMseqs2 server and **cached on disk** in `msa_cache_dir` (default `<install>/data/MSAs/colabfold_mmseqs`), then reused on every later run.

**How it works**
- The MSA is downloaded for the **full sequence** in each FASTA and saved as `<msa_cache_dir>/<header>.a3m`, where `<header>` is the sequence id (the first whitespace-delimited token of the FASTA header). The domain slice (`fragment_slice_coords` / `target_slice_coords`) is applied *after* loading, in memory.
- On each run, for every input fragfoldx looks for `<msa_cache_dir>/<header>.a3m`. If it exists, the download is **skipped** and the cached MSA is reused; otherwise it is downloaded and written there.

**What this buys you**
- Because the *full-length* MSA is cached and sliced in memory, one cached MSA is reused across runs that fragment **different regions** of the same protein, or that pair it with different targets — the alignment is downloaded only once.
- The cache is shared by every run that points at the same `msa_cache_dir`, so aim your configs at one shared directory to avoid re-downloading.

**Gotchas**
- ⚠️ **The cache key is the FASTA header, not the sequence.** Two different sequences that share a header id would collide on the same cache file. fragfoldx guards against this: on a cache hit it checks the cached MSA's query sequence against your input sequence and **raises an error** if they differ, telling you to rename the input header (recommended) or delete the stale cached `.a3m`. Still, give every protein a unique, descriptive header (e.g. `>ftsZ_ecoli`) to avoid the error in the first place.
- New downloads need internet access to the MMseqs2 server. On a compute node with no outbound network, pre-populate the cache from a login node first — run the pipeline once, or call `colabfold_tools.colabfold_batch_MSA_wrapper` directly (see the [Python API](#python-api)) — then the offline run reuses it.
- To force a fresh MSA, delete `<header>.a3m` (and the sibling ColabFold files) from `msa_cache_dir`.
- `msa_cache_dir` is itself resolved against the base (the config file's directory, or `--root`) like the other paths, unless an absolute path is provided.

---

## Output

```
output_directory/
├── fragfold_params.yaml     # copy of the resolved run parameters
├── structure_scores.csv     # per-structure scores (contacts, iptm, weighted contacts, …)
├── position_plot.html       # interactive plot: score vs. fragment position
├── position_plot.png        # static version of the same plot
├── input_files/             # generated a3m inputs for ColabFold
└── predictions/             # ColabFold output (.pdb structures + score .json files)
```

**`structure_scores.csv`** has one row per predicted structure. The key columns are:
- `n_contacts` — number of residue-residue contacts across the interface (within `contact_distance_cutoff`).
- `iptm` — AlphaFold's interface predicted TM-score (interface confidence, 0–1).
- `weighted_contacts` — `n_contacts × iptm`. **This is the main signal**, and what the position plots show vs. fragment position; a peak suggests a fragment that binds confidently with a large interface.

It also records `fragment_start`/`fragment_end`/`fragment_center` (the fragment's position, in your input base), `rank`, and the per-contact list. See the [coordinate note](#configuration) above for how residues are numbered in the `contacts` column.

To re-score without re-predicting, change `structure_score_params`, set `warn_output_exists: false`, and rerun. Existing predictions are kept (files already present are not regenerated); the score CSV and plots are overwritten.

---

## Python API

The pipeline is just a few steps you can also drive yourself: download MSAs → build custom a3m inputs → run ColabFold → score. MSA manipulation uses the companion [a3mcat](https://github.com/jacksonh1/a3mcat) package.

```python
from fragfoldx import colabfold_tools
import a3mcat

# 1. download an MSA from the ColabFold MMseqs2 server
#    (the output a3m is named after the FASTA header, e.g. ">ftsZ" -> ftsZ.a3m)
colabfold_tools.colabfold_batch_MSA_wrapper(
    input_file="input_sequence.fasta",
    output_dir="MSA_directory/",
)

# 2. build a custom input: slice a fragment MSA and concatenate it onto a target MSA
msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/input_sequence.a3m")
target_msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/target_msa.a3m")
fragment_msa = msa[0:100]               # first 100 residues
combined_msa = target_msa + fragment_msa
combined_msa.save("combined_msa.a3m")

# 3. run a structure prediction (pair-mode is always unpaired)
colabfold_tools.colabfold_batch_wrapper(
    input_file_or_directory="combined_msa.a3m",
    output_dir="predictions/",
    weights="alphafold2_ptm",
    extra_args="--num-models 3",   # optional: any extra colabfold_batch flags
)
```

See the [a3mcat docs](https://github.com/jacksonh1/a3mcat) for more on manipulating a3m files.

---

## Roadmap

- [ ] more scoring functions (e.g. DockQ against a reference PDB)
- [ ] renumber predicted PDBs to source-protein coordinates after prediction
- [ ] more examples and Python API documentation
