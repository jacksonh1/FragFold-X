# fragfold3

Predict how short peptide **fragments** of a protein bind to a full-length **receptor** protein, using AlphaFold2 (via ColabFold). fragfold3 slides a window along a protein sequence and, for each fragment, co-folds it against the receptor and scores the predicted interface. A peak in the score vs. fragment position points to a likely binding fragment.

This is a reimplementation of the original [FragFold](https://github.com/swanss/FragFold) with a rewritten, more flexible codebase. A few functions are adapted from the original and credited in their docstrings (in `fragfold3/tools/pdb_tools.py` and `fragfold3/structure_scoring/weighted_contacts.py`). Written by Jackson C. Halpin.

**What's different from the original:**
- multiple domains/chains as the receptor
- flexible fragmentation (sliding window of any stride, or fixed overlap)
- a clean Python API — generate inputs, run predictions, and score structures from a single script (see [Python API](#python-api))

---

## Installation

**Requirements:** a CUDA-capable GPU (for predictions), git, and conda/mamba.

```bash
git clone https://github.com/jacksonh1/FragFold3.git
cd FragFold3
make create_environment      # creates the `fragfold3` conda env and installs the package (editable)
conda activate fragfold3
```

Then install **localColabFold** (required to run structure predictions) by following its [instructions](https://github.com/YoshitakaMo/localcolabfold), and point fragfold3 at it (next section).

### Point fragfold3 at ColabFold

fragfold3 needs the path to your `colabfold_batch` executable and the ColabFold data directory. Set them either way:

- **Edit `fragfold3/executables.yaml`** (used by default in an editable install):
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
fragfold3 --input_params path/to/config.yaml
```

The `fragfold3` command is installed on your PATH (equivalently, run the script directly: `python fragfold3/scripts/run_fragfold3.py --input_params path/to/config.yaml`). See `examples/example1/` for a complete, runnable example.

Or use it from Python:

```python
import fragfold3
params = fragfold3.load_config("path/to/config.yaml")
fragfold3.fragfold3_pipeline(params)
```

### Running on a SLURM cluster

To run the many fragment predictions in parallel as separate SLURM jobs, add `--colabfold_scheduler slurm`:

```bash
fragfold3 --input_params config.yaml --colabfold_scheduler slurm --colabfold_max_jobs_allowed 2
```

fragfold3 submits each prediction with `sbatch` and keeps the queue topped up to `--colabfold_max_jobs_allowed` until all are done. Because it monitors until completion, you'll usually submit the `fragfold3` command itself as a job too — see `examples/example2_slurm/run_example2.sh`.

**sbatch settings** come from a JSON file. By default that's `fragfold3/job_schedulers/colabfold_sbatch_params.json` (edit it in place, since the install is editable). To use a different file, either set the `FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE` environment variable or pass `--colabfold_sbatch_param_file /path/to/file.json`. Only one file is used — they are **not** merged — with precedence: CLI flag > environment variable > default file. Example file:

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

# receptor: what each fragment is folded against (one or more domains/chains)
receptor_fastas:
  - data/receptor_domain1.fasta
  - data/receptor_domain2.fasta
receptor_slice_coords:
  - [12, 145]
  - [50, 210]

indexing_base: "1"                    # are the coords above 1-based or 0-based?
use_fragment_msa: true
use_receptor_msas: true
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

**Paths** are absolute or relative to the current directory, unless you pass `--root DIR` on the CLI, in which case relative paths are resolved against `DIR`.

### Parameters

| Parameter | Type | Purpose | Default |
| --- | --- | --- | --- |
| `fragment_source_fasta` | str | Protein that is sliced into fragments | required |
| `fragment_slice_coords` | [int, int] | Region of the fragment protein to fragment (inclusive) | required |
| `receptor_fastas` | list[str] | Receptor sequence(s) each fragment is folded against. Repeat a file for homo-oligomers | required |
| `receptor_slice_coords` | list[[int, int]] | Slice for each receptor entry (same length as `receptor_fastas`) | required |
| `fragment_length` | int | Width of each fragment window, in residues | `30` |
| `fragmentation_method` | `sliding_window` or `overlap` | How fragments are generated (see below) | `sliding_window` |
| `stride` | int | Step between windows — **`sliding_window` mode only** | `1` |
| `overlap_length` | int | Residues each window overlaps the previous — **`overlap` mode only** | `15` |
| `indexing_base` | `"1"` or `"0"` | Whether the coords above are 1- or 0-based. Outputs honor this base | `"1"` |
| `use_fragment_msa` | bool | Use a fragment MSA (`true`) or just the single query sequence (`false`, faster/less accurate) | `true` |
| `use_receptor_msas` | bool | Use receptor MSAs (vs. single sequences) | `true` |
| `msa_cache_dir` | str | Where downloaded MSAs are cached (reused across runs; keyed by FASTA header) | `config.MSA_CACHE_DIR` |
| `output_directory` | str | Output folder | `fragfold3_output` |
| `model_weights` | str | ColabFold model preset: `alphafold2`, `alphafold2_ptm`, `alphafold2_multimer_v1`/`_v2`/`_v3`, or `deepfold_v1` | `alphafold2_ptm` |
| `extra_colabfold_params` | str | Extra `colabfold_batch` flags appended verbatim to the command, e.g. `"--num-models 3"`. (Pair-mode is always `unpaired` — see note below) | `""` |
| `colabfold_batch` | str | Path to `colabfold_batch` (overrides global config) | `config.COLABFOLD_BATCH` |
| `colabfold_data` | str | Path to ColabFold data dir (overrides global config) | `config.COLABFOLD_DATA` |
| `structure_score_params` | dict | Scoring settings (see below) | defaults below |
| `reference_pdb` | str | Optional reference structure; currently just copied to the output (scoring not yet implemented) | `null` |
| `overwrite` | bool | Regenerate the input a3m files (only). To redo everything, delete the output dir | `true` |
| `warn_output_exists` | bool | Error out if the output dir already exists | `true` |

> ColabFold's pair-mode is always `unpaired` and cannot be changed — fragfold3 builds its inputs with `a3mcat`, which produces unpaired MSAs.

**Fragmentation methods:**
- `sliding_window` — windows of `fragment_length` residues, stepping by `stride` (e.g. `stride: 1` is every position; `stride: 5` is every 5th).
- `overlap` — consecutive windows of `fragment_length` that overlap by `overlap_length` residues.

In both cases the last window is shifted back to stay in-bounds, so it may overlap the previous one more than requested.

**`structure_score_params`:**
- `contact_distance_cutoff` (float, default `4.0`) — max atom-atom distance (Å) counted as a contact.
- `chain_group_a` / `chain_group_b` (list[str] or `null`) — the two sides of the interface to count contacts between; set both or neither. **If left null (default), the last chain (the fragment) is one side and all other chains (the receptor) are the other** — i.e. fragment-vs-receptor contacts.
- `n_processes` (int or `null`) — parallel workers for scoring (`null` = use all available / the SLURM allocation).

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
from fragfold3 import colabfold_tools
import a3mcat

# 1. download an MSA from the ColabFold MMseqs2 server
#    (the output a3m is named after the FASTA header, e.g. ">ftsZ" -> ftsZ.a3m)
colabfold_tools.colabfold_batch_MSA_wrapper(
    input_file="input_sequence.fasta",
    output_dir="MSA_directory/",
)

# 2. build a custom input: slice a fragment MSA and concatenate it onto a receptor MSA
msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/input_sequence.a3m")
receptor_msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/receptor_msa.a3m")
fragment_msa = msa[0:100]               # first 100 residues
combined_msa = receptor_msa + fragment_msa
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
