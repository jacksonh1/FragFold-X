# fragfold3

A tool to predict how short fragments of proteins bind to full length proteins.

This project was written by Jackson C. Halpin.

The idea for the project comes from the original ([FragFold](https://github.com/swanss/FragFold)).A small number of functions were copied or adapted/modified from the original FragFold codebase ([FragFold](https://github.com/swanss/FragFold)) and these cases are noted in the function docstrings where applicable. These are in the files `fragfold3/tools/pdb_tools.py` and `fragfold3/structure_scoring/weighted_contacts.py`

some of the code here I have written for other projects over the years and adapted for this project.

# overview

This is essentially a reimplementation of FragFold, but with a few added features and rewritten source code. The main improvements are:
- ability to use multiple domains as the "receptor" for fragment binding
- different ways to generate fragments (e.g. generate a fragment every N residues rather than every single residue, or split a protein into fragments which overlap by N residues)
- completely refactored codebase that you may or may not find easier to use

The original FragFold pipeline is recreated here. The potential benefit of this code over the original is flexibility. The [a3mcat](https://github.com/jacksonh1/a3mcat) tool allows you to generate custom inputs for alphafold predictions programmatically using python. The `fragfold3` module contains wrappers for running colabfold predictions from within python, as well as utilities for scoring predicted peptide-domain structures. So you can generate custom inputs for structure prediction, execute those predictions, and analyze the resulting structures all from within python (even in a single script if you want). This is useful for large scale or custom protein structure prediction tasks where you want to generate inputs and run predictions in a single pipeline. See the [custom applications](#custom-applications) section below for more details and examples.

## upcoming features:
- [ ] add more scoring functions - particularly dockQ for the reference pdb
- [ ] add a custom example and more demos
- [ ] add more documentation for python module usage
- [ ] add explanation of wrappers
- [ ] explain root issue

# Installation

**TL;DR**
- run the following commands:
  ```bash
  git clone https://github.com/jacksonh1/FragFold3.git
  cd FragFold3
  make create_environment
  conda activate fragfold3
  ```
- edit the `./fragfold3/executables.yaml` file to set paths to colabfold installation on your system (or set environment variables as described below)
- edit the `fragfold3/job_schedulers/default_colabfold_sbatch_params.json` file if you plan to use SLURM job scheduling for colabfold predictions


## Prerequisites
- CUDA-capable GPU (recommended for ColabFold predictions)
- Git
- Conda (technically optional but required to use the makefile for installation in the following instructions)


## Install Dependencies

1. **Clone the repository:**
```bash
git clone https://github.com/jacksonh1/FragFold3.git
cd FragFold3
```

2. **Create environment, install dependencies, install package in "editable mode":**
```bash
make create_environment
```

3. **Install localColabFold (required for structure predictions):**
- see instructions at [link](https://github.com/YoshitakaMo/localcolabfold)


4. **configure paths to colabfold:**
- see [Configuring External Executable Paths](#configuring-external-executable-paths) section below for details on how to set up paths to the colabfold executables.


## Configuring External Executable Paths
You can configure the paths to the colabfold installed locally on your machine in two ways:

### 1. Using `executables.yaml`

- The file `./fragfold3/executables.yaml` contains the default paths to required executables.
- After cloning the repository, you can edit this file to set the correct paths for your system.
- Example:
  ```yaml
  colabfold_batch: "/path/to/colabfold_batch"
  colabfold_data: "/path/to/colabfold_data"
  ```
- the paths specified in this file will be used by default when running FragFold3 as long as fragfold3 is installed in editable mode (i.e. using `pip install -e .` or `make create_environment` as shown above).

### 2. Using Environment Variables
- You can override the paths in `executables.yaml` by setting environment variables before running FragFold3.
- Supported variables:
    - `COLABFOLD_BATCH`: Path to the ColabFold batch executable
    - `COLABFOLD_DATA`: Path to ColabFold data directory
    - (Other tools can be similarly overridden if supported in code)
- Example:
  ```bash
  export COLABFOLD_BATCH=/custom/path/to/colabfold_batch
  export COLABFOLD_DATA=/custom/path/to/colabfold_data
  fragfold3 --input_params config.yaml
  ```
You can set these environment variables in your shell profile (e.g., `.bashrc`, `.zshrc`, `.zshenv`, etc.) for persistent configuration if you want
for example:
```bash
echo 'export COLABFOLD_BATCH=/custom/path/to/colabfold_batch' >> ~/.bashrc
echo 'export COLABFOLD_DATA=/custom/path/to/colabfold_data' >> ~/.bashrc
```


# Usage

## TL;DR

**main command line interface** - `fragfold3`
- `fragfold3 --input_params path/to/config.yaml`
- see [Input Configuration (YAML Parameters)](#input-configuration-yaml-parameters) for details on the configuration file
- see `./examples/example1/` for an example

**if you're using the SLURM job manager** - use the `colabfold_scheduler` option to run predictions in parallel:
- `fragfold3 --input_params path/to/config.yaml --colabfold_scheduler slurm`
- see `./examples/example2_slurm/run_example2.sh` for an example sbatch script
- see [distributing colabfold predictions with a job manager](#distributing-colabfold-predictions-with-a-job-manager) for details and configuration options


## Basics

FragFold3 can be run using the main script with a YAML configuration file:

```bash
fragfold3 --input_params path/to/config.yaml
```

The `fragfold3` command will be available in your path automatically after installation. But the corresponding python script can also be run directly:

```bash
python "./fragfold3/scripts/run_fragfold3.py" --input_params path/to/config.yaml
```

Alternatively, you can use FragFold3 as a Python module in your own scripts:

```python
import fragfold3
parameters = fragfold3.load_config("path/to/config.yaml")
fragfold3.fragfold3_pipeline(parameters)
```


## distributing colabfold predictions with a job manager
Currently the only supported job manager is SLURM.

The colabfold predictions can be distributed using a job manager by specifying the `--colabfold_scheduler` argument when running fragfold3. For example, to use SLURM:

```bash
fragfold3 --input_params path/to/config.yaml --colabfold_scheduler slurm
```

The colabfold predictions will be submitted as separate jobs to the SLURM job scheduler using the `sbatch` command. The script works by monitoring the number of running jobs and submitting new jobs as others finish, until all predictions are complete. As such, you probably also want to run the fragfold3 command itself as a job in SLURM, so that it can continue monitoring until all predictions are complete. This can be done using a simple sbatch script like that used in example 2 (`./examples/example2_slurm/run_example2.sh`).

Before running, make sure to configure the sbatch parameters to match your system and requirements. You can do this in one of three ways:
1. **Using the default sbatch parameters file**: If you do not provide a custom sbatch parameters file, fragfold3 will use the default sbatch parameters file located at `fragfold3/job_schedulers/default_colabfold_sbatch_params.json`. You can modify this file directly if you installed fragfold3 in editable mode and it will use your modified version.
2. **Using the FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE environment variable**: You can set the environment variable `FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE` to point to your custom sbatch parameters file. This way, every time you run fragfold3 with the `--colabfold_scheduler slurm` option, it will use your custom default sbatch parameters file.
3. **Using a custom sbatch parameters file as a cli input**: You can create your own sbatch parameters file and provide it using the `--colabfold_sbatch_param_file` argument when running fragfold3.

Parameter precedence:
- if `--colabfold_sbatch_param_file` is provided as a cli argument, it takes precedence over the environment variable and the default file
- if the environment variable `FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE` is set, it takes precedence over the default file.<br>
parameters are not merged - only one source is used based on the precedence above

sbatch parameter file format:
The sbatch parameters file for colabfold jobs (options 1-3 above) should be in a json file format. As an example, here is the default file located in `fragfold3/job_schedulers/default_colabfold_sbatch_params.json`
```json
{
    "--job-name": "colabfold",
    "-c": 5,
    "--partition": "pi_keating",
    "--output": "logs/colabfold_%A_%a.out",
    "--error": "logs/colabfold_%A_%a.err",
    "--gres": "gpu:l40s:1",
    "--mem": 20000
}
```



alternatives or other ways to parallelize:
- use `fragfold3` to generate many parameter files, then run them in parallel using a job manager
- general note - colabfold inference uses basically 100 % of a GPU's processing power, so running multiple jobs in parallel on a single GPU is not really beneficial and will use more GPU memory. However, you can run multiple jobs in parallel on multiple GPUs.


# Configuration (YAML Parameters)

The input YAML file defines all parameters needed for a FragFold3 run. Below is a representative configuration:

```yaml
fragment_source_fasta: data/example_fragment.fasta
fragment_slice_coords: [5, 34]
receptor_fastas:
  - data/receptor_domain1.fasta
  - data/receptor_domain2.fasta
receptor_slice_coords:
  - [12, 145]
  - [50, 210]
stride: 5
fragment_length: 30
use_fragment_msa: true
use_receptor_msas: true
msa_cache_dir: data/MSAs/colabfold_mmseqs
output_directory: examples/example1/output
colabfold_batch: /opt/localcolabfold/colabfold_batch
colabfold_data: /opt/localcolabfold/data
fragmentation_method: sliding_window
indexing_base: "0"
model_weights: alphafold2_ptm
extra_colabfold_params:
  pairmode: unpaired
  num_models: 3
structure_score_params:
  contact_distance_cutoff: 4.5
  chain_group_a: ["A"]
  chain_group_b: ["B"]
  n_processes: 4
```

Coordinates are 1-based whenever `indexing_base: "1"`. **Important note** - FragFold3 converts them to 0-based regardless of what you put here. It just converts the coordinates to 0-based internally. This is reflected in the output files however. So if you provide 1-based coordinates, **the output files will be in standard pythonic 0-based coordinates**.<br>
Paths are interpreted as absolute or relative to the working directory unless `--root` is provided to the cli.

## Parameter Descriptions

| Parameter | Type | Purpose | Default |
| --- | --- | --- | --- |
| `fragment_source_fasta` | str | Fragment sequence used to seed slicing and MSA lookups | required |
| `fragment_slice_coords` | [int, int] | Start and end residue (per `indexing_base`) for the fragment | required |
| `receptor_fastas` | list[str] | One or more receptor sequences to pair with the fragment. The same fasta file can be included multiple times for homo-oligomers | required |
| `receptor_slice_coords` | list[[int, int]] | Slice coordinates for each receptor entry (same length as `receptor_fastas`) | required |
| `stride` | int | Step size between fragment windows when generating fragments | `1` |
| `fragment_length` | int | Width of each fragment window in residues | `30` |
| `use_fragment_msa` | bool | Use a fragment MSA (true) or just the query sequence (false) | `true` |
| `use_receptor_msas` | bool | Use receptor MSAs when constructing paired inputs | `true` |
| `msa_cache_dir` | str | Directory where downloaded MSAs are cached | from `config.MSA_CACHE_DIR` |
| `output_directory` | str | Root output folder for generated inputs and predictions | `fragfold3_output` |
| `colabfold_batch` | str | Path to the `colabfold_batch` executable | from `config.COLABFOLD_BATCH` |
| `colabfold_data` | str | Path to ColabFold data directory | from `config.COLABFOLD_DATA` |
| `fragmentation_method` | str (`overlap` or `sliding_window`) | Strategy for generating fragments | `sliding_window` |
| `indexing_base` | str (`"1"` or `"0"`) | Declares whether slice coordinates are 1-based or 0-based. Note - outputs are always 0-based standard pythonic coordinates | `"1"` |
| `model_weights` | str | ColabFold model preset to run (`alphafold2`, `alphafold2_ptm`, etc.) | `alphafold2_ptm` |
| `extra_colabfold_params` | dict | Additional keyword arguments forwarded to `colabfold_batch_wrapper` | `{}` |
| `structure_score_params` | dict | Settings for downstream scoring routines (see below) | defaults shown |
| `reference_pdb` | str | Optional structure used for alignment-based scoring. No scoring is implemented yet, but it does copy the reference pdb to the output directory | `null` |
| `overwrite` | bool | regenerates the input a3m files only! To regenerate all outputs, delete the output directory | `true` |
| `warn_output_exists` | bool | raise error if output directory already exists | `true` |

**Structure scoring parameters**
- `contact_distance_cutoff` (float, default `4.0`): maximum atom distance for a contact.
- `chain_group_a` and `chain_group_b` (list[str] or `null`): define residue chain groupings; provide both or neither.
- `n_processes` (int or `null`): parallel worker count for scoring utilities.

**Path handling**
- if `--root` is provided as a cli argument, all relative paths in the input YAML file are interpreted as relative to the provided root directory. If `--root` is not provided, all paths are interpreted as absolute or relative to the current working directory.


# Output

FragFold3 generates the following output structure:
```
output_directory/
├── fragfold_params.yaml      # Copy of input parameters
├── position_plot.html        # Plot of weighted contacts vs. position in fragment
├── position_plot.png         # Plot of weighted contacts vs. position in fragment
├── structure_scores.csv      # Scoring results for predicted structures
├── input_files/              # Generated MSA files for ColabFold
│   ├── fragment1_vs_receptor.a3m
│   └── fragment2_vs_receptor.a3m
└── predictions/              # ColabFold structure predictions
    ├── fragment1_vs_receptor.a3m  # MSA files for ColabFold
    ├── fragment1_vs_receptor*_unrelaxed_rank_00*_*.pdb
    ├── fragment1_vs_receptor*_scores_rank_00*_*.json
    └── ...
```


if you would like to rerun the output calculations without rerunning the structure predictions, you can just change the `structure_score_params` in the input YAML file and rerun the pipeline with warn_output_exists set to false in the yaml. It will not regenerate files that already exist in the output, but will overwrite the existing output files.

# custom applications
The fragfold3 pipeline boils down to a couple of very basic steps: 
- downloading a3m files from colabfold mmSEQS server
- manipulation of the msas to generate custom inputs for structure prediction (i.e. slicing out fragment msas and combining them with receptor msas)
  - uses the `a3mcat` tool that I created for this purpose (see [a3mcat documentation](https://github.com/jacksonh1/a3mcat))
- running colabfold structure predictions on those inputs
- scoring the resulting structures using various scoring functions


You can use the tools in fragfold3 to perform the same steps for your own custom applications if you want.

### basic example

To download MSAs from the colabfold mmSEQS server:
```python
from fragfold3 import colabfold_tools
colabfold_tools.colabfold_batch_MSA_wrapper(
    input_file="input_sequence.fasta",
    output_dir="MSA_directory/",
)
```

see [a3mcat documentation](https://github.com/jacksonh1/a3mcat) for manipulating a3m files. Here's a very basic example:
```python
import a3mcat
msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/input_sequence.a3m") 
# whatever msa you downloaded in previous step (it should be named after the header in the fasta file)

# modify msa as desired using a3mcat functions
receptor_msa = a3mcat.MSAa3m.from_a3m_file("MSA_directory/receptor_msa.a3m")
fragment_msa = msa[0:100]  # first 100 residues as fragment
combined_msa = receptor_msa + fragment_msa  # concatenate receptor and fragment msas
combined_msa.save("combined_msa.a3m")
```

to run structure predictions using ColabFold:
```python
from fragfold3 import colabfold_tools
colabfold_tools.colabfold_batch_wrapper(
    input_file_or_directory="combined_msa.a3m",
    output_dir="path/to/predictions/",
    weights="alphafold2_ptm",
    pairmode="unpaired",
)
```
