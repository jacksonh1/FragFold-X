# fragfold3

A tool to predict how short fragments of proteins bind to full length proteins.

This project was written by Jackson C. Halpin.

directory structure was based on CCDS template (https://cookiecutter-data-science.drivendata.org/)

A small number of functions were copied or adapted/modified from the original FragFold codebase ([FragFold](https://github.com/swanss/FragFold)) and these cases are noted in the function docstrings where applicable. The only cases of this are in the files `fragfold3/tools/pdb_tools.py` and `fragfold3/structure_scoring/weighted_contacts.py`

some of the code here I have written for other projects over the years and adapted for this project.

# overview

This is essentially a reimplementation of FragFold, but with a few improvements and rewritten source code. The main improvements are:
- ability to use multiple domains as the "receptor" for fragment binding
- different ways to generate fragments (e.g. generate a fragment every N residues rather than every single residue, or split a protein into fragments which overlap by N residues)
- completely refactored codebase that you may or may not find easier to use

The original FragFold pipeline is recreated here. The main benefit of this code over the original is flexibility. The [a3mcat](https://github.com/jacksonh1/a3mcat) tool allows you to generate custom inputs for alphafold predictions programmatically using python. The `fragfold3` module contains wrappers for running colabfold predictions from within python, as well as utilities for scoring predicted peptide-domain structures. So you can generate custom inputs for structure prediction, execute those predictions, and analyze the resulting structures all from within python (even in a single script if you want). This is useful for large scale protein structure prediction tasks where you want to generate inputs and run predictions in a single pipeline.

## upcoming features:
- [ ] add more scoring functions
- [ ] add a custom example and more demos
- [ ] add more documentation for python module usage
- [ ] add explanation of wrappers
- [ ] explain root issue

## Installation

### Prerequisites
- CUDA-capable GPU (recommended for ColabFold predictions)
- Git
- Conda (technically optional but required to use the makefile for installation in the following instructions)

### Install Dependencies

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


### Configuring External Executable Paths
You can configure the paths to the colabfold installed locally on your machine in two ways:

#### 1. Using `executables.yaml`

- The file `./fragfold3/executables.yaml` contains the default paths to required executables.
- After cloning the repository, you can edit this file to set the correct paths for your system.
- Example:
  ```yaml
  colabfold_batch: "/path/to/colabfold_batch"
  colabfold_data: "/path/to/colabfold_data"
  ```
- the paths specified in this file will be used by default when running FragFold3 as long as fragfold3 is installed in editable mode (i.e. using `pip install -e .` or `make create_environment` as shown above).

#### 2. Using Environment Variables
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


## Usage

### TL;DR

**main command line interface** - `fragfold3`
- `fragfold3 --input_params path/to/config.yaml`
- see [Input Configuration (YAML Parameters)](#input-configuration-yaml-parameters) for details on the configuration file
- see `./examples/example1/` for an example

**if you're using the SLURM job manager** - use the `colabfold_scheduler` option to run predictions in parallel:
- `fragfold3 --input_params path/to/config.yaml --colabfold_scheduler slurm`
- see `./examples/example2_slurm/run_example2.sh` for an example sbatch script
- see [distributing colabfold predictions with a job manager](#distributing-colabfold-predictions-with-a-job-manager) for details and configuration options


### Basics

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


### distributing colabfold predictions with a job manager
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
- if the environment variable `FRAGFOLD3_COLABFOLD_SBATCH_PARAM_FILE` is set, it takes precedence over the default file
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


## Configuration (YAML Parameters)

The input YAML file defines all parameters needed for a FragFold3 run. Below is a complete example with explanations:

```yaml
```

### Parameter Descriptions


## Output

FragFold3 generates the following output structure:
```
output_directory/
├── fragfold_params.yaml      # Copy of input parameters
├── input_files/              # Generated MSA files for ColabFold
│   ├── fragment1_vs_receptor.a3m
│   └── fragment2_vs_receptor.a3m
└── predictions/              # ColabFold structure predictions
    ├── fragment1_vs_receptor*_unrelaxed_rank_001_*.pdb
    └── *.json
```


## custom applications
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

## Project Organization

```
├── LICENSE                <- Open-source license if one is chosen
├── Makefile               <- Makefile with convenience commands like `make create_environment` or `make requirements`
├── README.md              <- The top-level README for developers using this project.
├── data
│   ├── processed          <- The final, canonical data sets for modeling.
│   └── raw                <- The original, immutable data dump.
│
├── notebooks              <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the date in yyyy-mm-dd format, and a short `_` delimited description, e.g.
│                         `01-2025-05-13-data_exploration.ipynb`. "DE" equals data exploration.
│
├── pyproject.toml         <- Project configuration file with package metadata for
│                         lir_proteome_screen_pssm and configuration for tools like black
│
├── ref_materials          <- powerpoints + old code and such to draw from
│
├── reports                <- Generated analysis as HTML, PDF, LaTeX, ppt, etc.
│   └── figures            <- Generated graphics and figures to be used in reporting
│
├── requirements.txt       <- The requirements file for reproducing the analysis environment, e.g.
│                           generated with `pip freeze > requirements.txt`
│
├── processing_scripts     <- scripts for processing data. This is where the tools in `lir_proteome_screen_pssm` are actually applied.
│
└── fragfold3   <- Source code for use in this project.
    │
    ├── __init__.py        <- Makes fragfold3 a Python module
    │
    ├── config.py          <- Store useful global variables and configuration
    │
    └── ...                <- ...
```
