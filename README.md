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
- completely refactored codebase with a more modular design

The original FragFold pipeline is recreated here, however the main utility of this code over the original is its flexibility. The [a3mcat](https://github.com/jacksonh1/a3mcat) tool allows you to generate custom inputs for alphafold predictions programmatically using python. The `fragfold3` module contains wrappers for running colabfold predictions from within python, as well as utilities for scoring predicted peptide-domain structures. So you can generate custom inputs for structure prediction, execute those predictions, and analyze the resulting structures all from within python (even in a single script if you want). This is useful for large scale protein structure prediction tasks where you want to generate inputs and run predictions in a single pipeline.

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


#### Configuring External Executable Paths
You can configure the paths to the colabfold installed locally on your machine in two ways:

##### 1. Using `executables.yaml`

- The file `./fragfold3/executables.yaml` contains the default paths to required executables.
- After cloning the repository, you can edit this file to set the correct paths for your system.
- Example:
  ```yaml
  colabfold_batch: "/path/to/colabfold_batch"
  colabfold_data: "/path/to/colabfold_data"
  ```
- the paths specified in this file will be used by default when running FragFold3 as long as fragfold3 is installed in editable mode (i.e. using `pip install -e .` or `make create_environment` as shown above).

##### 2. Using Environment Variables
- You can override the paths in `executables.yaml` by setting environment variables before running FragFold3.
- Supported variables:
    - `COLABFOLD_BATCH`: Path to the ColabFold batch executable
    - `COLABFOLD_DATA`: Path to ColabFold data directory
    - (Other tools can be similarly overridden if supported in code)
- Example:
  ```bash
  export COLABFOLD_BATCH=/custom/path/to/colabfold_batch
  export COLABFOLD_DATA=/custom/path/to/colabfold_data
  run_fragfold3 --input_params config.yaml
  ```
You can set these environment variables in your shell profile (e.g., `.bashrc`, `.zshrc`, `.zshenv`, etc.) for persistent configuration if you want


## Usage - command line

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

### distributing colabfold predictions with a job manager
Currently the only supported job manager is SLURM.

[explain sbatch parameters file]

[explain where in the source code you could add support for other job managers]

[explain generally how it works]

```bash
```

alternatives or other ways to parallelize:
- use `fragfold3` to generate many parameter files, then run them in parallel using a job manager
- generate many parameter files
- general note - colabfold inference uses basically 100 % of a GPU's processing power, so running multiple jobs in parallel on a single GPU is not really beneficial and will use more GPU memory. However, you can run multiple jobs in parallel on multiple GPUs.
- use the `fragfold3.main_pipeline....` function to generate a list of colabfold commands that need to be run, then run them in parallel using a job manager or other parallelization method.



## Input Configuration (YAML Parameters)

The input YAML file defines all parameters needed for a FragFold3 run. Below is a complete example with explanations:

```yaml
```

### Parameter Descriptions


### Output Structure

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

## Usage - Python module
Alternatively, you can use FragFold3 as a Python module in your own scripts:

```python
import fragfold3
parameters = fragfold3.load_config("path/to/config.yaml")
fragfold3.fragfold3_pipeline(parameters)
```

### custom applications
To generate custom a3m inputs for structure prediction, you can use the `a3mcat` tool (see [a3mcat documentation](https://github.com/jacksonh1/a3mcat))

To download MSAs from the colabfold mmSEQS server:
```python
from fragfold3 import ...
```

to run structure predictions using ColabFold:
```python
from fragfold3 import ...
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
