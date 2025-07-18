# fragfold3

A tool to predict how short fragments of proteins bind to full length proteins

This project was written by Jackson C. Halpin.

directory structure was based on CCDS template (https://cookiecutter-data-science.drivendata.org/)

A small number of functions were adapted from the original FragFold codebase ([FragFold](https://github.com/swanss/FragFold)) and these cases are noted in the function docstrings where applicable. The only cases of this are in the files `fragfold3/tools/pdb_tools.py` and `fragfold3/structure_scoring/weighted_contacts.py`


## TODO/NOTES
- I am unhappy with the table/yaml organization I think? What is good is that they are decoupled from the main pipeline

- you can change to 0-based indexing
or you could include 1-based or 0based as a parameter in the yaml file.

- get glycosylation sites from a csv or annotation file, rather than structures
- colabfold predicted structures as reference pdbs?
    - Issue with disordered loops, etc. i.e. the colabfold model reference is not a reliable reference at those positions.
    - Could do a plddt normalized score. This would be a fragment to full length prediction DockQ score, normalized by full length prediction plddt score.

GPU multiprocessing and multiprocessing for structure scoring should be handled separately.

## Installation

### Prerequisites
- CUDA-capable GPU (recommended for ColabFold predictions)
- Git
- Conda (technically optional but required to use the makefile for installation in the following instructions)

### Install Dependencies

1. **Clone the repository:**
```bash
git clone <repository-url>
cd FragFold3
```

2. **Create environment, install dependencies, install package in "editable mode":**
```bash
make create_environment
```

3. **Install localColabFold (required for structure predictions):**
- see instructions at [link]()


4. **configure paths to colabfold:**
- see [Configuring External Executable Paths](#configuring-external-executable-paths) section below for details on how to set up paths to the colabfold executables.


#### Configuring External Executable Paths
You can configure the paths to the colabfold installed locally on your machine in two ways:

##### 1. Using `executables.yaml`

- The file `fragfold3/executables.yaml` contains the default paths to required executables.
- After cloning the repository, you can edit this file to set the correct paths for your system.
- Example:
  ```yaml
  colabfold_batch: "/path/to/colabfold_batch"
  colabfold_data: "/path/to/colabfold_data"
  ```

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


## Usage

### Basic Usage

FragFold3 can be run using the main script with a YAML configuration file:

```bash
run_fragfold3 --input_params path/to/config.yaml
```

### Input Configuration (YAML Parameters)

The input YAML file defines all parameters needed for a FragFold3 run. Below is a complete example with explanations:

```yaml
```

#### Parameter Descriptions


### Example Workflows


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

for notebooks, consider the following naming convention:
![image](./docs/image.png)

