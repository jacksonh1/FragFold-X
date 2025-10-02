from pathlib import Path

from dotenv import load_dotenv
from loguru import logger
from importlib_resources import files
import yaml
import pandas as pd
import os

# Load environment variables from .env file if it exists
load_dotenv()


# If tqdm is installed, configure loguru with tqdm.write
# https://github.com/Delgan/loguru/issues/135
try:
    from tqdm import tqdm
    # Safely remove existing handlers
    try:
        logger.remove(0)
    except ValueError:
        # Handler 0 might not exist if another package already configured loguru
        pass
    logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True)
except ModuleNotFoundError:
    pass

# CHIMERAX_ALIGNMENT_SCRIPT = files("fragfold3.local_data").joinpath(
    # "align_structures_chimera.py"
# )

# ==============================================================================
# // import executables.yaml file to get the paths of external executables
# ==============================================================================
executables_file = [
    i for i in files("fragfold3").iterdir() if i.name == "executables.yaml"
][0]
if executables_file is None:
    raise FileNotFoundError("executables.yaml file not found")
else:
    with executables_file.open() as f:
        EXECUTABLES = yaml.safe_load(f)
COLABFOLD_BATCH = os.environ.get("COLABFOLD_BATCH", EXECUTABLES["colabfold_batch"])
COLABFOLD_DATA = os.environ.get("COLABFOLD_DATA", EXECUTABLES["colabfold_data"])
# COLABFOLD_BATCH = EXECUTABLES["colabfold_batch"]
# COLABFOLD_DATA = EXECUTABLES["colabfold_data"]
# MAFFT_EXECUTABLE = EXECUTABLES["mafft"]
# CD_HIT_EXECUTABLE = EXECUTABLES["cd_hit"]
# CHIMERAX_EXECUTABLE = EXECUTABLES["chimerax"]
# USALIGN_EXECUTABLE = EXECUTABLES["USalign"]


# ==============================================================================
# // main project paths
# ==============================================================================
PROJ_ROOT = Path(__file__).resolve().parents[1]
PROJ_ROOT = Path(os.environ.get("PROJ_ROOT", str(PROJ_ROOT)))
logger.info(f"fragfold3 PROJ_ROOT path is: {PROJ_ROOT}")
DATA_DIR = PROJ_ROOT / "data"
MSA_CACHE_DIR = DATA_DIR / "MSAs/colabfold_mmseqs"  # directory for storing MSA files

# ==============================================================================
# // data loading
# ==============================================================================

# ==============================================================================
# // other environment variables
# ==============================================================================

COLABFOLD_PDB_PREDICTION_FILENAME_REGEX = r"(?P<name>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"
PDB_FILENAME_REGEX = r"(?P<fragment_protein>.+)-(?P<fragment_start>\d+)to(?P<fragment_end>\d+)_vs_(?P<receptor_proteins>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"


# version
# try:
#     from fragfold3 import __version__
# except ImportError:
#     __version__ = "0.0.0"


# PDB_FILENAME_REGEX = r"(?P<fragment_protein>.+)-(?P<fragment_start>\d+)to(?P<fragment_end>\d+)_vs_(?P<receptor_proteins>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.+\.pdb"
# PDB_ALIGNED_FILENAME_REGEX = r"(?P<fragment_protein>.+)-(?P<fragment_start>\d+)to(?P<fragment_end>\d+)_vs_(?P<receptor_proteins>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d\.pdb"

