import copy
from pathlib import Path
from typing import Any, Literal
from attrs import asdict, define, field, fields, validators
# from typing import Literal, Optional, Tuple, Any, Union
from fragfoldx import config
import yaml
from loguru import logger
import fragfoldx.tools.sequence_utils as seq_utils


def _normalize_extra_colabfold_params(value: Any) -> Any:
    """Backwards-compat for `extra_colabfold_params`.

    The field used to be a dict (keys like `pairmode` and `extra_args`); it is now a plain
    string of extra `colabfold_batch` flags. An old dict value is converted to that string by
    keeping `extra_args` only — pair-mode is always `unpaired` now, and any other keys never
    worked as colabfold flags, so they are dropped (with a warning). Strings (the new format)
    pass through unchanged; other types are left for the validator to reject.
    """
    if isinstance(value, dict):
        extra = value.get("extra_args", "")
        dropped = {k: v for k, v in value.items() if k != "extra_args"}
        if dropped:
            logger.warning(
                "`extra_colabfold_params` dict format is deprecated; converted to the string "
                f"{extra!r}. Ignored keys {dropped} — pair-mode is always 'unpaired', and other "
                "colabfold_batch flags must be passed as a raw string (e.g. '--num-models 3')."
            )
        return extra
    return value


@define
class StructureScoreParameters:
    """
    Parameters for scoring structures.
    """
    contact_distance_cutoff: float = field(default=4.0, converter=float)
    chain_group_a: list[str] | None = field(default=None)
    chain_group_b: list[str] | None = field(default=None)
    # pdb_filename_regex: str | Path = field(default=config.PDB_FILENAME_REGEX)
    # chain_groups: list[list[str]] | None = field(init=False, default=None)
    n_processes: int | None = field(default=None)

    def __attrs_post_init__(self):
        chain_groups = [self.chain_group_a, self.chain_group_b]
        if any(chain_groups) and not all(chain_groups):
            raise ValueError("If one chain group is defined, both must be defined.")
        # if chain_groups[0] is None and chain_groups[1] is None:
        #     self.chain_groups = None
        # else:
        #     self.chain_groups = chain_groups


@define
class Fragfold3Params:
    """
    """
    fragment_source_fasta: str | Path
    fragment_slice_coords: tuple[int, int] = field(factory=tuple)
    receptor_fastas: list[str | Path] = field(factory=list)
    receptor_slice_coords: list[tuple[int, int]] = field(factory=list)
    output_directory: str | Path = field(default="fragfoldx_output")
    overwrite: bool = field(default=True, converter=bool)
    warn_output_exists: bool = field(default=True, converter=bool)
    # clear_output_directory: bool = field(default=False, converter=bool)
    colabfold_batch: str | Path = field(default=config.COLABFOLD_BATCH)
    colabfold_data: str | Path = field(default=config.COLABFOLD_DATA)
    model_weights: str = field(
        default="alphafold2_ptm",
        validator=validators.in_(
            [
                "alphafold2",
                "alphafold2_ptm",
                "alphafold2_multimer_v1",
                "alphafold2_multimer_v2",
                "alphafold2_multimer_v3",
                "deepfold_v1",
            ]
        )
    )
    extra_colabfold_params: str = field(
        default="",
        converter=_normalize_extra_colabfold_params,
        validator=validators.instance_of(str),
    )
    use_fragment_msa: bool = field(default=True, converter=bool)
    use_receptor_msas: bool = field(default=True, converter=bool)
    msa_cache_dir: str | Path = field(default=config.MSA_CACHE_DIR)
    fragmentation_method: str = field(
        default="sliding_window",
        validator=validators.in_(["overlap", "sliding_window"])
    )
    fragment_length: int = field(default=30, converter=int)
    stride: int = field(default=1, converter=int)
    overlap_length: int = field(default=15, converter=int)
    indexing_base: Literal["1", "0"] = field(
        default="1",
        validator=validators.in_(["1", "0"]),
        converter=lambda x: str(x) # type: ignore
    )
    reference_pdb: str | Path | None = field(default=None)
    structure_score_params: StructureScoreParameters | None = field(default=StructureScoreParameters())

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        d = copy.deepcopy(d)
        ssp_dict = d.pop("structure_score_params", {})
        # Backwards-compat: drop structure_score_params keys that are no longer fields (e.g. the
        # removed `chain_groups`), so old saved fragfold_params.yaml files still load.
        valid = {f.name for f in fields(StructureScoreParameters)}
        unknown = {k: v for k, v in ssp_dict.items() if k not in valid}
        if unknown:
            logger.warning(
                f"Ignoring unrecognized structure_score_params keys (deprecated/removed): "
                f"{sorted(unknown)}"
            )
            ssp_dict = {k: v for k, v in ssp_dict.items() if k in valid}
        return cls(
            structure_score_params=StructureScoreParameters(**ssp_dict),
            **d
        )

    def __attrs_post_init__(self):
        """
        Post-initialization processing to convert Path-like attributes to Path objects.
        """
        self.fragment_source_fasta = Path(self.fragment_source_fasta)
        self.receptor_fastas = [Path(f) for f in self.receptor_fastas]
        self.msa_cache_dir = Path(self.msa_cache_dir)
        self.output_directory = Path(self.output_directory)
        if self.reference_pdb is not None:
            self.reference_pdb = Path(self.reference_pdb)
        self.colabfold_batch = Path(self.colabfold_batch)
        if not self.colabfold_batch.exists():
            logger.warning(f"colabfold_batch path {self.colabfold_batch} does not exist. Using {config.COLABFOLD_BATCH}.")
            self.colabfold_batch = config.COLABFOLD_BATCH
        self.colabfold_data = Path(self.colabfold_data)
        if not self.colabfold_data.exists():
            logger.warning(f"colabfold_data path {self.colabfold_data} does not exist. Using {config.COLABFOLD_DATA}.")
            self.colabfold_data = config.COLABFOLD_DATA
        # self._USalign_executable = Path(self._USalign_executable)


    def convert_paths2abs(self, root: str | Path):
        """
        Convert all Path-like attributes to absolute paths relative to the given root.
        """
        self.fragment_source_fasta = root / Path(self.fragment_source_fasta)
        self.receptor_fastas = [root / Path(f) for f in self.receptor_fastas]
        self.msa_cache_dir = root / Path(self.msa_cache_dir)
        if self.reference_pdb is not None:
            self.reference_pdb = root / Path(self.reference_pdb)
        self.output_directory = root / Path(self.output_directory)

    def convert_paths2relative(self,  root: str | Path):
        """
        Convert all Path-like attributes to relative paths from the given root.
        """
        self.fragment_source_fasta = Path(self.fragment_source_fasta).resolve().relative_to(root)
        self.receptor_fastas = [Path(f).resolve().relative_to(root) for f in self.receptor_fastas]
        self.msa_cache_dir = Path(self.msa_cache_dir).resolve().relative_to(root)
        if self.reference_pdb is not None:
            self.reference_pdb = Path(self.reference_pdb).resolve().relative_to(root)
        self.output_directory = Path(self.output_directory).resolve().relative_to(root)

    def to_writable_dict(self) -> dict[str, Any]:
        """
        Convert the Fragfold3Params object to a dictionary with writable values.
        """
        d = {}
        for k, v in asdict(self).items():
            # if k in ["_colabfold_data", "_colabfold_batch"]:
            #     # Skip private attributes
            #     continue
            if isinstance(v, Path):
                d[k] = str(v)
            elif isinstance(v, list):
                d[k] = []
                for i in v:
                    if isinstance(i, Path):
                        d[k].append(str(i))
                    elif isinstance(i, tuple):
                        d[k].append([int(p) for p in i])  # convert tuple to list
                    else:
                        d[k].append(i)
            elif isinstance(v, tuple):
                d[k] = [int(p) for p in v]
            else:
                d[k] = v
        return d

    def save(self, filename: str | Path):
        """
        Save the Fragfold3Params object to a YAML file.
        """
        filename = Path(filename)
        if not filename.suffix:
            filename = filename.with_suffix(".yaml")
        with open(filename, "w") as f:
            yaml.dump(self.to_writable_dict(), f, default_flow_style=False)


    def print_params(self):
        for k, v in asdict(self).items():
            if isinstance(v, dict):
                print(f"{k}:")
                for k2, v2 in v.items():
                    print(" ", f"{k2}:", v2)
                continue
            if isinstance(v, list):
                print(f"{k}:")
                for v2 in v:
                    print(" ", v2)
                continue
            print(f"{k}:", v)


# import typer
# app = typer.Typer()


def resolve_slice_coords(
    slice_coords: tuple[int, int], fasta_file: str | Path, base: int
) -> tuple[int, int]:
    """
    Resolve the ``-1`` "to-the-end" sentinel against the FASTA sequence length and
    validate, keeping the coordinates in the user's indexing base.

    Coordinates stay in ``base`` (0 or 1) — they are NOT converted to 0-based here. The
    conversion to 0-based happens only at the slicing boundary in the pipeline
    (``prepare_input_a3ms``). Both ends are inclusive. ``start == -1`` means "the first
    residue" (``base``) and ``end == -1`` means "the last residue" (``len(seq) - 1 + base``).

    Raises ValueError if the (resolved) coordinates fall outside ``[base, len(seq)-1+base]``.
    """
    faimporter = seq_utils.FastaImporter(fasta_file)
    seqs = faimporter.import_as_list()
    if len(seqs) != 1:
        raise ValueError(
            f"Input fasta file {fasta_file} must contain exactly one sequence."
        )
    seq_len = len(seqs[0])
    first = base
    last = seq_len - 1 + base
    start, end = slice_coords[0], slice_coords[1]
    if start == -1:
        start = first
    if end == -1:
        end = last
    if not (first <= start <= end <= last):
        raise ValueError(
            f"slice coordinates {(slice_coords[0], slice_coords[1])} (indexing base "
            f"{base}) are out of range [{first}, {last}] for {fasta_file} "
            f"(sequence length {seq_len})."
        )
    return (start, end)


def load_config(
    config_file: str | Path | None = None,
    root: Path | None = None,
    executables_file: str | Path | None = None,
    **extra_parameters
) -> Fragfold3Params:
    """
    Loads the configuration from a YAML file or user provided dictionary
    and return a Fragfold3Params object. If both `config_file` and `parameters`
    are None, raise a ValueError. If both `config_file` and `**extra_parameters` or `executables_file` are
    provided, they will be combined, with the following priority for any duplicate keys:
    `**extra_parameters` > `executables_file` > `config_file`.

    For example, if `config_file` contains:
    ```yaml
    fragment_length: 30
    ```
    and `**extra_parameters` contains:
    ```python
    fragment_length=50
    ```
    then the resulting Fragfold3Params object will have `fragment_length` set to 50.



    Parameters:
    -----------
    config_file: str | Path | None
        Path to the YAML configuration file.
    root: Path | None
        Root directory to convert relative paths to absolute paths. If None, paths are assumed to be absolute.
    executables_file: str | Path | None
        Path to a YAML file containing paths to required executables. If provided, these will be
        merged with the parameters from `config_file` and `extra_parameters`.
    **extra_parameters: dict
        Additional parameters to override those in `config_file` and `executables_file`.

    Returns:
    --------
    Fragfold3Params
        The configuration parameters as a Fragfold3Params object.
    """
    # if config_file is None and  is None:
        # raise ValueError("Either config_file or parameters must be provided.")
    if config_file is not None:
        config_file = Path(config_file)
        with open(config_file, "r") as f:
            config_dict = yaml.safe_load(f)
    else:
        config_dict = {}
    if executables_file is not None:
        executables_file = Path(executables_file)
        with open(executables_file, "r") as f:
            exec_dict = yaml.safe_load(f)
        config_dict.update(exec_dict)
    config_dict.update(extra_parameters)
    if len(config_dict) == 0:
        raise ValueError("No configuration parameters provided.")
    param_ob = Fragfold3Params.from_dict(config_dict)
    # Resolve the -1 "to-the-end" sentinel and bounds-check the slice coordinates, keeping
    # them in the user's indexing base. We do NOT convert to 0-based and do NOT mutate
    # indexing_base here: the object stays self-consistent (coords always match
    # indexing_base) and saved params round-trip. The 0-based conversion needed for slicing
    # happens locally in prepare_input_a3ms. Reading the FASTA needs absolute paths, so this
    # is wrapped in the same abs/relative path conversion as before.
    base = int(param_ob.indexing_base)
    if root is not None:
        param_ob.convert_paths2abs(root=root)
    param_ob.fragment_slice_coords = resolve_slice_coords(
        param_ob.fragment_slice_coords, param_ob.fragment_source_fasta, base
    )
    for i, receptor_slice_coord in enumerate(param_ob.receptor_slice_coords):
        param_ob.receptor_slice_coords[i] = resolve_slice_coords(
            receptor_slice_coord, param_ob.receptor_fastas[i], base
        )
    if root is not None:
        param_ob.convert_paths2relative(root=root)
    return param_ob


