import fragfoldx.pipeline.main_pipeline as main_pipeline
import argparse
from pathlib import Path
import fragfoldx.config as config
import fragfoldx.pipeline.parameters as ffparams
from loguru import logger
import fragfoldx.job_schedulers.slurm_job_submitter as slurm_job_submitter
import json


def parse_args():
    parser = argparse.ArgumentParser(description="Run FragFold-X pipeline.")
    parser.add_argument(
        "--input_params",
        type=str,
        required=True,
        help="Path to the input parameters YAML file.",
    )
    parser.add_argument(
        '--root',
        default=None,
        help="Override the base directory that relative paths in the input parameters file are resolved against. If not provided, the base defaults to the directory containing the --input_params file (absolute paths are used as-is)."
    )
    parser.add_argument(
        '--colabfold_scheduler',
        type=str,
        choices=['slurm', 'none'],
        default='none',
        help="Job scheduler to use for distributing colabfold predictions on HPC. Use 'none' to run without a scheduler."
    )
    parser.add_argument(
        '--colabfold_max_jobs_allowed',
        type=int,
        default=32,
        help="Maximum number of colabfold jobs allowed to run simultaneously. Only used if --colabfold_scheduler is set to 'slurm'."
    )
    parser.add_argument(
        '--colabfold_sbatch_param_file',
        type=str,
        default=config.COLABFOLD_SBATCH_PARAM_FILE,
        help=f"Path to the colabfold sbatch parameters JSON file. If not provided, uses {config.COLABFOLD_SBATCH_PARAM_FILE}"
    )
    # parser.add_argument(
    #     '--version',
    #     action='version',
    #     version=f'FragFold-X version {config.__version__}',
    #     help="Show the version of FragFold-X."
    # )
    args = parser.parse_args()
    if args.root is not None:
        args.root = Path(args.root)
    return args


def main():
    args = parse_args()
    input_params_path = Path(args.input_params)

    if not input_params_path.exists():
        raise FileNotFoundError(f"Input parameters file {input_params_path} does not exist.")

    param_ob = ffparams.load_config(config_file=input_params_path, root=args.root)
    if args.colabfold_scheduler == 'slurm':
        if not Path(args.colabfold_sbatch_param_file).exists():
            raise FileNotFoundError(f"SLURM sbatch parameters file {args.colabfold_sbatch_param_file} does not exist.")
        logger.info(f"Using SLURM scheduler with sbatch parameters file: {args.colabfold_sbatch_param_file}")
        colab_slurm_submitter = slurm_job_submitter.SlurmJobSubmitter(
            sbatch_param_file=args.colabfold_sbatch_param_file
        )
        logger.info(f"SLURM sbatch parameters: {colab_slurm_submitter._create_sbatch_param_str()}")
        # load_config already recorded the base on param_ob.base_dir (from --root or the config
        # file's directory), so the pipeline resolves paths from that; no need to pass root again.
        main_pipeline.fragfoldx_pipeline_scheduler(
            params=param_ob,
            max_jobs_allowed=args.colabfold_max_jobs_allowed,
            job_submitter=colab_slurm_submitter,
        )
    elif args.colabfold_scheduler == 'none':
        main_pipeline.fragfoldx_pipeline(
            params=param_ob,
        )


if __name__ == "__main__":
    main()