import fragfoldx.pipeline.main_pipeline as main_pipeline
import argparse
from pathlib import Path
import fragfoldx.pipeline.parameters as ffparams
from loguru import logger


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
        default=None,
        help="Override colabfold_slurm.max_jobs_allowed (max simultaneous colabfold jobs). If omitted, the value from the config / cluster defaults is used. Only used with --colabfold_scheduler slurm."
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
        # CLI override is applied to the params object BEFORE the pipeline writes the run record,
        # so the effective throttle value is recorded in fragfold_params.yaml.
        if args.colabfold_max_jobs_allowed is not None:
            param_ob.colabfold_slurm.max_jobs_allowed = args.colabfold_max_jobs_allowed
        logger.info(
            f"Using SLURM scheduler; colabfold job name {param_ob.colabfold_slurm.job_name}, "
            f"max_jobs_allowed={param_ob.colabfold_slurm.max_jobs_allowed}"
        )
        # load_config already recorded the base on param_ob.base_dir (from --root or the config
        # file's directory), so the pipeline resolves paths from that; no need to pass root again.
        main_pipeline.fragfoldx_pipeline_scheduler(params=param_ob)
    elif args.colabfold_scheduler == 'none':
        main_pipeline.fragfoldx_pipeline(
            params=param_ob,
        )


if __name__ == "__main__":
    main()