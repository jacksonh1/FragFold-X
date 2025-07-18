import fragfold3.pipeline.main_pipeline as main_pipeline
import argparse
from pathlib import Path
import fragfold3.config as config
import fragfold3.pipeline.parameters as ffparams


def parse_args():
    parser = argparse.ArgumentParser(description="Run FragFold3 pipeline.")
    parser.add_argument(
        "--input_params",
        type=str,
        required=True,
        help="Path to the input parameters YAML file.",
    )
    parser.add_argument(
        '--root',
        default=None,
        help=f"Root directory. Path that all relative paths in the input parameters file are relative to. If not provided, uses absolute paths."
    )
    parser.add_argument(
        '--scheduler',
        type=str,
        choices=['slurm', 'none'],
        default='slurm',
        help="Job scheduler to use for distributing predictions on HPC. Use 'none' to run without a scheduler."
    )
    parser.add_argument(
        '--max_jobs_allowed',
        type=int,
        default=32,
        help="Maximum number of jobs allowed to run simultaneously. Only used if --scheduler is set to 'slurm'."
    )
    # parser.add_argument(
    #     '--version',
    #     action='version',
    #     version=f'FragFold3 version {config.__version__}',
    #     help="Show the version of FragFold3."
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
    if args.scheduler == 'slurm':
        main_pipeline.fragfold3_pipeline_scheduler(
            params=param_ob,
            root=args.root,
            max_jobs_allowed=args.max_jobs_allowed,
        )
    elif args.scheduler == 'none':
        main_pipeline.fragfold3_pipeline(
            params=param_ob,
            root=args.root,
        )


if __name__ == "__main__":
    main()