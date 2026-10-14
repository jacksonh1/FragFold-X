import subprocess
import time
import os
import re
from typing import TYPE_CHECKING
from loguru import logger

if TYPE_CHECKING:
    from fragfoldx.pipeline.parameters import SlurmParameters


class SlurmJobSubmitter:
    """Submit and throttle colabfold prediction jobs on SLURM from a typed SlurmParameters.

    The sbatch flags come entirely from `slurm_params` (typed fields + `extra`). Throttling
    counts only this run's jobs via an exact match on `slurm_params.job_name` (generated unique
    per run), so concurrent runs never cross-throttle.
    """

    def __init__(self, slurm_params: "SlurmParameters", logger=logger):
        self.slurm_params = slurm_params
        self.logger = logger

    def _create_sbatch_param_str(self) -> str:
        par_list = []
        for k, v in self.slurm_params.sbatch_flag_dict().items():
            if k.startswith("--"):
                par_list.append(f"{k}={v}")
            elif k.startswith("-"):
                par_list.append(f"{k} {v}")
            else:
                # Unreachable: SlurmParameters validates extra keys start with '-'. Fail loudly
                # rather than silently drop the flag if that invariant is ever broken.
                raise ValueError(f"sbatch flag {k!r} must start with '-'.")
        return " ".join(par_list)

    def submit(self, command: str, run=True, sbatch_params_str: str | None = None) -> str:
        """Submit a single sbatch command with the configured parameters."""
        if sbatch_params_str is None:
            sbatch_params_str = self._create_sbatch_param_str()
        full_command = f"sbatch {sbatch_params_str} {command}".strip()
        if run:
            subprocess.run(full_command, shell=True, check=True)
        return full_command

    def count_jobs(self) -> int:
        """Count this run's pending+running jobs, matched exactly on the generated job name."""
        job_name_filter = self.slurm_params.job_name
        if not job_name_filter:
            raise ValueError("SlurmParameters.job_name is empty; cannot count jobs.")
        squeue_cmd = [
            "squeue",
            "-u",
            os.environ["USER"],
            "-t",
            "pending,running",
            "-o",
            "%.18i %.100j",  # job id and job name
        ]
        squeue_out = subprocess.check_output(squeue_cmd).decode()
        lines = squeue_out.strip().split("\n")[1:]
        total = 0
        for line in lines:
            parts = line.strip().split(maxsplit=1)
            if len(parts) != 2:
                continue
            job_id, job_name = parts
            if job_name != job_name_filter:
                continue
            m = re.search(r"_\[(.*)\]", job_id)
            if m:
                for part in m.group(1).split(","):
                    if "-" in part:
                        a, b = map(int, part.split("-"))
                        total += b - a + 1
                    else:
                        total += 1
            else:
                total += 1
        return total

    def watch_and_submit(
        self,
        commands: list[str],
        max_jobs_allowed: int,
        sleep_time: int = 5,
    ):
        """Submit jobs as queue slots free up, keeping at most `max_jobs_allowed` of this run's
        jobs in flight. Returns once every command has been submitted (not when they finish)."""
        os.makedirs(self.slurm_params.log_dir, exist_ok=True)
        start = 0
        total_tasks = len(commands)
        while start < total_tasks:
            current_jobs = self.count_jobs()
            available_slots = max_jobs_allowed - current_jobs
            if available_slots > 0:
                end = min(start + available_slots, total_tasks)
                for cmd in commands[start:end]:
                    if self.logger is not None:
                        self.logger.info(f"Submitting: {cmd}")
                    self.submit(cmd)
                start = end
            time.sleep(sleep_time)
