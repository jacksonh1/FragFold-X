import subprocess
import time
import os
from typing import Callable
from functools import partial
import json
import re
from loguru import logger
import fragfoldx.config as config

# from loguru import Logger


class SlurmJobSubmitter:
    """
    A class to handle Slurm job submission with sbatch parameters.
    """

    def __init__(
        self,
        sbatch_param_file: str | None = None,
        extra_sbatch_params: dict | None = None,
        logger=logger,
    ):
        """
        Initializes the SlurmJobSubmitter with optional sbatch parameters from a file or a string.
        extra_sbatch_params overwrite the sbatch parameters from the file if a specific sbatch parameter is present in both.
        i.e. if sbatch_param_file contains "--time":"1:00:00" and extra_sbatch_params contains "--time":"2:00:00", the final sbatch parameters will be "--time":"2:00:00".
        """
        self.sbatch_param_file = sbatch_param_file
        self.extra_sbatch_params = extra_sbatch_params
        self.sbatch_params = self._combine_sbatch_params()
        self.logger = logger

    def _read_sbatch_param_file_json(self) -> dict:
        """
        Reads the sbatch parameters from a file and returns them as a dictionary with keys as sbatch parameters and values as their settings.
        """
        if self.sbatch_param_file is None:
            raise ValueError("No sbatch parameter file provided.")
        if not os.path.exists(self.sbatch_param_file):
            raise FileNotFoundError(f"File {self.sbatch_param_file} does not exist.")
        with open(self.sbatch_param_file, "r") as f:
            sbatch_params = json.load(f)
        return sbatch_params

    def _combine_sbatch_params(self) -> dict:
        """
        Creates a string of sbatch parameters from the file and/or string.
        """
        params = {}
        if self.sbatch_param_file is not None:
            params = self._read_sbatch_param_file_json()
        if self.extra_sbatch_params is not None:
            params.update(self.extra_sbatch_params)
        return params

    def _create_sbatch_param_str(self) -> str:
        par_list = []
        for k, v in self.sbatch_params.items():
            if k.startswith("--"):
                par_list.append(f"{k}={v}")
            elif k[0] == "-":
                par_list.append(f"{k} {v}")
        return " ".join(par_list)

    def submit(self, command: str, run=True, sbatch_params_str: str | None = None) -> str:
        """
        Submits a single sbatch command with the given parameters.
        """
        if sbatch_params_str is None:
            sbatch_params_str = self._create_sbatch_param_str()
        full_command = f"sbatch {sbatch_params_str} {command}".strip()
        if run:
            subprocess.run(full_command, shell=True, check=True)
            return full_command
        return full_command

    def count_jobs(self, count_only_job_name: bool = True) -> int:
        if count_only_job_name:
            if "--job-name" in self.sbatch_params:
                job_name_filter = self.sbatch_params["--job-name"]
            else:
                raise ValueError(
                    "No job name provided in sbatch parameters so cannot count jobs with `count_only_job_name=True`."
                )
        else:
            job_name_filter = None
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
            if job_name_filter is not None and job_name_filter not in job_name:
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
        count_only_job_name: bool = True,
    ):
        """
        Watches the job queue and submits new jobs as slots become available.
        """
        start = 0
        total_tasks = len(commands)
        while start < total_tasks:
            current_jobs = self.count_jobs(count_only_job_name=count_only_job_name)
            available_slots = max_jobs_allowed - current_jobs
            if available_slots > 0:
                end = start + available_slots
                if end >= total_tasks:
                    end = total_tasks
                if start <= end:
                    cmds = commands[start:end]
                    for cmd in cmds:
                        if self.logger is not None:
                            self.logger.info(f"Submitting: {cmd}")
                        self.submit(cmd)
                    start = end
            time.sleep(sleep_time)


colabfold_sbatch_submitter = SlurmJobSubmitter(
    sbatch_param_file=config.COLABFOLD_SBATCH_PARAM_FILE, logger=logger
)
