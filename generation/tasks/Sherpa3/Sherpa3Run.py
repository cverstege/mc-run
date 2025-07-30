import os

import law
import luigi
from generation.framework.htcondor import HTCondorWorkflow
from generation.framework.tasks import GenerationScenarioConfig, GenRivetTask
from generation.framework.utils import run_command, set_environment_variables
from law.decorator import safe_output
from law.logger import get_logger
from luigi.util import inherits

from .Sherpa3Integrate import Sherpa3Integrate

logger = get_logger(__name__)


@inherits(GenerationScenarioConfig)
class Sherpa3Run(GenRivetTask, HTCondorWorkflow, law.LocalWorkflow):
    """
    Use the prepared grids in Herwig-cache to generate HEP particle collision \
    events
    """

    # allow outputs in nested directory structure
    output_collection_cls = law.NestedSiblingFileCollection

    # configuration variables
    start_seed = luigi.IntParameter(
        default=42,
        description="Start seed for the MC generator. "
        "The first job will use this seed, while the nth job is using this seed + n.",
    )
    number_of_jobs = luigi.IntParameter(
        default=1,
        description="Number of individual generation jobs. Each will generate statistically independent events.",
    )
    events_per_job = luigi.IntParameter(
        default=10000, description="Number of events generated in each job."
    )

    exclude_params_req = HTCondorWorkflow.exclude_params_req | {
        "setupfile",
        # "number_of_jobs",
        "events_per_job",
        "start_seed",
    }

    def workflow_requires(self):
        # Each job requires the sherpa setup to be present
        return {
            "Sherpa3Integrate": Sherpa3Integrate.req(self, _exclude={"branch"}),
        }

    def create_branch_map(self):
        return {i: self.start_seed + i for i in range(self.number_of_jobs)}

    def remote_path(self, *path):
        parts = (self.__class__.__name__, self.campaign, self.mc_setting) + path
        return os.path.join(*parts)

    def output(self):
        dir_number = int(self.branch) / 1000
        return self.remote_target(
            "{DIR_NUMBER}/job{JOB_NUMBER}.root".format(
                DIR_NUMBER=int(dir_number),
                JOB_NUMBER=str(self.branch),
            )
        )

    @safe_output
    def run(self):
        _num_events = str(self.events_per_job)
        seed = int(self.branch_data)

        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        # get the prepared Sherpack and runfiles and unpack them
        with self.workflow_input()["Sherpa3Integrate"].localize("r") as _file:
            os.system("tar -xzf {} -C {}".format(_file.path, tmp_dir.abspath))

        # Keep total events constant, distribute across CPUs
        cpu_count = int(os.environ.get("PYTHON_CPU_COUNT", os.cpu_count() / 4))
        events_per_cpu = max(1, int(_num_events) // cpu_count)

        # set environment variables
        sherpa_env = set_environment_variables(
            os.path.expandvars("$ANALYSIS_PATH/setup/setup_sherpa3.sh")
        )
        _sherpa_exec = ["mpirun", "-n", str(cpu_count), "Sherpa"]
        _sherpa_args = [
            "-R {SEED}".format(SEED=seed),
            "-e {NEVENTS}".format(NEVENTS=events_per_cpu),
        ]

        if self.mc_setting == "withNP":
            _gen_opts = []
        elif self.mc_setting == "NPoff":
            _gen_opts = ["--framgentation None", "--mi-handler None"]
        elif self.mc_setting == "Hadoff":
            _gen_opts = ["--framgentation None"]
        elif self.mc_setting == "MPIoff":
            _gen_opts = ["--mi-handler None"]
        else:
            raise ValueError("Unknown mc_setting: {}".format(self.mc_setting))

        run_command(
            _sherpa_exec + _sherpa_args + _gen_opts, env=sherpa_env, cwd=tmp_dir.abspath
        )
        root_files = tmp_dir.glob("*.root*")

        with self.output().localize("w") as output:
            hadd_command = ["hadd", "-f", output.abspath] + root_files
            run_command(
                hadd_command,
                env=sherpa_env,
                cwd=tmp_dir.abspath,
            )
