import os

import law
import luigi
from generation.framework.htcondor import HTCondorWorkflow
from generation.framework.tasks import GenerationScenarioConfig, GenRivetTask
from generation.framework.utils import run_command, set_environment_variables
from generation.tasks.HerwigRun import HerwigRun
from generation.tasks.Sherpa3 import Sherpa3Run
from generation.tasks.SherpaRun import SherpaRun
from law.logger import get_logger
from luigi.util import inherits

from .Rivet4Build import Rivet4Build

logger = get_logger(__name__)


generator_class_map = {
    "herwig": HerwigRun,
    "sherpa": SherpaRun,
    "sherpa3": Sherpa3Run,
}


@inherits(GenerationScenarioConfig)
class Rivet4Run(GenRivetTask, HTCondorWorkflow, law.LocalWorkflow):
    """
    Analyze generated HEPMC files with Rivet and create YODA files
    """

    # allow outputs in nested directory structure
    output_collection_cls = law.NestedSiblingFileCollection

    # configuration variables
    files_per_job = luigi.IntParameter(
        default=10,
        description="Number of HepMC files analyzed per Rivet job. \
                Rivet is very fast analyzing HepMC files, so a sufficient high number should be given. \
                At the same time don't overdo it, since the files might be quite large and fill the scratch space.",
    )  # from RunRivet
    rivet_analyses = luigi.ListParameter(
        default=["MC_XS", "MC_WEIGHTS"],
        description="List of IDs of Rivet analyses to run.",
    )
    mc_generator = luigi.Parameter(
        default="herwig",
        description="Name of the MC generator used for event generation.",
    )

    # dummy parameter for run step
    number_of_gen_jobs = luigi.IntParameter()

    exclude_params_req = HTCondorWorkflow.exclude_params_req | {
        "files_per_job",
    }

    def workflow_requires(self):
        # Every task requires the Rivet analyses to be compiled
        return {"analyses": Rivet4Build.req(self, _exclude={"branch"})}

    def create_branch_map(self):
        # each analysis job analyzes a chunk of HepMC files
        gen_class = generator_class_map.get(str(self.mc_generator).lower())
        if gen_class is None:
            raise ValueError("Unknown MC generator: {}".format(self.mc_generator))
        return gen_class.req(
            self, number_of_jobs=self.number_of_gen_jobs
        ).get_all_branch_chunks(self.files_per_job)

    def requires(self):
        # each branch task requires existent HEPMC files to analyze
        gen_class = generator_class_map.get(str(self.mc_generator).lower())
        if gen_class is None:
            raise NotImplementedError(
                "Unknown MC generator: {}".format(self.mc_generator)
            )
        return gen_class.req(
            self,
            number_of_jobs=self.number_of_gen_jobs,
            branch=-1,
            branches=self.branch_data,
        )

    def remote_path(self, *path):
        parts = (
            self.__class__.__name__,
            str(self.mc_generator).lower(),
            self.campaign,
            self.mc_setting,
        ) + path
        return os.path.join(*parts)

    def output(self):
        dir_number = int(self.branch) / 1000
        return self.remote_target(
            "_".join(sorted(self.rivet_analyses)),
            "{DIR_NUMBER}/{INPUT_FILE_NAME}job{JOB_NUMBER}.yoda".format(
                DIR_NUMBER=int(dir_number),
                INPUT_FILE_NAME=str(self.campaign),
                JOB_NUMBER=str(self.branch),
            ),
        )

    def run(self):
        # branch data
        _map = {
            "Dijet_3_lowpt": "Dijet_3",
            "Dijet_3_highpt": "Dijet_3",
        }
        _rivet_analyses = list(self.rivet_analyses)
        _mapped_analyses = [
            _map.get(analysis, analysis) for analysis in _rivet_analyses
        ]

        # actual payload:
        print("=======================================================")
        print("Running Rivet analyses on HEPMC files ")
        print("=======================================================")

        # set environment variables
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()
        rivet_env = set_environment_variables("$ANALYSIS_PATH/setup/setup_rivet4.sh")

        # identify and get the compiled Rivet analyses
        logger.info(
            "Shared object Rivet files: {}".format(
                self.workflow_input()["analyses"]["collection"].targets.values()
            )
        )
        for so_file in self.workflow_input()["analyses"]["collection"].iter_existing():
            local_path = os.path.join(tmp_dir.abspath, so_file.basename)
            so_file.copy_to_local(local_path)

        # identify and get the HEPMC files for analyzing
        logger.info("Input events: {}".format(self.input()["collection"]))
        input_paths = self.input()["collection"].uri(base_name="root")

        yoda_out = tmp_dir.child("histos.yoda", type="f")
        cmd = ["rivet", "--pwd", f"--histo-file={yoda_out.abspath}"]
        cmd += [f"--analysis={_rivet_analysis}" for _rivet_analysis in _mapped_analyses]
        cmd += input_paths
        run_command(
            cmd,
            env=rivet_env,
            cwd=tmp_dir.abspath,
        )

        self.output().move_from_local(yoda_out)

        print("=======================================================")
