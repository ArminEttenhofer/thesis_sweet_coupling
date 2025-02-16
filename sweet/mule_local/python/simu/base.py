from abc import ABC, abstractmethod

from mule.JobMule import JobGeneration
from mule.JobParallelizationDimOptions import JobParallelizationDimOptions

from mule.sdc import getParamsSDC


class Benchmark(ABC):

    def __init__(self, verbose=True):
        self.job = JobGeneration()
        self.verbose = True
        self.pSpace = None
        self.pTime = None


    def setSimulation(self, tEnd,
                      compileMode="release", gui="disable", verbosity=0,
                      timeLimit=None):
        job = self.job

        job.compile.mode = "release"
        job.compile.gui = "disable"
        job.runtime.verbosity = verbosity

        job.runtime.max_simulation_time = tEnd

        timeLimit = 60*60*1 if timeLimit is None else timeLimit
        job.parallelization.max_wallclock_seconds = timeLimit

        job.unique_id_filter = [
            "compile",
            "runtime.max_simulation_time",
            "runtime.benchmark",
            "compile_cart2d",
        ]
        job.compilecommand_in_jobscript = False


    def setOutput(self, dt, fileMode='bin'):
        job = self.job

        job.runtime.output_timestep_size = dt
        job.runtime.output_file_mode = fileMode


    @abstractmethod
    def setBenchmark(self):
        print("LALALALALALALALALALALA")


    def setTimeStepping(self, timeStepping, dt, **kwargs):
        job = self.job
        self.pTime = None

        job.runtime.timestepping_method = timeStepping
        job.runtime.timestep_size = dt

        if "SDC" in timeStepping:

            paramsSDC = kwargs.get("paramsSDC", getParamsSDC())
            job.runtime.sdc_params = paramsSDC
            job.runtime.sdc_file = "params_SDC.sweet"

            nProcTime = kwargs.get("nProcTime", 1)

            if nProcTime > 1:
                pTime = JobParallelizationDimOptions("time")
                pTime.num_cores_per_rank = nProcTime
                pTime.num_threads_per_rank = nProcTime
                pTime.num_ranks = 1
                job.runtime.sdc_parallel = 1
                pTime.print()

                self.pTime = pTime
            else:
                job.runtime.sdc_parallel = 0


    def writeJobDir(self, dirName=None):
        job = self.job

        # Setup parallelization
        parConfs = [p for p in [self.pSpace, self.pTime] if p is not None]
        if len(parConfs) > 0:
            job.setup_parallelization(parConfs)
            if self.verbose:
                print("Parallelization output:")
                job.parallelization.print()

        # Generate job directory
        jobDir = "job_bench_"+(job.getUniqueID() if dirName is None else dirName)
        job.gen_jobscript_directory(jobDir)

        # SDC specific
        if "SDC" in job.runtime.timestepping_method:
            paramsSDC = job.runtime.sdc_params
            paramsSDC.writeToFile(job.job_dirpath+"/"+job.runtime.sdc_file)

        # Generate compilation command
        job.write_compilecommands(
            content_accum=job.get_compilecommands_accum())
