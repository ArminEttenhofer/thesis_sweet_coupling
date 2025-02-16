from mule.JobParallelizationDimOptions import JobParallelizationDimOptions
from mule.simu.base import Benchmark


class GalewskyBenchmark(Benchmark):

    def setBenchmark(self, nModes, nProcSpace=1):
        job = self.job

        job.compile.program = "programs/PDE_SWESphere2D"


        job.compile.cart2d_spectral_space = \
            "enable" if job.compile.gui == "enable" else "disable"
        job.compile.cart2d_spectral_dealiasing = "disable"
        job.compile.sphere2d_spectral_space = "enable"
        job.compile.sphere2d_spectral_dealiasing = "enable"

        job.runtime.benchmark_name = "galewsky"
        job.runtime.space_res_spectral = nModes
        job.runtime.space_res_grid = -1

        job.compile.rexi_thread_parallel_sum = "disable"

        job.runtime.instability_checks = 0

        pSpace = JobParallelizationDimOptions("space")
        pSpace.num_cores_per_rank = nProcSpace
        pSpace.num_threads_per_rank = nProcSpace
        pSpace.num_ranks = 1

        self.pSpace = pSpace

        job.runtime.sh_setup_num_threads = pSpace.num_threads_per_rank
        job.parallelization.withNestedOpenMPNumThreads = False
        job.compile.threading = "omp"
        job.parallelization.core_oversubscription = False
        job.parallelization.core_affinity = "compact"
