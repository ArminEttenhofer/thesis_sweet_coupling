#! /usr/bin/env python3
import glob
import sys

from mule.postprocessing.JobData import JobData, InfoError


class JobsData(InfoError):

    def __init__(
        self,
        job_dirs_pattern = 'job_bench_*',
        job_dirs = [],
        verbosity = 10,
        withReferenceJobs = True,
        useUniqueId = True
    ):
        """
        Load the data of all jobs in the current folder

        Parameters:
        -----------
            job_dirs: string / list
                if string:
                    Pattern with wildcards to autodetect job directories, e.g. 'job_bench_*'
        """
        InfoError.__init__(self, 'JobsData')

        self.withReferenceJobs = withReferenceJobs
        self.verbosity = verbosity
        self.useUniqueId = useUniqueId

        # If job_dirs is just a string, it's a search pattern
        if len(job_dirs) > 0:
            pass

        elif isinstance(job_dirs_pattern, str):
            job_dirs = glob.glob(job_dirs_pattern)

        else:
            raise Exception("Unknown type for job_dirs")

        # Load raw job data
        self.__load_job_raw_data(job_dirs)

        # Create flattened data
        self.__create_flattened_data()

    @property
    def jobNames(self):
        return list(self.__jobs_data.keys())

    def __getitem__(self, item) -> JobData:
        try:
            if isinstance(item, int):
                return self.__jobs_data[self.jobNames[item]]
            else:
                return self.__jobs_data[item]
        except (KeyError, IndexError):
            raise ValueError(f'{item} not in JobsData, only got {self.jobNames}')


    def __load_job_raw_data(
            self,
            job_dirs = []
    ):
        """
        Parse all output.out files and extract all kind of job output information

        Return a dictionary with content from the job directories
        {
            [name of job directory] :
            {
                #
                # Dictionary with data from job generation
                # (read from [jobdir]/jobgeneration.pickle)
                #
                'jobgeneration':
                {
                    'compile': [...],
                    'runtime': [...],
                    'parallelization': [...],
                    'platforms_platform': [...],
                    'platform_resources': [...],
                },
                'output':
                {
                    'SimulationBenchmarkTimings.main': [value],
                    'SimulationBenchmarkTimings.main_simulationLoop': [value],
                    [...]
                }
            }
        """

        self.__jobs_data = {}
        for job_dir in job_dirs:
            if self.verbosity > 5:
                self.info("")
                self.info("Processing '"+job_dir+"'")


            job = JobData(job_dir, verbosity=self.verbosity)
            flattened_data = job.get_flattened_data()

            if not self.withReferenceJobs:
                if 'jobgeneration.reference_job' in flattened_data:
                    if flattened_data['jobgeneration.reference_job']:
                        continue

            if 'jobgeneration.job_unique_id' not in flattened_data:
                # Be backward compatible
                self.__jobs_data[job_dir] = job
            else:
                if self.useUniqueId:
                    job_unique_id = job.get_flattened_data()['jobgeneration.job_unique_id']
                    self.__jobs_data[job_unique_id] = job
                else:
                    self.__jobs_data[job_dir] = job



    def get_jobs_data(self):
        """
        Return the raw job information data
        Warning: This storage format is likely to change!
        """
        return self.__jobs_data



    def __create_flattened_data(self):
        """
        Return a dictionary with a flattened hierarchy

        This joins all hierarchical structures with a '.'

        E.g. the value
            self.jobs_raw_data['jobgeneration'].parallelization.num_cores_per_socket
        is stored in
            retval['jobgeneration.parallelization.num_cores_per_socket']
        """

        self.__flattened_jobs_data = {}

        # For all jobs
        for job_key, job_data in self.__jobs_data.items():
            self.__flattened_jobs_data[job_key] = job_data.get_flattened_data()



    def get_flattened_data(self):
        return self.__flattened_jobs_data


    def __str__(self):
        retstr = ""
        for key, value in self.__flattened_jobs_data.items():
            retstr += str(key)+"\n"
        return retstr


if __name__ == "__main__":

    verbosity = 10
    verbosity = 0

    from mule.postprocessing.JobsDataConsolidate import JobsDataConsolidate

    if True:
        j = JobsData(verbosity=verbosity)
#        c = JobsDataConsolidate(verbosity=verbosity)
        d = j.get_flattened_data()

        for jobdir, job_data in d.items():
            print("*"*80)
            print("Data for '"+jobdir+"'")
            print("*"*80)
            for key, value in job_data.items():
                print(key+" => "+str(value))

    sys.exit(1)

    if False:
        j = JobsData(verbosity=verbosity)

        jobs_raw_data = j.get_jobs_raw_data()
        for key, values in jobs_raw_data.items():
            print(key)
            print(values)
            print("")

    #if True:
    if False:
        j = JobsData(verbosity=verbosity)
        d = j.get_flattened_data()

        for jobdir, job_data in d.items():
            print("")
            print("Data for '"+jobdir+"'")
            for key, value in job_data.items():
                print(key+" => "+str(value))

    #if True:
    if False:
        j = JobsData(verbosity=verbosity)
        d = j.get_flattened_data()

        # Just print the first item
        #if False:
        if True:
            for jobdir, job_data in d.items():
                print("")
                print("Data for '"+jobdir+"'")
                for key, value in job_data.items():
                    print(key+" => "+str(value))
                break

        g = j.create_groups(['runtime.timestepping_method'])

        #for group_id, group_jobs in g.items():
        for group_id in sorted(groups):
            group_jobs = groups[group_id]

            print(group_id)
            for jobdir, jobdata in group_jobs.items():
                print(" + "+jobdir+"\t"+jobdata['output.benchmark_timings.main_simulationloop'])


    #if True:
    if False:
        j = JobsData(verbosity=verbosity)

        data_table = j.create_data_table(
                ['runtime.timestepping_method'],
                'parallelization.num_threads_per_rank',
                'output.benchmark_timings.main_simulationloop'
            )

        print("Data table:")
        j.print_data_table(data_table)
