#! /usr/bin/env python3

import sys
import math

from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
import mule.utils as utils

# Group together  similar time stepping methods
groups = ["runtime.timestepping_method"]


###########################################################
# User configuration ends here ############################
###########################################################


mule_plotting_usetex(False)


# Load all
jobs_data = JobsData("./job_bench_*", verbosity=0)




print("")
print("Groups:")

c = JobsDataConsolidate(jobs_data)
job_groups = c.create_groups(groups)

"""
Compute speedup and store to job data
"""
base_timing = None
for key, value in jobs_data.get_jobs_data().items():
    data = value.get_flattened_data()
    if int(data["parallelization.num_threads_per_rank"]) == 1:
        base_timing = float(data["output.simulation_benchmark_timings.main_timestepping"])

if base_timing == None:
    raise Exception("Base timing for one thread not found")

if base_timing == None:
    raise Exception("Base timing for one thread not found")

for key, value in jobs_data.get_jobs_data().items():
    data = value.get_flattened_data()
    data["output.speedup"] = base_timing/float(data["output.simulation_benchmark_timings.main_timestepping"])



"""
Prepare job data for plotting:
 * Iterate over all jobs
 * For each job:
   * Search for reference file job matching the tag "ref_file_tag"
   * Use this to lookup the error information
   * Write back the particular error information to a particular job tag
"""
JobsData_GroupsCleanupPostprocessed(job_groups, [], pickle_file_default_prefix="sphere2d_data_norms_grid_space_")


"""
The plotting data is now available at known dictionary entires (tagnames_y)

We are now ready to plot all the data.
"""

if True:

    params = []
    params += [
            {
                "tagname_x": "parallelization.num_threads_per_rank",
                "xlabel": "Number of threads",

                "tagname_y": "output.simulation_benchmark_timings.main_timestepping",
                "ylabel": "Wallclock time",

                "title": "Number of threads vs. wallclock time",

                "xscale": "linear",
                "yscale": "linear",
            },
        ]

    params += [
            {
                "tagname_x": "parallelization.num_threads_per_rank",
                "xlabel": "Number of threads",

                "tagname_y": "output.speedup",
                "ylabel": "Speedup",

                "title": f"Speedup",

                "xscale": "linear",
                "yscale": "linear",
            },
        ]

    for param in params:

        tagname_x = param["tagname_x"]
        xlabel = param["xlabel"]

        tagname_y = param["tagname_y"]
        ylabel = param["ylabel"]

        title = param["title"]

        xscale = param["xscale"]
        yscale = param["yscale"]

        print("*"*80)
        print("Processing tag "+tagname_x)
        print("*"*80)



        if True:
            """
            Plotting format
            """
            d = JobsData_GroupsPlottingScattered(
                    job_groups,
                    tagname_x,
                    tagname_y,
                )

            fileid = "output_plotting_"+tagname_x.replace(".", "-").replace("_", "-")+"_vs_"+tagname_y.replace(".", "-").replace("_", "-")


            p = Plotting_ScatteredData()


            def fun(p):
                from matplotlib import ticker
                from matplotlib.ticker import FormatStrFormatter

                plt.tick_params(axis="x", which="minor")
                # Only integer numbers
                #p.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
                #p.ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
                # Only one digit
                #p.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
                #p.ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))

                #p.ax.xaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))

                for tick in p.ax.xaxis.get_minor_ticks():
                    tick.label1.set_fontsize(8) 


                #plt.tick_params(axis="y", which="minor")
                #p.ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1e"))
                #p.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
 
                #p.ax.yaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))

                #for tick in p.ax.yaxis.get_minor_ticks():
                #    tick.label1.set_fontsize(6) 



            annotate_text_template = "{:.1f} / {:.3f}"
            p.plot(
                    data_plotting = d.get_data_float(),
                    xlabel = xlabel,
                    ylabel = ylabel,
                    title = title,
                    xscale = xscale,
                    yscale = yscale,
                    #annotate = True,
                    #annotate_each_nth_value = 3,
                    #annotate_fontsize = 6,
                    #annotate_text_template = annotate_text_template,
                    legend_fontsize = 8,
                    grid = True,
                    outfile = fileid+".pdf",
                    lambda_fun = fun,
                )

            print("Data plotting:")
            d.print()
            d.write(fileid+".csv")

        print("Info:")
        print("    NaN: Errors in simulations")
        print("    None: No data available")
