#! /usr/bin/env python3

import sys
import math
import numpy as np

from mule.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

# switch latex format on: 
#   output latex table + save plot with latex font
latex = False
# switch plotting off to pass unit tests
plotting = False

# Group together similar time stepping methods
groups = ["runtime.idString"]
color_palette = mcolors.TABLEAU_COLORS

###########################################################
# User configuration ends here ############################
###########################################################

if latex:
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino"],
        "font.variant": "normal"
    })

tagnames_y = []
tag_cleanup_info = []

tagnames_y += [
    #f"generic_data_norms.res_norm_l1",
    #f"generic_data_norms.res_norm_l2",
    f"generic_data_norms.res_norm_linf",
]

tag_cleanup_info += [
    {"ref_file_starts_with": f"output_", "tag_src": "res_norm_l1", "tag_dst": f"generic_data_norms.res_norm_l1"},
    {"ref_file_starts_with": f"output_", "tag_src": "res_norm_l2", "tag_dst": f"generic_data_norms.res_norm_l2"},
    {"ref_file_starts_with": f"output_", "tag_src": "res_norm_linf", "tag_dst": f"generic_data_norms.res_norm_linf"},
]


j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)

# Cleanup postprocessed data
JobsData_GroupsCleanupPostprocessed(job_groups, tag_cleanup_info, pickle_file_default_prefix="generic_data_norms_")


for tagname_y in tagnames_y:
    print("*"*80)
    print("Processing tagname "+tagname_y)
    print("*"*80)

    tagname_x = "runtime.timestep_size"

    d = JobsData_GroupsPlottingScattered(
            job_groups,
            tagname_x,
            tagname_y,
        )
    
    colors = list(color_palette.keys())

    for group_name, group_data in d.get_data_float().items():
    
        if latex:
            print(r"\hline", r"\multicolumn{3}{|c|}{" + group_name.replace('_', r'\_') + r"}  \\ \hline")
        else:
            print("*"*80)
            print("Group: "+group_name)

        dtSizes = group_data['x_values']
        errors = group_data['y_values']
        if plotting:
            # plot methods according to order
            if (group_name.find("ETDSDC") != -1):
                group_name = group_name.replace('_LEGENDRE','')
                order = min(int(group_name[group_name.find("N")+1]), int(group_name[group_name.find("M")+1]) + 1) # min(N, M+1)
                plt.loglog(dtSizes, errors, '*:', mfc='w', markersize=5, color=colors[order-1], label=group_name)
            elif (group_name.find("LRK") != -1):
                order = int(group_name[group_name.find("o=") + 2])
                plt.loglog(dtSizes, errors, 'd-', markersize=5, mfc='w', color=colors[order-1], label=group_name)
            else:
                order = int(group_name[group_name.find("o=") + 2])
                plt.loglog(dtSizes, errors, 'o-', markersize=5, mfc='w', color=colors[order-1], label=group_name)

        prev_value = None
        conv = "-"
        convergence_order = None
        for (x, y) in sorted(zip(group_data["x_values"], group_data["y_values"]), key=lambda x: -x[0]):

            if prev_value == None:
                conv = '-----'
            else:
                conv = f"{np.log2(prev_value/y):.2f}"

            if latex:
                print(f"{x:.5e}\t & \t{y:.16e}\t & {conv} ", r"\\ \hline")
            else:
                print(f"\t{x:.5e}\t=>\t{y:.16e}\tconvergence: {conv}\t")
            prev_value = y

    if len(sys.argv) <= 1:
        print("*"*80)
        print("Convergence tests successful")
        print("*"*80)

    if plotting:
        plt.grid()
        plt.xlabel(r'$\Delta{t}$')
        plt.ylabel(r'$L_\infty$ error')
        plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0,fontsize="8")
        plt.title("Dahlquist's test problem")
        figname = "conv_ETDRK_ETDSDC_LRK.png"
        plt.savefig(figname, bbox_inches='tight', dpi=600)
