#! /usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
# protting runtime and convergence

from mule.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.variant": "normal"
})

color_palette = mcolors.TABLEAU_COLORS
colors = list(color_palette.keys())

# Group together  similar time stepping methods
groups = ['runtime.idString']

# Create plots for these variables
vars_ = ["phi_pert"]#, "vrt", "div"]

###########################################################
# User configuration ends here ############################
###########################################################

tagnames_y = []
tag_cleanup_info = []

for i in vars_:
    tagnames_y += [
        #f"sphere2d_data_diff_prog_{i}.res_norm_l1",
        #f"sphere2d_data_diff_prog_{i}.res_norm_l2",
        f"sphere2d_data_diff_prog_{i}.res_norm_linf",
    ]

    tag_cleanup_info += [
        #{"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_l1", "tag_dst": f"sphere2d_data_diff_prog_{i}.res_norm_l1"},
        #{"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_l2", "tag_dst": f"sphere2d_data_diff_prog_{i}.res_norm_l2"},
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_linf", "tag_dst": f"sphere2d_data_diff_prog_{i}.res_norm_linf"},
    ]


j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)



# Cleanup postprocessed data
JobsData_GroupsCleanupPostprocessed(job_groups, tag_cleanup_info, pickle_file_default_prefix="sphere2d_data_norms_grid_space_")


for tagname_y in tagnames_y:
    print("*"*80)
    print("Processing tagname "+tagname_y)
    print("*"*80)
    # processing convergence
    tagname_x = 'runtime.timestep_size'
    
    # Use plotting format to create (x/y) data
    d = JobsData_GroupsPlottingScattered(
            job_groups,
            tagname_x,
            tagname_y
        )

    for group_name, group_data in d.get_data_float().items():
        print("*"*80)
        print("Group: "+group_name)
        order = min(int(group_name[group_name.find("_N")+2]), int(group_name[group_name.find("_M")+2]) + 1) # min(N, M+1)

        if (group_name.find("ETDSDC") != -1):
            plt.loglog(group_data['x_values'], group_data['y_values'], '*:', mfc='w', markersize=8, color=colors[order-1], label=group_name[3:])
        else:
            plt.loglog(group_data['x_values'], group_data['y_values'], 'o:', mfc='w', markersize=6, color=colors[order-1], label=group_name[3:])
        prev_y = None
        prev_x = None
        conv = "----"

        for (x, y) in sorted(zip(group_data['x_values'], group_data['y_values']), key=lambda x: -x[0]):

            if prev_y == None or np.isnan(prev_y) or prev_y > 10e3: # what could a good value for phi-pert be? eps * 12000 or?
                conv = "----"
            else:
                rate_x = prev_x / x
                conv = f"{np.emath.logn(rate_x, prev_y/y):.2f}"

            print(f"\t{x}\t=>\t{y:.4e}\tconvergence: {conv}\t expected: " + str(order))
            prev_y = y
            prev_x = x

        if len(sys.argv) <= 1:
            print("[OK]")

    if len(sys.argv) <= 1:
        print("*"*80)
        print("Convergence tests successful")
        print("*"*80)
    
    plt.grid()
    plt.xlabel(r'$\Delta{t}$')
    plt.ylabel(r'$L_\infty$ error')
    plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0, fontsize="8")
    plt.title("Galewsky test case")
    figname = "conv_ETDSDC_order_verification.png"
    plt.savefig(figname, bbox_inches='tight', dpi=600)

