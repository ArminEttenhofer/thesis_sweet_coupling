#! /usr/bin/env python3
import multiprocessing
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
from fontTools.misc.roundTools import nearestMultipleShortestRepr
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
import mule.utils as utils
from mule.postprocessing.Sphere2DDataOperators import Sphere2DDataOperators
from mule.postprocessing.Sphere2DDataSpectral import Sphere2DDataSpectral
import mule.postprocessing.Sphere2DDataOperators


def get_data_phys(filename):
    # sphere_data = SphereDataSpectral(filename, setup_physical=False)
    # ops = SphereDataOperators(file_info=sphere_data.file_info)
    # data_phys = ops.spec2phys(sphere_data.data_spectral)
    # file_info = sphere_data.file_info
    #
    # file_info['lons'] = ops.lons
    # file_info['lats'] = ops.lats

    sphere_data = Sphere2DDataSpectral(filename, setup_grid=True)
    data_phys = sphere_data.data_grid
    file_info = sphere_data.file_info

    if "h_pert" in filename:
        data_phys += 1000

    # Convert to degree
    file_info['lons'] = file_info['lons'] / (2 * np.pi) * 360
    file_info['lats'] = file_info['lats'] / (2 * np.pi) * 360

    return data_phys, file_info


def cmin_cmax_post(cmin, cmax):
    """
    Make cmin/cmax symmetric for vorticity or div fields
    """
    if "vrt_" in sys.argv[1] or "div_" in sys.argv[1]:
        print("Ensuring symmetry of cmin/cmax")
        a = max(abs(cmin), abs(cmax))
        cmin = -a
        cmax = a

    return cmin, cmax

cmin = None
cmax = None

def plot(input_filepath, cmin, cmax, data_phys, file_info):
    input_filepath_noext = os.path.splitext(input_filepath)[0]
    input_filename_noext = os.path.basename(input_filepath_noext)


    if cmin == None:
        """
        We assume that only one file will be processed for the plotting.
        Hence, we determine the min/max only from this file.
    
        We do this here since we only want to load the data once
        """
        cmin = np.min(data_phys)
        cmax = np.max(data_phys)

        cmin, cmax = cmin_cmax_post(cmin, cmax)

    if "prog_h_pert" in input_filepath:
        cmin = 900
        cmax = 1100

    # Clear plot
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(16, 8))

    # Locations of ticks
    xtickslocs = np.arange(data_phys.shape[1]) + 0.5
    ytickslocs = np.arange(data_phys.shape[0]) + 0.5

    # Labels of ticks
    xticklabels = file_info['lons']
    yticklabels = file_info['lats']

    xticklabels = np.array([round(_, 1) for _ in xticklabels])
    yticklabels = np.array([round(_, 1) for _ in yticklabels])

    assert len(xtickslocs) == len(xticklabels)
    assert len(ytickslocs) == len(yticklabels)

    if True:
        """
        Cleanup ticks so that there are only Nx ticks
        """
        Nx = 16
        N = len(xticklabels)
        step = max(1, N // Nx)
        r = np.arange(Nx, dtype=int) * step

        Ny = 8
        N = len(yticklabels)
        step = max(1, N // Ny)
        r = np.arange(Ny, dtype=int) * step

    # Make pixel centered around integer coordinates
    extent = [-0.5, data_phys.shape[1] - 0.5, data_phys.shape[0] - 0.5, -0.5]
    cmap = "viridis"
    if "prog_h_pert" in input_filepath:
        cmap = "seismic"

    extent=(0, 360, -90, 90)
    imhandler = ax.imshow(data_phys, interpolation="nearest", extent=extent, cmap=cmap, vmin=cmin, vmax=cmax)

    if 'vrt' in input_filename_noext:
        e = 2e-5
        ax.contour(data_phys, levels=np.arange(e, e * 50, e), linestyles='solid', linewidths=0.2, colors='black')
        ax.contour(data_phys, levels=np.arange(-e * 50, 0, e), linestyles='dashed', linewidths=0.2, colors='black')
    else:
        e = 2e-5
        # ax.contour(data_phys, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)

    # Fontsize
    fontsize = 18

    # Colorbar
    # ax.clim(cmin, cmax)
    cbar = fig.colorbar(imhandler, format="%dm", shrink=1)
    cbar.ax.tick_params(labelsize=fontsize)

    ax = fig.gca()

    rectangle = plt.Rectangle((360*67500/100000, -90*15000/25000), 360*15000/100000, 180*15000/50000, linewidth=2, edgecolor='green', facecolor='none')
    plt.gca().add_patch(rectangle)

    dx = 360/data_phys.shape[1]
    dy = 180/data_phys.shape[0]

    rectangle2 = plt.Rectangle((360*67500/100000 + 2*dx, -90*15000/25000 + 2 * dy), 360*15000/100000 - 4 * dx, 180*15000/50000 - 4 * dy, linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(rectangle2)

    plt.xticks(fontsize=fontsize)
    ax.set_xlabel("Longitude (Degrees)", fontsize=fontsize)
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d°'))

    plt.yticks(fontsize=fontsize)
    ax.set_ylabel("Latitude (Degrees)", fontsize=fontsize)
    ax.yaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d°'))

    pos = input_filename_noext.find('_t')
    time = input_filename_noext[pos + 2:]
    time = str(round(float(time)))

    ax.set_title("Depth (m) | t=" + time + "s", fontsize=fontsize)

    x1, x2, x3 = 50/100*360, 75/100*360, 95/100*360

    offset = (5,5)

    # ax.scatter(x1, 0, color='black', label='Point A')
    # ax.annotate('A', (x1, 0), textcoords="offset points", xytext=offset, ha='left', size="x-large")
    #
    # ax.scatter(x2, 0, color='black',label='Point B')
    # ax.annotate('B', (x2, 0), textcoords="offset points", xytext=offset, ha='left', size="x-large")
    #
    # ax.scatter(x3, 0, color='black',label='Point C')
    # ax.annotate('C', (x3, 0), textcoords="offset points", xytext=offset, ha='left', size="x-large")
    #

    output_filepath = input_filename_noext + ".png"

    fig.tight_layout()

    print("Writing to " + str(output_filepath))
    fig.savefig(output_filepath, dpi=300, bbox_inches='tight')

if len(sys.argv) > 2:
    """
    We have more than one file.
    As a first step, we determine the min/max of the data to plot it nicely
    """
    print("Preprocessing to get min/max")
    for input_file in sys.argv[1:]:
        (data_phys, file_info) = get_data_phys(input_file)

        _ = np.min(data_phys)
        if cmin == None:
            cmin = _
        else:
            cmin = min(_, cmin)

        _ = np.max(data_phys)
        if cmax == None:
            cmax = _
        else:
            cmax = max(_, cmax)


    cmin, cmax = cmin_cmax_post(cmin, cmax)
    print(f" ++ cmin: {cmin}")
    print(f" ++ cmax: {cmax}")

pool = multiprocessing.Pool(processes=16)
tasks = []
for input_filepath in sys.argv[1:]:
    input_filepath_noext = os.path.splitext(input_filepath)[0]

    # Full path without extension
    if input_filepath == input_filepath_noext:
        raise Exception("Missing file extension")

    (data_phys, file_info) = get_data_phys(input_filepath)

    tasks.append((input_filepath, cmin, cmax, data_phys, file_info))

results = pool.starmap(plot, tasks)
#
# for input_filepath in sys.argv[1:]:
#     plot(input_filepath, cmin, cmax)


