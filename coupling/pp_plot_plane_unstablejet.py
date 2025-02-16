#! /usr/bin/python3
# Plot unstable jet fields
# 
# --------------------------------------
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import ThreadPool

import matplotlib
from matplotlib.transforms import BboxBase

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


import numpy as np
import sys

# Figure definitions
fontsize = 18
figsize = (9, 7)

def plot(filename):
    # Load data
    print(filename)

    if 'prog_mag' in filename:
        filename_u = filename.replace('prog_mag', 'prog_u')
        data = np.loadtxt(filename_u)
        filename_v = filename.replace('prog_mag', 'prog_v')
        data2 = np.loadtxt(filename_v)
        data = np.square(data)
        data2 = np.square(data2)
        data = data + data2
        data = np.sqrt(data)

    else:
        data = np.loadtxt(filename)


    if 'prog_h_pert' in filename:
        # Get full depth
        depth = 1000
        #print(f"Summing full depth with h_0={depth}")
        depth = float(depth)
        data = data + depth

    # Get max/min
    cmin = np.amin(data)
    cmax = np.amax(data)

    # Set physical grid for axis
    x_min = 0
    x_max = 100 # In km scale
    y_min = 0
    y_max = 100
    n = data.shape[0]
    cell_size = x_max / n

    # Start plotting
    fig = plt.figure(figsize=figsize)

    # Contour levels for fields
    extent = (x_min, x_max, y_min, y_max)

    # Color plot
    # plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto')
    # plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
    # plt.imshow(data, interpolation='nearest', origin='lower')

    # Colorbar
    #if 'diag_vort' in filename:
        # plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto',
        #            cmap=plt.get_cmap('seismic'))
        #
        #
        # # Fix max and min for vorticity
        # plt.clim(cmin, cmax)
        # cmin = -5e-5
        # cmax = 5e-5
        # s = 2e-5
        # eta_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))
        # plt.clim(cmin, cmax)
        # cref = max(abs(cmin), abs(cmax))
        # plt.clim(-cref, +cref)
        # cbar = plt.colorbar(format='%.0e')
    if 'prog_h_pert' in filename:
        plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto',
                   cmap=plt.get_cmap('seismic'))
        plt.clim(900, 1100)
        cbar = plt.colorbar(format="%.0fm")
    #elif 'spec' in filename:
        # data = data + 1
        # data = np.log(data)
        # plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto',
        #            cmap=plt.get_cmap('jet'))
        # cmin = np.amin(data)
        # cmax = np.amax(data)
        # plt.clim(cmin, cmax)
        # cbar = plt.colorbar()
    elif 'prog_mag' in filename:
        plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto',
                   cmap=plt.get_cmap('seismic'))
        plt.clim(-10, +10)
        cbar = plt.colorbar()
    else:
        plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto',
                   cmap=plt.get_cmap('seismic'))
        plt.clim(-10, +10)
        cbar = plt.colorbar()

    rectangle = plt.Rectangle(((135 + 2) * cell_size, (85 + 2) * cell_size), (30 - 4) * cell_size, (30 - 4) * cell_size, linewidth=1, edgecolor='black', facecolor='none')
    plt.gca().add_patch(rectangle)
    rectangle = plt.Rectangle(((135) * cell_size, (85) * cell_size), (30) * cell_size, (30) * cell_size, linewidth=1, edgecolor='green', facecolor='none')
    plt.gca().add_patch(rectangle)


    # Contour lines (black)
    if 'diag_vort' in filename:
        pass
    #		plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
    # plt.contour(x,y,data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
    # plt.contourf(x, y, data, vmin=cmin, vmax=cmax, levels=eta_contour_levels)
    elif 'prog_h_pert' in filename:
        pass
    # plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
    # plt.contour(x,y, data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
    else:
        if cmin != cmax:
            pass
    # plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)


    # Set tittle
    title = ""
    if 'diag_vort' in filename:
        title += "Vorticity "
        cbar.set_label('1/s', rotation=270, labelpad=+20, size=fontsize)
    if 'prog_h_pert' in filename:
        title += "Depth (m) "
        # cbar.set_label('m', rotation=270, labelpad=+20, size=fontsize)
    if 'prog_u' in filename:
        title += "x-Velocity (m/s) "
        cbar.set_label('m/s', rotation=270, size=fontsize)
    if 'prog_v' in filename:
        title += "y-Velocity (m/s) "
        cbar.set_label('m/s', rotation=270, size=fontsize)
    if 'prog_mag' in filename:
        title += "Velocity Magnitude (m/s) "
        cbar.set_label('m/s', rotation=270, size=fontsize)

    cbar.ax.tick_params(labelsize=fontsize)

    # Method
    #print("Methods")
    pos1 = filename.find('_tsm_')
    pos2 = filename.find('_tso')
    method1 = filename[pos1 + 5:pos2]
    #print(method1)

    if method1 == "l_cn_na_sl_nd_settls":
        method1 = "SL-SI-SETTLS"
    elif method1 == "l_rexi_na_sl_nd_settls":
        method1 = "SL-EXP-SETTLS"
    elif method1 == "l_rexi_na_sl_nd_etdrk":
        method1 = "SL-ETD2RK"
    elif method1 == "l_rexi_n_etdrk":
        method1 = "ETD2RK"
    elif method1 == "ln_erk":
        if 'ref' in filename:
            method1 = "REF"
        else:
            method1 = "RK-FDC"

    #title += method1
    #title += " "

    title += '| t='
    pos1 = filename.find('output')
    name = filename[pos1:]
    pos2 = name.find('_t')
    pos3 = filename.find('.csv')
    time = filename[pos1 + pos2 + 2:pos3]
    time = float(time)
    time = round(time)

    title += str(time)
    title += 's'

    #print(title)
    plt.title(title, fontsize=fontsize)

    # Axis
    ax = plt.gca()
    ax.xaxis.set_label_coords(0.5, -0.075)


    #ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    #plt.grid(which='both', axis='both', alpha=1, linewidth=0.1)
    plt.xticks(fontsize=fontsize, minor=False)
    plt.yticks(fontsize=fontsize, minor=False)

    # plt.xticks(labelsx, fontsize=fontsize)
    plt.xlabel("x (km)", fontsize=fontsize)

    # plt.yticks(labelsy, fontsize=fontsize)
    plt.ylabel("y (km)", fontsize=fontsize)

    # Save file as eps
    outfilename = filename.replace('.csv', '.png')
    outfilename = outfilename.replace('/output_', '')
    outfilename = outfilename.replace('output', '')
    outfilename = outfilename.replace('../', '')

    if 'prog_u' in filename:
        outfilename = 'u/' + outfilename
    elif 'prog_v' in filename:
        outfilename = 'v/' + outfilename
    elif 'prog_mag' in filename:
        outfilename = 'mag/' + outfilename


    x1, x2, x3 = 50, 75, 95

    offset = (5,5)

    ax.scatter(x1, 50, color='black', label='Point A')
    ax.annotate('A', (x1, 50), textcoords="offset points", xytext=offset, ha='left', size="x-large")

    ax.scatter(x2, 50, color='black',label='Point B')
    ax.annotate('B', (x2, 50), textcoords="offset points", xytext=offset, ha='left', size="x-large")

    ax.scatter(x3, 50, color='black',label='Point C')
    ax.annotate('C', (x3, 50), textcoords="offset points", xytext=offset, ha='left', size="x-large")



    # print(outfilename)
    # plt.savefig(outfilename, dpi=300)
    #
    # ax.set_xlim(35, 95)
    # ax.set_ylim(20,80)

    # arrow = FancyArrowPatch((42, 30), (49, 45), color="black", mutation_scale=10)
    # ax.add_patch(arrow)
    #
    # arrow = FancyArrowPatch((25, 70), (38, 55), color="black", mutation_scale=10)
    # ax.add_patch(arrow)

    #print(outfilename)
    plt.savefig(outfilename, dpi=300, transparent=False, bbox_inches='tight', \
                pad_inches=0)

    # plt.show()
    plt.close()

# Plot in parallel
pool = multiprocessing.Pool(processes=16)
mags = [s.replace('prog_u', 'prog_mag') for s in sys.argv[1:] if 'prog_u' in s]
pool.map(plot, sys.argv[1:] + mags)

