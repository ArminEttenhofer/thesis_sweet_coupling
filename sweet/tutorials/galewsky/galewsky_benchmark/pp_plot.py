#! /usr/bin/env python3

import sys
import os


import matplotlib.pyplot as plt
import numpy as np
from mule.postprocessing.Sphere2DDataSpectral import Sphere2DDataSpectral
import mule.plot_config as pc


def get_data_grid(filename):
    from mule.postprocessing.Sphere2DDataSpectral import Sphere2DDataSpectral

    sphere2d_data = Sphere2DDataSpectral(filename, setup_grid=True)
    data_grid = sphere2d_data.data_grid
    file_info = sphere2d_data.file_info

    # Convert to degree
    file_info['lons'] = file_info['lons']/(2*np.pi)*360
    file_info['lats'] = file_info['lats']/(2*np.pi)*360


    return (data_grid, file_info)



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

if len(sys.argv) > 2:
    """
    We have more than one file.
    As a first step, we determine the min/max of the data to plot it nicely
    """
    print("Preprocedding to get min/max")
    for input_file in sys.argv[1:]:
        print(" + "+input_file)

        (data_grid, file_info) = get_data_grid(input_file)

        _ = np.min(data_grid)
        if cmin == None:
            cmin = _
        else:
            cmin = min(_, cmin)

        _ = np.max(data_grid)
        if cmax == None:
            cmax = _
        else:
            cmax = max(_, cmax)

        print(f" ++ cmin: {cmin}")
        print(f" ++ cmax: {cmax}")

    cmin, cmax = cmin_cmax_post(cmin, cmax)




for input_filepath in sys.argv[1:]:

    input_filepath_noext = os.path.splitext(input_filepath)[0]

    data = Sphere2DDataSpectral(input_filepath, setup_grid=True)
    data_grid = data.data_grid

    # Shortcut to latitudes
    lats = data.file_info['lats']
    lons = data.file_info['lons']

    if cmin == None:
        """
        We assume that only one file will be processed for the plotting.
        Hence, we determine the min/max only from this file.

        We do this here since we only want to load the data once
        """
        cmin = np.min(data_grid)
        cmax = np.max(data_grid)

        cmin, cmax = cmin_cmax_post(cmin, cmax)

    # Setup figure
    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(6*s,3*s))

    ps = pc.PlotStyles()

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Longitude (degrees)")

    # We need to map this grid to a regular grid
    from scipy.interpolate import RegularGridInterpolator

    # Source mesh
    x = lons/(2*np.pi)*360
    y = lats/(2*np.pi)*360
    grid_x, grid_y = np.meshgrid(y, x, indexing='ij')

    # Create interpolator
    interp = RegularGridInterpolator((y, x), data_grid, bounds_error=False, fill_value=None)

    # Target mesh
    N = 4*data_grid.shape[0]
    xi = np.linspace(0, 360, N*2)
    yi = np.linspace(-90, 90, N)
    xi_mesh, yi_mesh = np.meshgrid(yi, xi, indexing='ij')

    # Create a flattened array of target grid points
    points = np.column_stack((xi_mesh.ravel(), yi_mesh.ravel()))

    # Interpolate the values at the target grid points
    data_ip = interp(points)

    # Reshape the interpolated values to match the target grid dimensions
    data_ip = data_ip.reshape(xi_mesh.shape)

    # Make pixel centered around integer coordinates
    pixel_angle_x = 360/data_ip.shape[1]
    pixel_angle_y = 180/data_ip.shape[0]

    ax.set_xticks([30 + 50*i for i in range(8)])
    ax.set_yticks([-90+30*i for i in range(0, 7)])

    extent = [0-pixel_angle_x, 360-pixel_angle_x, -90-pixel_angle_y, 90-pixel_angle_y]

    if True:
        imhandler = ax.imshow(np.flip(data_ip,axis=0), cmap="viridis", vmin=cmin, vmax=cmax, extent=extent)

        # Colorbar
        cbar = fig.colorbar(imhandler, ax=ax)
        cbar.ax.tick_params() 


    if True:
        levels = np.array([2*1e-5*i for i in range(-20, 0)])
        ax.contour(data_ip, levels=levels, linestyles='dashed', linewidths=0.2, colors='black', extent=extent)

        levels = np.array([2*1e-5*i for i in range(1, 20)])
        ax.contour(data_ip, levels=levels, linestyles='solid', linewidths=0.2, colors='black', extent=extent)


    # Fontsize
    ax.set_title(input_filepath_noext.split("/")[-1])

    fig.tight_layout()

    output_filepath = input_filepath_noext+".png"
    print("Saving to '"+output_filepath+"'")
    pc.savefig(output_filepath, dpi=600)

