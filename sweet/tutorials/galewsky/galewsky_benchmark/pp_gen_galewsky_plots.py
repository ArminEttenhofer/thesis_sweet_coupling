#! /usr/bin/env python3

import sys
import os

import numpy as np
import glob

from mule.JobData import *
from mule.postprocessing.Sphere2DDataOperators import Sphere2DDataOperators
from mule.postprocessing.Sphere2DDataSpectral import Sphere2DDataSpectral
import mule.plot_config as pc

data = {}

job_dirs = glob.glob("job_*")
if len(job_dirs) > 1:
    print(jobs)
    raise Exception("Multiple job directories found")
job_dir = job_dirs[0]


variables = ["phi_pert", "vrt", "div"]
for var in variables:
    file = job_dir+"/output_prog_"+var+"_t00000000000.00000000.sweet"
    sphere2d_data = Sphere2DDataSpectral(file, setup_grid=True)
    data[var] = sphere2d_data



jobdata = JobData(jobdir=job_dir).get_flattened_data()

r = float(jobdata["output.shack.sphere2d.sphere2d_radius"])
g = float(jobdata["output.shack.pde_swesphere2d.gravitation"])
h0 = float(jobdata["output.shack.pde_swesphere2d.h0"])
omega = float(jobdata["output.shack.pde_swesphere2d.sphere2d_rotating_coriolis_omega"])




ops = Sphere2DDataOperators(file_info=data[variables[0]].file_info, rsphere2d=r)

class Dummy:
    def __init__(self):
        self.data = None
        self.file_info = None

data["u"] = Dummy()
data["v"] = Dummy()
(data["u"].data_grid, data["v"].data_grid) = ops.vrtdiv2uv(data["vrt"].data_spectral, data["div"].data_spectral)


# Shortcut to latitudes
lats = data[variables[0]].file_info['lats']
lons = data[variables[0]].file_info['lons']






def plot_velocity_profile():
    """
    Plot velocity profile
    """

    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(2*s,5*s))

    ps = pc.PlotStyles()
    pstyle = ps.getNextStyle()

    pstyle["linestyle"] = "solid"
    pstyle["linewidth"] = 0.5
    pstyle["marker"] = None
    pstyle["color"] = "black"

    ax.plot(data["u"].data_grid[:,0], lats/(2*np.pi)*360, **pstyle, label="geopotential")

    pstyle["linestyle"] = "dashed"
    pstyle["linewidth"] = 0.3
    pstyle = {
            "linestyle": "dashed",
            "linewidth": 0.3,
            "color": "black",
        }

    for y in [np.pi/7, np.pi/2-np.pi/7]:
        ax.hlines(xmin=0, xmax=90, y=y/(2*np.pi)*360, **pstyle)

    ax.set_ylim([20, 70])
    ax.set_yticks([20 + 5*i for i in range(10)])

    ax.set_xlim([0, 90])
    ax.set_xticks([20*i for i in range(5)])

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Zonal Wind (m/s)")
    fig.tight_layout()

    outfile = "output_u.pdf"
    print("Saving to '"+outfile+"'")
    pc.savefig(outfile)



def plot_geopotential_profile():
    """
    Plot geopotential
    """

    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(2*s,5*s))

    ps = pc.PlotStyles()
    pstyle = ps.getNextStyle()

    pstyle["linestyle"] = "solid"
    pstyle["linewidth"] = 0.5
    pstyle["marker"] = None
    pstyle["color"] = "black"

    ax.plot(data["phi_pert"].data_grid[:,0], lats/(2*np.pi)*360, **pstyle, label="geopotential")
    ax.set_ylim([20, 70])
    ax.set_yticks([20 + 5*i for i in range(10)])

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Geopotential")
    ax.legend()
    fig.tight_layout()

    outfile = "output_geopotential.pdf"
    print("Saving to '"+outfile+"'")
    pc.savefig(outfile)




def plot_height_profile():
    """
    Plot height
    """

    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(2*s,5*s))

    ps = pc.PlotStyles()
    pstyle = ps.getNextStyle()

    pstyle["linestyle"] = "solid"
    pstyle["linewidth"] = 0.5
    pstyle["marker"] = None
    pstyle["color"] = "black"

    ax.plot(data["phi_pert"].data_grid[:,0]/g + h0, lats/(2*np.pi)*360, **pstyle, label="geopotential")
    ax.set_ylim([20, 70])
    ax.set_yticks([20 + 5*i for i in range(10)])

    ax.set_xlim([9000, None])

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Height Field (m)")
    fig.tight_layout()

    outfile = "output_h_profile.pdf"
    print("Saving to '"+outfile+"'")
    pc.savefig(outfile)



def plot_height_perturbation_profile():
    """
    Plot height perturbation
    """

    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(2*s,5*s))

    ps = pc.PlotStyles()
    pstyle = ps.getNextStyle()

    pstyle["linestyle"] = "solid"
    pstyle["linewidth"] = 0.5
    pstyle["marker"] = None
    pstyle["color"] = "black"

    data_grid = data["phi_pert"].data_grid
    data_pert = data_grid-np.repeat(np.expand_dims(data_grid[:,0], axis=1), data_grid.shape[1], axis=1)

    ilon = len(lons)//2

    ax.plot(data_pert[:,ilon]/g, lats/(2*np.pi)*360, **pstyle, label="geopotential")
    ax.set_ylim([20, 70])
    ax.set_yticks([20 + 5*i for i in range(10)])

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Height perturbation")
    ax.legend()
    fig.tight_layout()

    outfile = "output_h_perturb.pdf"
    print("Saving to '"+outfile+"'")
    pc.savefig(outfile)


def plot_height_perturbation_latlon():
    """
    Contour of perturbation
    """

    s = 0.7
    fig, ax = pc.setup(scale=1, figsize=(2.3*s,5*s))

    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Longitude (degrees)")

    data_grid = data["phi_pert"].data_grid
    data_pert = data_grid-np.repeat(np.expand_dims(data_grid[:,0], axis=1), data_grid.shape[1], axis=1)
    data_pert /= g

    # We need to map this grid to a regular grid
    from scipy.interpolate import RegularGridInterpolator

    # Source mesh
    x = lons/(2*np.pi)*360
    y = lats/(2*np.pi)*360
    grid_x, grid_y = np.meshgrid(y, x, indexing='ij')

    # Create interpolator
    interp = RegularGridInterpolator((y, x), data_pert, bounds_error=False, fill_value=None)

    # Target mesh
    N = 4*data_pert.shape[0]
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

    ax.set_xticks([180 + 50*i for i in range(-1, 2)])
    ax.set_yticks([5*i for i in range(-22,23)])

    extent = [0-pixel_angle_x, 360-pixel_angle_x, -90-pixel_angle_y, 90-pixel_angle_y]

    #if True:
    if False:
        #imhandler = ax.imshow(data_ip, cmap="viridis", extent=extent)
        imhandler = ax.imshow(np.flip(data_ip,axis=0), cmap="viridis", extent=extent)

    if True:
        levels = [10*i for i in range(1, 10)]
        ax.contour(data_ip, levels=levels, linestyles='solid', linewidths=0.2, colors='black', extent=extent)

    ax.set_ylim([20, 70])
    ax.set_xlim([-50+180, 50+180])

    fig.tight_layout()

    outfile = "output_h_perturb_latlon.pdf"
    print("Saving to '"+outfile+"'")
    pc.savefig(outfile)



def plot_vorticity_latlon():
    """
    Contour of vorticity on lat-lon field (Fig. 4)
    """

    plots = [
                {
                    "file": job_dir+"/output_prog_vrt_t00000345600.00000000.sweet",
                    "outputfile": "output_vrt_096_hours.pdf",
                },
                {
                    "file": job_dir+"/output_prog_vrt_t00000432000.00000000.sweet",
                    "outputfile": "output_vrt_120_hours.pdf",
                },
                {
                    "file": job_dir+"/output_prog_vrt_t00000518400.00000000.sweet",
                    "outputfile": "output_vrt_144_hours.pdf",
                },
            ]

    for p in plots:
        data = Sphere2DDataSpectral(p["file"], setup_grid=True)

        data_grid = data.data_grid

        s = 0.7
        fig, ax = pc.setup(scale=1, figsize=(6*s,2.5*s))

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
        ax.set_yticks([10*i for i in range(1, 9)])

        extent = [0-pixel_angle_x, 360-pixel_angle_x, -90-pixel_angle_y, 90-pixel_angle_y]

        #if True:
        if False:
            #imhandler = ax.imshow(data_ip, cmap="viridis", extent=extent)
            imhandler = ax.imshow(np.flip(data_ip,axis=0), cmap="viridis", extent=extent)

        if True:
            levels = np.array([2*1e-5*i for i in range(-20, 0)])
            ax.contour(data_ip, levels=levels, linestyles='dashed', linewidths=0.2, colors='black', extent=extent)

            levels = np.array([2*1e-5*i for i in range(1, 20)])
            ax.contour(data_ip, levels=levels, linestyles='solid', linewidths=0.2, colors='black', extent=extent)

        ax.set_ylim([10, 80])
        ax.set_xlim([0, 360])

        fig.tight_layout()

        outfile = p["outputfile"]
        print("Saving to '"+outfile+"'")
        pc.savefig(outfile)



plot_velocity_profile()
plot_geopotential_profile()
plot_height_profile()
plot_height_perturbation_profile()
plot_height_perturbation_latlon()
plot_vorticity_latlon()


