#! /usr/bin/env python3

import numpy as np
import sys
import math

import scipy
from scipy.interpolate import RectBivariateSpline

from mule.postprocessing.Sphere2DData import *


class Sphere2DDataNormsGridSpace:

    def __init__(self, filename_a = None, filename_b = None, i_output_file = None, grid_output = False, verbosity=0, output_prefix=""):

        if filename_b != None:
            self.compute_diff(filename_a, filename_b, grid_output, verbosity, output_prefix)

        if i_output_file is not None:
            self.write_file(i_output_pickle_file, verbosity, output_prefix)


    ## interpolate from fine to coarse grid
    def interpolateSolution(self,data_ref, data_cmp):

        (ny_ref, nx_ref) = data_ref.shape
        (ny_cmp, nx_cmp) = data_cmp.shape

        multiplier_j = (ny_ref)/(ny_cmp)
        multiplier_i = (nx_ref)/(nx_cmp)

        #Comparison via interpolation
        #print("Interpolation")
        # A-grid REFERENCE (file1) - sweet outputs only A grids physical space
        dx_ref=1.0/(nx_ref)
        dy_ref=1.0/(ny_ref)

        x_ref = np.arange(0, 1, dx_ref)
        x_ref = np.linspace(0, 1, nx_ref, endpoint=False)

        y_ref = np.arange(0, 1, dy_ref)
        y_ref = np.linspace(0, 1, ny_ref, endpoint=False)

        x_ref += dx_ref/2
        y_ref += dy_ref/2
        X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

        #Creat cubic interpolation of reference file
        print(y_ref.shape, x_ref.shape, data_ref.shape)
        interp_spline = RectBivariateSpline(y_ref, x_ref, data_ref)

        #A-grid cmp file (file2)
        dx_cmp=1.0/nx_cmp
        dy_cmp=1.0/ny_cmp

        x_cmp = np.arange(0, 1, dx_cmp)
        x_cmp = np.linspace(0, 1, nx_cmp, endpoint=False)

        y_cmp = np.arange(0, 1, dy_cmp)
        y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)

        x_cmp += dx_cmp/2
        y_cmp += dy_cmp/2
        X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

        #Get reduced reference resolution
        data_ref_low = interp_spline(y_cmp, x_cmp)

        return data_ref_low


    def compute_diff(
            self,
            filename_a,
            filename_b,
            grid_output,
            verbosity = 0,
            output_prefix = ""
        ):
        if grid_output:
            file_a = Sphere2DData(filename_a, setup_grid = False)
            file_b = Sphere2DData(filename_b, setup_grid = False)
        else:
            file_a = Sphere2DData(filename_a, setup_grid=True)
            file_b = Sphere2DData(filename_b, setup_grid=True)

        self.norm_l1_value = 0.0
        self.norm_l2_value = 0.0
        self.norm_linf_value = 0.0
        self.norm_rms_value = 0.0

        size_ref_j = file_a.data_grid.shape[0]
        size_ref_i = file_a.data_grid.shape[1]
        size_cmp_j = file_b.data_grid.shape[0]
        size_cmp_i = file_b.data_grid.shape[1]

        multiplier_j = (size_ref_j+1)/(size_cmp_j+1)
        multiplier_i = (size_ref_i+1)/(size_cmp_i+1)


        verbosity=6
        if verbosity > 5:
            print(output_prefix+"Dimensions of reference solution: ", size_ref_i, size_ref_j)
            print(output_prefix+"Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)

        if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
            if verbosity > 5:
                print(output_prefix+"Grids are not aligned")
                print(output_prefix+"Try to use (TODO) interpolation script")
                print(output_prefix+"Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
                print(output_prefix+"Multipliers: ", multiplier_i, multiplier_j)

                file_a.data_grid = self.interpolateSolution(file_a.data_grid, file_b.data_grid)
                multiplier_i = 1
                multiplier_j = 1
            #######raise Exception("Grids not properly aligned")

        multiplier_j = int(multiplier_j)
        multiplier_i = int(multiplier_i)

        if verbosity > 5:
            print(output_prefix+"Using multipliers (int): ", multiplier_i, multiplier_j)

        for j in range(0, size_cmp_j):
            for i in range(0, size_cmp_i):
                value = file_b.data_grid[j,i]-file_a.data_grid[j*multiplier_j,i*multiplier_i]

                # http://mathworld.wolfram.com/L1-Norm.html
                self.norm_l1_value += abs(value)
                # http://mathworld.wolfram.com/L2-Norm.html
                self.norm_l2_value += value*value
                # http://mathworld.wolfram.com/L-Infinity-Norm.html
                self.norm_linf_value = max(abs(value), self.norm_linf_value)

                # http://mathworld.wolfram.com/Root-Mean-Square.html
                self.norm_rms_value += value*value

        self.N = size_cmp_i*size_cmp_j

        # Compute sqrt() for Euklidian L2 norm
        self.norm_l2_value = math.sqrt(self.norm_l2_value)

        # RMS final sqrt(N) computation
        self.norm_rms_value  = math.sqrt(self.norm_rms_value/self.N)

        # resolution normalized L1 value
        self.res_norm_l1_value = self.norm_l1_value/float(self.N)


    def print(self, prefix=""):
        print(f"{prefix}norm l1: {self.norm_l1_value}")
        print(f"{prefix}norm l2: {self.norm_l2_value}")
        print(f"{prefix}norm linf: {self.norm_linf_value}")
        print(f"{prefix}norm rms: {self.norm_rms_value}")
        print(f"{prefix}res norm l1: {self.res_norm_l1_value}")



    def write_file(
            self,
            picklefile,
            tagprefix = None,
            verbosity = 0,
            output_prefix = ""
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        import pickle

        if tagprefix != None:
            if tagprefix[-1] != ".":
                tagprefix += '.'
        else:
            tagprefix = ""

        pickle_data = {
            tagprefix+'N' : self.N,
            tagprefix+'norm_l1' : self.norm_l1_value,
            tagprefix+'norm_l2' : self.norm_l2_value,
            tagprefix+'norm_linf' : self.norm_linf_value,
            tagprefix+'norm_rms' : self.norm_rms_value,
        }

        # Write values for resolution neutral values into the .pickle files
        pickle_data.update({
            tagprefix+'res_norm_l1' : self.res_norm_l1_value,
            tagprefix+'res_norm_l2' : self.norm_rms_value,    # This is the RMS = L2
            tagprefix+'res_norm_linf' : self.norm_linf_value,    # No normalization required
            tagprefix+'res_norm_rms' : self.norm_rms_value,    # Already normalized
        })

        pickle_data['WARNING'] = "L1, L2 and RMS don't include scaling factors for different cell spacings around the sphere2d!!!"

        print(output_prefix+"writing picklefile: "+str(picklefile))

        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(pickle_data, f)
