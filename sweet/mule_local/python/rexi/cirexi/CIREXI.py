#! /usr/bin/env python3
#
# Author: Martin Schreiber <schreiberx@gmail.com>
# Date: 2017-08-16
#


import sys, os
import math
import cmath
import numpy as np

p = os.path.dirname(os.path.realpath(__file__))+'/../../'
sys.path.append(p)
from rexi.EFloat import *
from rexi.Functions import *
from rexi.REXICoefficients import *
sys.path.pop()




class CIREXI:
    """!
    Cauchy Contour Integration REXI method
    """

    #
    # Constructor
    # See setup(...) for documentation on parameters
    # 
    def __init__(
        self,
        efloat_mode = "float"
    ):
        self.efloat_mode = efloat_mode
        self.efloat = EFloat(efloat_mode)

        self.alphas = []
        self.betas = []
        self.unique_id_string = ""



    def setup(
        self,
        function_name,
        N: int,
        R: float = None,    # Radius
        Rx: float = None,   # X extension
        lambda_shift: complex = None,

        lambda_max_real: float = None,
        lambda_include_imag: float = None,
        pole_shift_value = None
    ):
        self.alphas = []
        self.betas = []
        self.function_name = ""

        if lambda_shift != None:
            self.setup_shifted_circle(function_name, N=N, R=R, Rx=Rx, lambda_shift=lambda_shift, pole_shift_value=pole_shift_value)

        elif lambda_max_real != None:
            self.setup_limited_circle(function_name, N=N, lambda_max_real=lambda_max_real, lambda_include_imag=lambda_include_imag, pole_shift_value=pole_shift_value)

        elif lambda_include_imag != None:
            self.setup_circle(function_name, N=N, R=lambda_include_imag, pole_shift_value=pole_shift_value)

        elif R != None:
            self.setup_shifted_circle(function_name, N=N, R=R, Rx=Rx, lambda_shift=0, pole_shift_value=pole_shift_value)

        else:
            raise Exception("Can't calculate circle radius, please provide imag/real limits.")

        coeffs = REXICoefficients()
        coeffs.alphas = self.alphas[:]
        coeffs.betas = self.betas[:]
        coeffs.gamma = 0
        coeffs.function_name = self.function_name

        coeffs.unique_id_string = self.getUniqueId()

        return coeffs



    def setup_shifted_circle(
        self,
        function_name,
        N: int,
        R: float,
        lambda_shift: complex,
        Rx: float = None,  # x extension of circle making it an ellipse
        pole_shift_value = None
    ):
        """
        Setup coefficients for Cauchy contour integral coefficients
        using a shifted circle with the shift given by lambda_shift
        which fulfills the constraints given by lambda_*

        Args:
            N (int):
                Number of points on Cauchy contour.

            R (float):
                Radius of circle for Cauchy contour

            pole_shift_value:
                value in [0;1]
                Shift of pole from current pole (=0) towards next pole (=1)
        """

        Ry = R
        if Rx is None:
            Rx = R

        self.function_name = function_name
        self.lambda_shift = lambda_shift

        self.fun = Functions(function_name, efloat_mode = self.efloat_mode)

        if pole_shift_value is None:
            pole_shift_value = 0.0

        self.pole_shift_value = pole_shift_value
 
        for j in range(N):
            x = j+pole_shift_value

            theta = x*np.pi*2.0/N

            """
            For ellipse:
            https://tutorial.math.lamar.edu/classes/calciii/vectorarclength.aspx

            fx(x) = Rx*cos(x*2*pi)
            fy(x) = Ry*sin(x*2*pi)

            fx'(x) = -2*pi*Rx*sin(x*2*pi)
            fy'(x) =  2*pi*Ry*cos(x*2*pi)

            |r'| = sqrt( fx'^2 + fy'^2 )
                 = sqrt( (Rx^2 + Ry^2) )*2*pi

            We see that this is a constant!
            """

            # Compute arclength to previous point
            # Sampling position of support point
            pos = Rx*np.cos(theta) + 1j*Ry*np.sin(theta)
            sigma_prime = -Rx*np.sin(theta) + 1j*Ry*np.cos(theta)

            # shifted position
            alpha = pos + lambda_shift
            beta = self.fun.eval(alpha)*sigma_prime
            beta *= 1j/N

            self.betas.append(beta)
            self.alphas.append(alpha)


        self.unique_id_string = "shic"
        self.unique_id_string += "_"+function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(self.efloat.re(lambda_shift))
        self.unique_id_string += "_"+str(self.efloat.im(lambda_shift))



    def setup_circle(
        self,
        function_name,
        N: int,
        R: float,
        pole_shift_value = None
    ):
        """
        Setup circle contour integral centered at origin

        Args:
            N (int):
                Number of quadrature poles

            R (float):
                Radius of circle
        """

        self.setup_shifted_circle(function_name, N, R, 0.0, pole_shift_value=pole_shift_value)

        self.unique_id_string = ""
        #self.unique_id_string = "circ"
        self.unique_id_string += function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(R)



    def setup_limited_circle(
        self,
        function_name: str,
        N: int,
        lambda_max_real: float,
        lambda_include_imag: float,
        pole_shift_value = None
    ):
        """
        Setup coefficients for Cauchy contour integral coefficients circle
        which fulfills the constraints given by lambda_*

        Args:
            N (int):
                Number of points on Cauchy contour.

            lambda_max_real (float):
                Maximum allowed real number on contour.

            lambda_include_imag (float):
                Include at least this imaginary value as part of the contour.
                Even, if the contour has to be enlarged
        """

        # If the maximal real number is larger than the max to be included imaginary number, set a lower value for the max real number
        if lambda_max_real > lambda_include_imag:
            lambda_max_real = lambda_include_imag

        x0 = lambda_max_real
        xm = lambda_include_imag
        r = (x0*x0 + xm*xm)/(2.0*x0)
        center = lambda_max_real-r

        self.setup_shifted_circle(function_name, N, r, center, pole_shift_value=pole_shift_value)

        self.unique_id_string = ""
        #self.unique_id_string = "limc"
        self.unique_id_string += function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(lambda_max_real)
        self.unique_id_string += "_"+str(lambda_include_imag)
        if self.pole_shift_value != 0:
            self.unique_id_string += "_p"+str(pole_shift_value)



    def getUniqueId(self):
        return "CIREXI_"+self.unique_id_string
