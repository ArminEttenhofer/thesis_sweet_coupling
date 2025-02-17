#! /usr/bin/env python3


import copy
import sys
import numpy
import mule.rexi.brexi.rk_co as rk_co

import mule.rexi.EFloat as ef
from mule.rexi.REXICoefficients import *
from mule.rexi.Functions import *




class BREXI:
    """
    B-REXI method

    Implements the diagonalization of the Butcher table method to get a REXI
    formulation
    """

    def __init__(
        self,
        efloat_mode = None,
    ):
        self.efloat = ef.EFloat(efloat_mode)



    def setup(self,
              N: int = 1,
              quadrature_method: str = "gauss_legendre",
              ) -> REXICoefficients:
        """
        Set up coefficients for BREXI method.

        Parameters
        ----------
        N : int, optional
            Number of quadrature points. The default is 1.
        quadrature_method : str, optional
            Quadrature points. The default is "gauss_legendre".

        Returns
        -------
        REXICoefficients
            DESCRIPTION.
        """
        self.N = N
        self.quadrature_method = quadrature_method
        
        if quadrature_method is not None:
            if quadrature_method.startswith('gauss_'):
                self.quadrature_method_short = quadrature_method[6:9].upper()
            else:
                self.quadrature_method_short = quadrature_method[0:3].upper()
        else:
            self.quadrature_method_short = "NP"

        co = rk_co.rk_co(N, quadrature_method)

        """
        Step 1) Diagonalization
        """

        """
        Compute Eigendecomposition with Eigenvalues in D
        """
        D, E = numpy.linalg.eig(co.A)

        """
        Compute diagonal-only W so that
            W * E^-1 * ONES = ONES
        with ONES a vector with only ones
        """
        #
        # E^-1 * ONES = eta
        # <=> E*eta = ONE
        # <=> eta = solve(E, ONE)
        #
        eta = numpy.linalg.solve(E, numpy.ones(N))

        #
        # W * eta = ONES
        # diag-only W
        # diag(W) = eta^-1
        # diag(W_inv) = eta
        #
        W_inv = numpy.diag(eta)
        b_tilde = co.b.T.dot(E.dot(W_inv))

        # Create modified/unified REXI form
        """
        Step 2) Get REXI formulation
        """

        gamma = 1.0 - numpy.sum(b_tilde/D)
        alphas = 1/D
        betas = -b_tilde/(D*D)

        self.coeffs = REXICoefficients()
        self.coeffs.function_name = "phi0"
        self.coeffs.efloat = self.efloat
        self.coeffs.alphas = alphas
        self.coeffs.betas = betas
        self.coeffs.gamma = gamma

        self.coeffs.unique_id_string = self.getUniqueId()

        return self.coeffs

    def getCoeffs(self):
        return self.coeffs


    def getUniqueId(self):
        return "BREXI_phi0_"+self.quadrature_method_short+"_"+str(self.N)



if __name__ == "__main__":
    numpy.set_printoptions(precision=20)

    for method in ['gauss_legendre', 'gauss_chebyshev_u', 'gauss_chebyshev_t', 'gauss_hermite', 'gauss_lobatto', 'gauss_jacobi']:

        print("Method: "+method)

        brexi = BREXI()

        N=16

        coeffs = brexi.setup(N=N, quadrature_method=method)
        filepre = 'brexi_'+method


        # Convert to floating point
        coeffs = coeffs.toFloat()


        print("Alphas:")
        for i in coeffs.alphas:
            print(" + "+str(i))
        print("")

        print("Betas:")
        for i in coeffs.betas:
            print(" + "+str(i))
        print("")

        print("Gamma:")
        print(coeffs.gamma)
        print("")


        if True:
        #if False:
            import matplotlib
            matplotlib.use('Agg')
            import numpy as np
            import matplotlib.pyplot as plt

            x = np.real(coeffs.alphas)
            y = np.imag(coeffs.alphas)

            plt.clf()
            plt.scatter(x, y)
            plt.savefig('output_'+filepre+'_alphas.pdf')

            plt.clf()
            plt.scatter(x, y, s=np.log(np.absolute(coeffs.betas))+1.0)
            plt.savefig('output_'+filepre+'_alphas_scaled.pdf')

            x = np.real(coeffs.betas)
            y = np.imag(coeffs.betas)

            plt.clf()
            plt.scatter(x, y)
            plt.savefig('output_'+filepre+'_betas.pdf')


            #
            # Error plot
            #
            plt.clf()
            function = Functions(
                function_name = "phi0",
                efloat_mode = "float"
            )

            test_range = [-N, N]
            num_test_samples = 4096

            max_error = 0
            xpts = np.linspace(test_range[0], test_range[1], num_test_samples)
            yvals = np.zeros(num_test_samples)

            for i in range(num_test_samples):
                x = xpts[i]
                lam = 1j*x

                y = function.eval(lam)
                yn = coeffs.eval(lam)

                err = np.abs(y-yn)

                yvals[i] = err

            plt.plot(xpts, yvals, 'r-')

            plt.xlim(test_range[0], test_range[1])
            plt.ylim(1e-12, 10)
            plt.yscale("log")

            plt.savefig('output_'+filepre+'_error_oscillatory.pdf')
