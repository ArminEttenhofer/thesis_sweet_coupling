#! /usr/bin/env python3

import os
import sys
import numpy as np
import copy

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d)
from EFloat import *
from Functions import *
sys.path.pop()

from mule.SWEETFileDict import SWEETFileDict

import math
import struct


class REXICoefficients:

    def __init__(self):

        self.efloat = None

        self.function_name = None

        self.alphas = []
        self.betas = []
        self.gamma = None

        self.unique_id_string = None
        
        self.normalized_steady_state = False

        self.normalization_successful = None
        return

    def len(self):
        return len(self.alphas)

    def getUniqueId(self):
        return self.unique_id_string

    def toFloat(self):
        c = REXICoefficients()

        c.efloat = EFloat("float")
        c.function_name = self.function_name

        if self.gamma == None:
            c.gamma = None
        else:
            if self.efloat != None:
                if self.efloat.abs(self.efloat.im(self.gamma)) > 1e-10:
                    print("WARNING")
                    print("WARNING")
                    print("WARNING")
                    print("WARNING: Imaginary value "+str(self.efloat.im(self.gamma))+" should be close to zero")
                    print("WARNING")
                    print("WARNING")
                    print("WARNING")

                c.gamma = float(self.efloat.re(self.gamma))
            else:
                c.gamma = float(self.gamma.real)

        c.alphas = [complex(a) for a in self.alphas]
        c.betas = [complex(b) for b in self.betas]

        c.unique_id_string = self.unique_id_string

        return c


    def getNextOrderPhi(self):
        """
        Return the next order phi function
        """
        d = copy.deepcopy(self)

        # Update phi function
        assert d.function_name[0:3] == "phi"
        phiid = int(d.function_name[3:])
        new_function_name = "phi"+str(phiid+1)

        d.unique_id_string = d.unique_id_string.replace(d.function_name, new_function_name)
        d.function_name = new_function_name
        d.betas = [d.betas[i]/d.alphas[i] for i in range(len(d.alphas))]


        # There's no gamma for higher order phi functions
        d.gamma = 0
        return d


    def write_file(self, filename, binary=True):
        """
        Use dictionary
        """
        sfd = SWEETFileDict()
        
        sfd["N"] = len(self.alphas)
        sfd["function_name"] = self.function_name
        if self.gamma != None:
            sfd["gamma"] = complex(self.gamma)
        else:
            sfd["gamma"] = 0+0j
        sfd["alphas"] = np.array(self.alphas)
        sfd["betas"] = np.array(self.betas)
        
        sfd.writeToFile(filename)

        

    def symmetric_reduction(self):
        """
        Exploit a symmetry of the poles in case that they are complex conjugate symmetric
        """
        
        def print_coeffs():
            for i in range(len(self.alphas)):
                print(f"alpha: {self.alphas[i]}\t\tbeta: {self.betas[i]}")
            print("")

        print_coeffs()
        print("Number of coefficients before symmetric reduction: "+str(len(self.alphas)))

        for i in range(len(self.alphas)):
            if i >= len(self.alphas):
                break

            # Search for corresponding conjugate one
            a1 = self.alphas[i]
            for j in range(i+1, len(self.alphas)):
                if j >= len(self.alphas):
                    break

                a2 = self.alphas[j]
                if np.isclose(np.real(a1), np.real(a2)) and np.isclose(np.imag(a1), -np.imag(a2)):
                    del self.alphas[j]
                    self.betas[i] *= 2.0

        print_coeffs()
        print("Number of coefficients after symmetric reduction: "+str(len(self.alphas)))
        
        self.symmetric_reduction_applied = True
        if not '_symred' in self.unique_id_string:
                self.unique_id_string += '_symred'
        
        print("*"*80)
        print("TODO: Not yet tested!")
        print("*"*80)


    def get_beta_filtered(
        self,
        beta_threshold,
        normalized = False, # Normalize by number of REXI terms
    ):
        """
        Filter out beta coefficients which beta amplitude is below a given threshold
        """

        c = REXICoefficients()

        c.efloat = self.efloat
        c.function_name = self.function_name

        c.unique_id_string = self.unique_id_string+"_betafiltered"

        c.gamma = self.gamma
        c.alphas = []
        c.betas = []

        """
        Normalization of filter by number of poles.

        For a larger number of poles, each beta coefficient is weighted less.
        """
        fac = 1.0
        if normalized:
            fac /= len(self.alphas)

        for i in range(len(self.alphas)):
            if abs(self.betas[i]) < beta_threshold*fac:
                continue

            c.alphas.append(self.alphas[i])
            c.betas.append(self.betas[i])

        return c



    def normalize_steady_state(
        self
    ):
        """
        Rescale beta coefficients so that they the solution converges to 1 for dt -> 0
        """
        
        if self.normalized_steady_state:
            raise Exception("Normalization already triggered")
        
        if self.function_name == None:
            raise Exception("Function name required for normalization")
            
        
        x = 0
        acc = 0
        for i in range(len(self.alphas)):
            acc += self.betas[i] / (x - self.alphas[i])

        fun = Functions(self.function_name, efloat_mode='float')     
        target_value = fun.eval(0)


        eps = 1e-10

        # Check if gamma is nearby 0
        if self.gamma != None:
            if abs(self.gamma) < eps:
                self.gamma = None

        if self.gamma is not None:
            raise Exception(f"Normalization with gamma '{self.gamma}' is too dangerous. Not further continuing.")

            """
            Deactivated, since it doesn't make that much sense
            """
            if abs(self.gamma-1.0) < eps:
                self.gamma = 1

            if self.gamma == 0:
                rescale = target_value/acc
                self.betas = [i*rescale for i in self.betas]

            else:
                if abs(acc) < eps:
                    # accummulator is nearby 0 => only change gamma value to steady state
                    self.gamma = target_value

                else:
                    rescale = (target_value - self.gamma)/acc
                    self.betas = [i*rescale for i in self.betas]

        rescale = target_value/acc
        self.betas = [i*rescale for i in self.betas]

        self.unique_id_string += "_nrm"
        

    def eval(self, x: float) -> float:
        """
        Evaluate REXI.

        Parameters
        ----------
        x : float
            Point to evaluate REXI.

        Returns
        -------
        retval : float
            Evaluation of REXI.

        """
        
        if self.gamma is not None:
            retval = self.gamma
        else:
            retval = 0

        for i in range(len(self.alphas)):
            retval += self.betas[i] / (x - self.alphas[i])

        return retval


if __name__ == "__main__":

    r = REXICoefficients()
    r.efloat = EFloat()
    r.gamma = 1e-12
    r.alphas.append(1j+2.0)
    r.betas.append(3j+4.0)

    r.write_file("/tmp/test.txt")

    r2 = r.toFloat()
