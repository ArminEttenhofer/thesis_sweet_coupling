#! /usr/bin/env python3


import numpy as np
import sys
from mule.SWEETFileDict import *


class ODE_Generic_Data:
    def __init__(self, filename = None, setup_grid=False, setup_spectral=False):
        self.data = None

        self._read_file(filename)

    def _read_file(self, filename):
        
        f = SWEETFileDict()
        f.readFromFile(filename)
        if f.dict['sweetMagicCode'] != "SWEET1505":
            raise Exception("Invalid magic code")
        
        if f.dict['dataType'] == "VectorComplex":
            self.data = f.dict['data']
        else:
            raise Exception("Unknown dataType '"+f.dict['dataType']+"'")