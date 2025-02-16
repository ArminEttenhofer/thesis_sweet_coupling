#! /usr/bin/env python3


import numpy as np
from mule.SWEETFileDict import SWEETFileDict


class ScalarData:
    """
    Read scalar data from .csv files:

        * with .csv ending

    """

    def __init__(self, filename = None):

        self.data = None

        if filename.endswith(".csv"):
            self._read_file_sweet(filename)
            ###self._read_file_csv(filename)

        elif filename.endswith(".sweet"):
            self._read_file_sweet(filename)

        else:
            raise Exception("Unknown file ending")

    def _read_file_csv(self, filename):

        self.data = np.loadtxt(filename)
        if self.data.size == 1:
            self.data = np.array([self.data])

    def _read_file_sweet(self, filename):

        d = SWEETFileDict(filename = filename)

        assert d["sweetMagicCode"] == "SWEET1505"
        assert d["dataType"] == "VectorComplex"

        self.data = d["data"]

