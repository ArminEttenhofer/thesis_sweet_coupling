#! /usr/bin/env python3


import numpy as np


class Sphere2DData:
    """
    Read sphere2d data from arbitrary sources:

        * with .sweet ending

        * with .csv ending

        * in grid/spectral space

    Warning: Not everything is supported, yet
    """

    def __init__(self, filename = None, setup_grid=False, setup_spectral=False):

        self.original_data_space = None

        self.data_grid = None
        self.data_spectral = None

        self.setup_grid = setup_grid
        self.setup_spectral = setup_spectral

        if filename.endswith(".sweet"):
            self._read_file_sweet(filename)
            setup_grid = True

        elif filename.endswith(".csv"):
            self._read_file_csv(filename)
            setup_grid = False
            ###raise Exception("TODO")

        else:
            raise Exception("Unknown file ending")

        if setup_grid:
            # Convert to grid space
            if self.data_grid == None:
                self._spectral_to_grid()

        if setup_spectral:
            # Convert to spectral space
            if self.data_spectral == None:
                self._grid_to_spectral()


    def _spectral_to_grid(self):
        from mule.postprocessing.Sphere2DDataOperators import Sphere2DDataOperators

        ops = Sphere2DDataOperators(file_info=self.file_info)
        self.data_grid = ops.spec2phys(self.data_spectral)

        self.file_info['lats'] = ops.lats
        self.file_info['lons'] = ops.lons

    def _grid_to_spectral(self):
        from mule.postprocessing.Sphere2DDataOperators import Sphere2DDataOperators

        raise Exception("TODO")

        ops = Sphere2DDataOperators(file_info=self.file_info)
        self.data_spectral = ops.phys2spec(self.data_grid)

        self.file_info['lats'] = ops.lats
        self.file_info['lons'] = ops.lons


    def _read_file_csv(self, filename):
        ## skip first row and column (lat-lon values)
        self.data_grid = np.loadtxt(filename)[1:, 1:]


    def _read_file_sweet(self, filename, setup_grid=True):
        f = open(filename, 'rb')
        content = f.read()
        f.close()

        self.file_info = {}

        acc = ""
        line_nr = 0
        fin_detected = False

        self.original_data_space = None

        i = 0
        for i in range(0, len(content)):
            c = content[i]

            if c == ord('\n'):
                if line_nr == 0:
                    if acc == "SWEET":
                        pass
                    else:
                        raise Exception("Header not detected, stopping")

                else:
                    s = "DATA_TYPE"
                    if acc[0:len(s)] == s:
                        self.file_info['data_type'] = acc[len(s)+1:]

                    s = "MODES_N_MAX"
                    if acc[0:len(s)] == s:
                        self.file_info['modes_n_max'] = int(acc[len(s)+1:])

                    s = "MODES_M_MAX"
                    if acc[0:len(s)] == s:
                        self.file_info['modes_m_max'] = int(acc[len(s)+1:])

                    s = "NUM_ELEMENTS"
                    if acc[0:len(s)] == s:
                        self.file_info['num_elements'] = int(acc[len(s)+1:])

                    s = "GRID_TYPE"
                    if acc[0:len(s)] == s:
                        self.file_info['grid_type'] = acc[len(s)+1:]

                    s = "FIN"
                    if acc[0:len(s)] == s:
                        #print("FIN detected")
                        fin_detected = True
                        i += 1
                        break

                acc = ""
                line_nr += 1

            else:
                acc += chr(c)

        if not fin_detected:
            raise Exception("FIN not detected in file")

        if self.file_info['data_type'] != "SH_DATA":
            print("ERROR: Only binary data supported")
            raise Exception(
                f"Wrong data type {self.file_info['data_type']} in binary file")

        if 0:
            print("*"*80)
            for key, value in self.file_info.items():
                print(str(key)+": "+str(value))
            print("*"*80)

        data = content[i:]

        self.data_spectral = np.frombuffer(data, dtype=np.cdouble)
