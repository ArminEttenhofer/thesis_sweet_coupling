import numpy as np

from mule.postprocessing.JobData import JobData
from mule.postprocessing.Sphere2DData import Sphere2DData
from mule.postprocessing.Sphere2DDataOperators import Sphere2DDataOperators


class Solution(object):
    TYPE = "SPHERE2D"

    def __init__(self, job:JobData, time=None):

        self.job = job
        self.rSphere = int(job.output['shack.Sphere2D.sphere2d_radius'])

        self._div:Sphere2DData = None
        self._vrt:Sphere2DData = None
        self._phi_pert:Sphere2DData = None

        self.ops:Sphere2DDataOperators = None
        self.time:float = None

        if time is not None:
            self.load(time)

    @property
    def times(self):
        self.job.getOutputFiles(majorKey='TIME')
        return self.job.outputTimes

    def load(self, time:float):
        """
        Load solution at a given simulation time

        Parameters
        ----------
        time : float, optional
            Simulation time. If not given, take the one already loaded.
        """
        files = self.job.getOutputFiles(majorKey='TIME')
        if time not in self.times:
            raise ValueError(f'cannot load solution for t={time}, '
                             f'got only {self.times}')

        self._div = Sphere2DData(files[time]['div'])
        self._vrt = Sphere2DData(files[time]['vrt'])
        self._phi_pert = Sphere2DData(files[time]['vrt'])

        self.ops = Sphere2DDataOperators(
            rsphere2d=self.rSphere, file_info=self._div.file_info)

        self.time = time

    def get(self, var:str, spectral:bool=True, time:float=None) -> np.ndarray:
        """
        Get the solution data for a given variable, eventually at given time,
        in spectral or grid space.

        Parameters
        ----------
        var : str
            Name of the variable ('div', 'vrt' or 'phi_pert').
        spectral : bool, optional
            Get data in spectral space. The default is True.
        time : float, optional
            Simulation time. If not given, take the one already loaded.

        Returns
        -------
        data : np.ndarray
            The data in array form (1D for spectral, 2D for grid space)
        """
        if time is not None:
            self.load(time)
        data:Sphere2DData = None
        if var == "div":
            data = self._div
        elif var == "vrt":
            data = self._vrt
        elif var == "phi_pert":
            data = self._phi_pert
        else:
            raise ValueError(f'cannot get {var} variable')
        if spectral:
            if data.data_spectral is None:
                data._grid_to_spectral()
            return data.data_spectral
        else:
            if data.data_grid is None:
                data._spectral_to_grid()
            return data.data_grid

    def computeKineticEnergySpectrum(self, time:float=None) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute the kinetic energy spectrum at a some simulation time,
        following Eq. (5) in [1]_.

        Parameters
        ----------
        time : float, optional
            Simulation time. If not given, take the one already loaded.

        Returns
        -------
        spectrum : np.ndarray
            The 1D energy spectrum values
        wavelengths : np.ndarray
            The associated wavelengths (in km)

        .. [1] Koshyk, J. N., & Hamilton, K. (2001).
               "The horizontal kinetic energy spectrum and spectral budget simulated by a
               high-resolution troposphere-stratosphere-mesosphere GCM."
               Journal of the Atmospheric Sciences, 58(4), 329–348.
               https://doi.org/10.1175/1520-0469(2001)058<0329:THKESA>2.0.CO;2.
        """
        if time is None and self.time is None:
            raise ValueError('no time loaded, need some')

        if time is not None:
            self.load(time)

        sh = self.ops
        vrt = self.get('vrt', spectral=True)
        div = self.get('div', spectral=True)

        def getArrayIndexByModes(n, m):
            return (m*(2*sh.ntrunc-m+1)>>1)+n

        spectrum = np.zeros(sh.lats.shape[0])
        for m in range(sh.ntrunc):
            idx = getArrayIndexByModes(m, m)
            for n in range(m, sh.ntrunc):
                v = vrt[idx]
                d = div[idx]

                if n != 0:
                    spectrum[n] += np.real(1/4*self.rSphere**2/(n*(n+1))*(
                        v*np.conj(v) + d*np.conj(d)))

                idx += 1

        modes = np.arange(len(spectrum))
        modes[0] = -1.0 # avoid div/0
        wavelengths = np.pi * 2.0 * self.rSphere / modes

        spectrum = spectrum[1:-1]
        wavelengths = wavelengths[1:-1]
        wavelengths /= 1e3

        return spectrum, wavelengths
