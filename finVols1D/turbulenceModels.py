"""
Various turbulence models to compute turbulent eddy viscosity
"""

from abc import ABC, abstractmethod
import numpy as np
from finVols1D.fv.fvTools import getGradCells


class turbulenceModel(ABC):

    def __init__(self, U):
        self._U = U  # velocity field


### MIXING LENGTH TURBULENCE MODEL ###
class mixingLength(turbulenceModel):

    def __init__(self, U, **kwargs):
        super(mixingLength, self).__init__(U)
        self._L = kwargs["length"]
        self._wallSide = kwargs["wall"]  # int, 0 or -1


    def _getMixingLength(self):
        """Return Prandtl mixing length, cell centers"""
        # get wall position
        self._Xwall = self._U.mesh.Xfaces[self._wallSide]
        # compute mixing length
        self._mL = self._L * np.abs(self._U.mesh.Xcells - self._Xwall)


    def mixingLength(self):
        self._getMixingLength()
        return self._mL
        

    def nut(self):
        """Return turbulent viscosity at cell centers"""
        self._getMixingLength()
        self._gradU = getGradCells(self._U)
        return self._gradU * self._mL


### SPALART ALLMARAS TURBULENCE MODEL ###


### K-OMEGA TURBULENCE MODEL ###

"""
class kOmega(turbulenceModel):

    def __init__(
            self,
            U,
            alpha=0.52,
            beta=None,
            betaStar=0.09,
            sigK=0.6,
            sigOm=0.5,
            sigD=0.125,
            **kwargs
    ):
        super(mixingLength, self).__init__(U)
        self._alpha = alpha
        self._beta = beta
        self._betaStar = betaStar,
        self._sigK = sigK
        self._sigOm = sigOm
        self._sigD = sigD
"""
