"""
finite volume mesh object
motion of mesh prescribing boundary faces motion
"""

import numpy as np
from fvEquations import fvEqn
from fvFields import fvField, surfaceField
from runTime import runTime

class fvMesh:
    
    def __init__(self, Xfaces, time):
        """
        Inputs:
        - Xfaces: ndarray, mesh faces coordinates
        - time: runTime
        """
        self.time = time
        self.Xfaces = Xfaces
        self.nFaces = len(Xfaces)
        self.nCells = self.nFaces - 1
        self.Xcells = self._getCellCenters()
        self._getCellWidths()
        self.dX0 = np.copy(self.dX)


    def _getCellCenters(self):
        """compute cell center from cell edges"""
        return 0.5 * (self.Xfaces[1:] + self.Xfaces[:-1])


    def _getCellWidths(self):
        """compute cell widths from cell edges"""
        self.dX = np.abs(self.Xfaces[1:] - self.Xfaces[:-1])


class dynamicFvMesh(fvMesh):

    def __init__(self, Xfaces, time):
        """
        Inputs:
        - Xfaces: ndarray, mesh faces coordinates
        - time: runTime
        """
        super(dynamicFvMesh, self).__init__(Xfaces, time)
        # field of cell center displacements
        self.dXc = fvField(
            name="dXc", mesh=self, time=self.time,
            values=np.zeros(self.nCells),
            bc0 = {"type":"fixedValue", "value":0.},
            bcN = {"type":"fixedValue", "value":0.}
        )
        # mesh diffusivity
        self.fvDiff = fvField(
            "meshDiff", self, time=self.time,
            values=np.ones(self.nCells))
        self.diff = surfaceField(
            "meshDiff", self, self.fvDiff)
        # matrix system to solve displacement
        self.eqn = fvEqn(self)
        # mesh velocity
        self.Umesh = fvField(
            "Umesh", self, time=self.time,
            values=np.zeros(self.nCells),
            bc0={"type":"fixedValue", "value":0.},
            bcN={"type":"fixedValue", "value":0.}
        )
        # flux due to mesh motion
        self.phiMesh = surfaceField(
            "phiMesh", self, self.Umesh)


    def meshMotion(self, dX0, dXN):
        """
        Inputs:
        - dX0: float, displacement of boundary 0
        - dXN: float, displacement of boundary N
        """
        self.dXc.bc0.update(dX0)
        self.dXc.bcN.update(dXN)
        self.eqn.addLaplacian(self.diff, self.dXc)
        self.dXc.update(self.eqn.solve())
        print("mesh motion, mean cell displacement: "
              + f"{np.mean(self.dXc.field)} m")
        self._updateMesh(dX0, dXN)
        self.eqn.reset()
    

    def _updateMesh(self, dX0, dXN):
        """
        update mesh when mesh is moving
        Inputs:
        - dX0: float, displacement of boundary 0
        - dXN: float, displacement of boundary N"""
        self.dX0 = np.copy(self.dX)
        # change cell positions
        self.Xcells += self.dXc[:]
        self.Xfaces[1:-1] = 0.5 * (self.Xcells[1:] + self.Xcells[:-1])
        self.Xfaces[0] += dX0
        self.Xfaces[-1] += dXN
        # compute new cell widths
        self._getCellWidths()
        self.Umesh.update(self.dXc[:] / (
            self.time.time - self.time.time_1))
        self.Umesh.bc0.update(dX0 / (
            self.time.time - self.time.time_1))
        self.Umesh.bcN.update(dXN / (
            self.time.time - self.time.time_1))
        self.phiMesh.update(self.Umesh)
        
