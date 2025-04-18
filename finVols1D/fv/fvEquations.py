"""
class to manage partial 1D differential equation
"""

import numpy as np
from finVols1D.fv.fvSchemes.divSchemes import divSchemeSelector
from finVols1D.fv.fvTools import linInterp, courantNo


class fvEqn:

    def __init__(self, mesh):
        """
        Inputs:
            mesh: Mesh
        """
        self._mesh = mesh
        # initialize matrix and source vector
        self.reset()
        # initialize limiter array, 0 is upwind scheme, 1 is linear, and 2 is downwind
        self._lim = np.zeros(self._mesh.nFaces)


    def addDdt(self, field, scheme=None):
        """
        Add a temporal derivative term in equation
        """
        if field.time==None:
            raise ValueError(
                "trying to add time scheme on stationary problem\n"
                + f"to run transient case, field {field.name} "
                + "must be associated to a time object -> "
                + "field(name, mesh, time, ..."
            )
        time = field.time
        for i in range(self._mesh.nCells):
            self._Amat[i, i] += self._mesh.dX[i] / (time.time - time.time_1)
            self._Bvec[i] += field.field[i] * self._mesh.dX0[i] / (time.time - time.time_1)


    def addRhoDdt(self, rho, field, scheme=None):
        """
        Add a temporal derivative term in equation
        """
        if field.time==None:
            raise ValueError(
                "trying to add time scheme on stationary problem\n"
                + f"to run transient case, field {field.name} "
                + "must be associated to a time object -> "
                + "field(name, mesh, time, ..."
            )
        time = field.time
        for i in range(self._mesh.nCells):
            self._Amat[i, i] += self._mesh.dX[i] * rho.field[i] / (time.time - time.time_1)
            self._Bvec[i] += rho.field0[i] * field.field[i] * self._mesh.dX0[i] / (time.time - time.time_1)

    
    def addDiv(self, phi, field, scheme="upwind"):
        """
        Add a divergence term to matrix system
        Inputs:
            phi: surfaceField
            field: fvField
            scheme: default is upwind
        """
        meanCo, maxCo, minCo = courantNo(phi, phi.time._dt)
        print(f"- Courant number, div({phi.name},{field.name}): mean={round(meanCo, 5)}"
              + f", max={round(maxCo, 5)}, min={round(minCo, 5)}")
        divScheme = divSchemeSelector(scheme)
        divScheme.addDiv(self, phi, field)
        
        if field.bc0.name == "cyclic":
            if phi[0] >= 0:
                self._Amat[0, -1] -= phi[0]
                self._Amat[-1, -1] += phi[0]
            else:
                self._Amat[0, 0] -= phi[0]
                self._Amat[-1, 0] += phi[0]
        else:
            field.bc0.correctBCdiv(self, phi)
            field.bcN.correctBCdiv(self, phi)


    def addRhoDiv(self, rho, phi, field, scheme="upwind"):
        """
        Add a divergence term to matrix system
        Inputs:
            rho: fvField, density field
            phi: surfaceField
            field: fvField
            scheme: default is upwind
        """
        meanCo, maxCo, minCo = courantNo(phi, phi.time._dt)
        print(f"- Courant number, div({phi.name},{field.name}): mean={round(meanCo, 5)}"
              + f", max={round(maxCo, 5)}, min={round(minCo, 5)}")
        divScheme = divSchemeSelector(scheme)
        divScheme.addRhoDiv(self, rho, phi, field)
        
        if field.bc0.name=="cyclic":
            if phi[0] >= 0:
                self._Amat[0, -1] -= phi[0]
                self._Amat[-1, -1] += phi[0]
            else:
                self._Amat[0, 0] -= phi[0]
                self._Amat[-1, 0] += phi[0]
        else:
            field.bc0.correctBCdiv(self, phi)
            field.bcN.correctBCdiv(self, phi)


    def addLaplacian(self, diff, field, scheme=None):
        """
        Add a laplacian term to matrix system
        Inputs:
            diff: surfaceField, diffusivity
            field: fvField
            scheme: ???
        """
        # internal field
        for i in range(1, self._mesh.nCells-1):
            self._Amat[i, i-1] -= diff[i] / (
                self._mesh.Xcells[i]-self._mesh.Xcells[i-1])
            self._Amat[i, i] += diff[i] / (
                self._mesh.Xcells[i]-self._mesh.Xcells[i-1])
            self._Amat[i, i] += diff[i+1] / (
                self._mesh.Xcells[i+1]-self._mesh.Xcells[i])
            self._Amat[i, i+1] -= diff[i+1] / (
                self._mesh.Xcells[i+1]-self._mesh.Xcells[i])
        self._Amat[0, 0] += diff[1] / (
            self._mesh.Xcells[1]-self._mesh.Xcells[0])
        self._Amat[0, 1] -= diff[1] / (
            self._mesh.Xcells[1]-self._mesh.Xcells[0])
        self._Amat[-1, -2] -= diff[-2] / (
            self._mesh.Xcells[-1]-self._mesh.Xcells[-2])
        self._Amat[-1, -1] += diff[-2] / (
            self._mesh.Xcells[-1]-self._mesh.Xcells[-2])
        # boundary conditions
        field.bc0.correctBClaplacian(self, diff)
        field.bcN.correctBClaplacian(self, diff)


    def addSource(self, field):
        """
        add explicit source term 
        Inputs:
            field: fvField, source term
        """
        self._Bvec[:] += field.field[:]

            
    def solve(self):
        """solve matrix system and return field values"""
        return np.linalg.solve(self._Amat, self._Bvec)


    def reset(self):
        """Reset the matrix system to zero"""
        self._Amat = np.zeros((self._mesh.nCells, self._mesh.nCells))
        self._Bvec = np.zeros(self._mesh.nCells)
