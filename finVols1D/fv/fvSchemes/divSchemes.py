"""
class to manage various divergence schemes
"""

from finVols1D.fv.fvTools import getGradCells


def divSchemeSelector(schemeName):
    """
    Inputs:
    - schemeName: str, name of divergence scheme
    """
    if schemeName=="linear":
        return linear()
    elif schemeName=="upwind":
        return upwind()
    elif schemeName=="linearUpwind":
        return linearUpwind()
    else:
        raise ValueError(f"scheme {schemeName} not available")


# - - - LINEAR SCHEME - - - #

class linear:

    def addDiv(self, eqn, phi, field):
        """
        Inputs:
        - eqn: fvEqn, equation to modify
        - phi: surfaceField, flux through faces
        - field: fvField, variable
        """
        mesh = field.mesh
        for i in range(1, mesh.nFaces-1):
            eqn._Amat[i-1, i] += phi[i] * mesh.dX[i-1] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i-1, i-1] += phi[i] * mesh.dX[i] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i, i] -= phi[i] * mesh.dX[i-1] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i, i-1] -= phi[i] * mesh.dX[i] / (mesh.dX[i] + mesh.dX[i-1])


    def addRhoDiv(self, eqn, rho, phi, field):
        """
        Inputs:
        - eqn: fvEqn, equation to modify
        - rho: fvField, density field
        - phi: surfaceField, flux through faces
        - field: fvField, variable
        """
        mesh = field.mesh
        for i in range(1, mesh.nFaces-1):
            eqn._Amat[i-1, i] += phi[i] * mesh.dX[i-1] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i-1, i-1] += phi[i] * mesh.dX[i] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i, i] -= phi[i] * mesh.dX[i-1] / (mesh.dX[i] + mesh.dX[i-1])
            eqn._Amat[i, i-1] -= phi[i] * mesh.dX[i] / (mesh.dX[i] + mesh.dX[i-1])


# - - - UPWIND SCHEME - - - #

class upwind:

    def addDiv(self, eqn, phi, field):
        """
        Inputs:
        - eqn: fvEqn, equation to modify
        - phi: surfaceField, flux through faces
        - field: fvField, variable
        """
        for i in range(1, field.mesh.nFaces-1):
            if phi[i] >= 0:
                eqn._Amat[i-1, i-1] += phi[i]
                eqn._Amat[i, i-1] -= phi[i]
            else:
                eqn._Amat[i, i] -= phi[i]
                eqn._Amat[i-1, i] += phi[i]


# - - - LINEAR-UPWIND SCHEME - - - #
                
class linearUpwind(upwind):

    def addDiv(self, eqn, phi, field):
        """
        Inputs:
        - eqn: fvEqn, equation to modify
        - phi: surfaceField, flux through faces
        - field: fvField, variable
        """
        super(linearUpwind, self).addDiv(eqn, phi, field)
        mesh = field.mesh
        grad = getGradCells(field)
        for i in range(1, mesh.nFaces-1):
            if phi[i]>=0:
                eqn._Bvec[i] += (mesh.Xfaces[i]-mesh.Xcells[i-1]) * grad[i-1]
                eqn._Bvec[i] -= (mesh.Xfaces[i]-mesh.Xcells[i-1]) * grad[i-1]
            elif phi[i]<0:
                eqn._Bvec[i] += (mesh.Xfaces[i]-mesh.Xcells[i-1]) * grad[i]
                eqn._Bvec[i] -= (mesh.Xfaces[i]-mesh.Xcells[i-1]) * grad[i]


# - - - LIMITED SCHEME - - - #

class limited(linear, upwind):

    def __init__(self, beta=0.):
        self._beta = beta


    def addDiv(self, eqn, phi, field):
        pass
