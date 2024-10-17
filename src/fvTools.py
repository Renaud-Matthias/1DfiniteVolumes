"""
Function to handle matrix for 1D finite volumes
"""

import numpy as np


def courantNo(phi, dt):
    """
    Inputs:
        phi: surfaceField
        dt: float, time step
    """
    Co = np.abs(phi[1:-1]) * dt / np.abs(phi.mesh.Xcells[1:] - phi.mesh.Xcells[:-1])
    meanCo, maxCo, minCo = np.mean(Co), np.max(Co), np.min(Co)
    return meanCo, maxCo, minCo
    

def linInterp(mesh, field):
    """
    return value of field phi on internal faces through a linear interpolation
    Inputs:
        mesh: fvMesh
        field: fvField, field to interpolate
    """
    phiLI = mesh.dX[1:] * field[:-1] + mesh.dX[:-1] * field[1:]
    phiLI /= (mesh.dX[1:] + mesh.dX[:-1])
    return phiLI


def getGradCells(mesh, phi, phiTop=0, phiBot=0):
    """
    compute gradient of field phi at cell centers from faces values
    linear interpolation is used to get phi on faces
    Inputs:
        mesh: fvMesh
        phi: surfaceField
    """
    # mesh data
    ncells = len(dzCells)
    # compute gradient
    phiF = linInterp(mesh, phi)  # phi value on faces, linear interpolation
    grad = np.zeros(mesh.nCells)  # gradient at cell centers
    grad[:-1] += phiF
    grad[1:] -= phiF
    # boundary condition
    grad[0] -= phiBot
    grad[-1] += phiTop
    grad /= dzCells
    return grad


def getGradFaces(mesh, phi):
    """compute gradient of field phi at internal faces
    from values of phi at neighbour cells"""
    # mesh data
    Zcells = mesh["zcells"]
    ncells = len(Zcells)
    # compute gradient
    grad = (phi[1:] - phi[:-1]) / (Zcells[1:] - Zcells[:-1])
    return grad
