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


def getGradCells(field):
    """
    compute gradient of field at cell centers from faces values
    linear interpolation is used to get phi on faces
    Inputs:
    - field: fvField
    """
    mesh = field.mesh
    # get surfaceField, linear interpolation and BC
    phi = fvFields.surfaceField("phi", field.mesh, field)
    
    grad = np.zeros(mesh.nCells)  # gradient at cell centers
    grad += phi[1:]
    grad -= phi[:-1]
    grad /= mesh.dX
    return grad


def getGradFaces(field):
    """
    compute gradient of field on mesh faces
    from values of phi at neighbour cells
    Inputs:
    - field: fvField
    """
    mesh = field.mesh
    grad = np.zeros(mesh.nFaces)
    #grad[1:-1] = (field[1:] - field[:-1]) / (Zcells[1:] - Zcells[:-1])
    return grad
