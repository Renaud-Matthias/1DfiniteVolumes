"""
Boundary conditions classes
"""

from abc import ABC, abstractmethod


class BCreator:

    def factoryMethod(BCdict, side):
        """
        Initialize correct boundary condition
        """
        if BCdict["type"] == "fixedValue":
            return fixedValueBC(BCdict, side)
        elif BCdict["type"] == "fixedGradient":
            return fixedGradientBC(BCdict, side)
        elif BCdict["type"] == "cyclic":
            return cyclicBC(BCdict, side)
        else:
            raise ValueError("boundary condition type not supported: " + BCdict["type"])

    
class fvBC(ABC):
    pass


class fixedValueBC(fvBC):

    def __init__(self, bcDict, side):
        self.name = "fixedValue"
        self._value = bcDict["value"]
        self._side = side  # value 0 if left BC, 1 if right BC


    def update(self, value):
        """change value of field at boundary"""
        self._value = value
    
        
    def correctBC(self, phi):
        """
        Inputs:
            phi: surfaceField
        """
        phi[self._side] = self._value


    def correctBCdiv(self, eqn, phi):
        """
        Inputs:
            eqn: fvEqn
            phi: surfaceField, advection velocity on faces
        """
        sign = 1.
        if self._side==0:
            sign = -1.
        eqn._Bvec[self._side] += sign * phi[self._side] * self._value


    def correctBClaplacian(self, eqn, diff):
        """
        Inputs:
            eqn: fvEqn
            diff: surfaceField, diffusivity on faces
        """
        eqn._Amat[self._side, self._side] += 2 * diff[self._side] / diff.mesh.dX[self._side]
        eqn._Bvec[self._side] += 2 * diff[self._side] * self._value / diff.mesh.dX[self._side]

        
class fixedGradientBC(fvBC):

    def __init__(self, bcDict, side):
        self.name = "fixedGradient"
        self._value = bcDict["value"]
        self._side = side


    def correctBC(self, phi):
        """
        Inputs:
            phi: surfaceField
        """
        phi[self._side] = phi.fvField0[self._side] + 0.5 * self._value * phi.mesh.Xcells[self._side]


    def correctBCdiv(self, eqn, phi):
        """
        Inputs:
            eqn: fvEqn
            phi: surfaceField, advection velocity on faces
            field: fvField
        """
        sign = 1.
        if self._side==0:
            sign = -1.
        eqn._Amat[self._side, self._side] += sign * phi[self._side]
        eqn._Bvec[self._side] += sign *0.5 * self._value * phi[self._side] * phi.mesh.dX[self._side]


    def correctBClaplacian(self, eqn, diff):
        """
        Inputs:
            eqn: fvEqn
            diff: surfaceField
            field: fvField
        """
        sign = 1.
        if self._side==0:
            sign = -1.            
        eqn._Bvec[self._side] += sign * diff[self._side] * self._value


class cyclicBC(fvBC):

    def __init__(self, bcDict, side):
        self.name = "cyclic"
        self._side = side


    def correctBC(self, phi):
        """
        Inputs:
            phi: surfaceField
        """
        # linear interpolation between first and last cells
        phi0 = phi.mesh.dX[-1] * phi.fvField0[0]
        phiN = phi.mesh.dX[0] * phi.fvField0[-1]
        phiCyclic = (phi0 + phiN) / (phi.mesh.dX[0] + phi.mesh.dX[-1])
        phi[0] = phiCyclic
        phi[-1] = phiCyclic
