"""
Class to manage finite volume fields
"""

import numpy as np
from fvTools import linInterp
import fvBoundaries as fvb


class fvField:
    
    def __init__(
            self,
            name,
            mesh,
            time=None,
            values=None,
            bc0=None,
            bcN=None
    ):
        """
        Inputs:
        - name: str, name of field
        - mesh: Mesh
        - time: runTime
        - values: ndarray, field cell center values
        - bc0: dict, entries for first boundary condition
             default, zeroGradient
        - bcN: dict, entries for second boundary condition
             default, zeroGradient
        """
        self.name = name
        self.mesh = mesh
        self.time = time
        self.field = np.zeros(self.mesh.nCells)
        self.field0 = None  # field at previous time step
        self.field00 = None  # field two time step before
        self._initialize(values=values)
        self._setBC(bc0, bcN)


    def update(self, values):
        """
        Inputs:
        - values: ndarray, field cell center values
        """
        self.field00 = self.field0
        self.field0 = self.field
        self.field = values

        
    def _initialize(self, values):
        if np.all(values)!=None:
            val = np.copy(values)
            if val.shape==(1,) or isinstance(values, float):
                self.field = val * np.ones(self.mesh.nCells)
            elif val.shape==self.field.shape:
                self.field = val
            else:
                raise ValueError(
                    "values has shape different from (1,) or (nCells,)"
                    + f"\n values.shape = {values.shape}")


    def _setBC(self, bc0Dict, bcNDict):
        """Set boundary conditions"""        
        if bc0Dict==None:
            self.bc0 = fvb.fixedGradientBC(
                {"type":"fixedValue", "value":0.}, side=0)
        else:
            self.bc0 = fvb.BCreator.factoryMethod(bc0Dict, side=0)

        if bcNDict==None:
            self.bcN = fvb.fixedGradientBC(
                {"type":"fixedValue", "value":0.}, side=-1)
        else:
            self.bcN = fvb.BCreator.factoryMethod(bcNDict, side=-1)
        # check number of cyclic boundary, must be 0 or 2
        if self.bc0.name=="cyclic" and self.bcN.name!="cyclic":
            raise ValueError(
                "boundaries with cylic condition need to be 2")
            
        

    def __add__(self, field):
        """addition fvField"""
        return self.field + field.field
    

    def __iter__(self):
        return iter(self.field)
        

    def __getitem__(self, index):
        return self.field[index]



class surfaceField:

    def __init__(self, name, mesh, fvField0):
        """
        Surface field on internal faces only
        Inputs:
        - name: str, name of surface field
        - mesh: Mesh
        - fvField0: fvField to interpolate
        """
        self.name = name
        self.mesh = mesh
        self.fvField0 = fvField0
        self.time = self.fvField0.time
        self.phi = np.zeros(self.mesh.nFaces)
        self.update(self.fvField0)


    def update(self, field):
        """
        update values of surfaceField
        Inputs:
        - field: fvField to interpolate
        """
        self.phi[1:-1] = linInterp(self.mesh, field.field)
        self.correctBC()


    def makeRelative(self):
        """
        correct flux with flux due to mesh motion
        """
        self.phi[:] -= self.mesh.phiMesh[:]


    def correctBC(self):
        self.fvField0.bc0.correctBC(self)
        self.fvField0.bcN.correctBC(self)


    def __add__(self, phi):
        """addition surfaceField"""
        return self.phi + phi.phi
    
        
    def __iter__(self):
        return iter(self.phi)
        

    def __getitem__(self, index):
        return self.phi[index]

    
    def __setitem__(self, index, value):
        self.phi[index] = value
