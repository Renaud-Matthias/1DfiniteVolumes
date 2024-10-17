"""
A simple stationnary heat diffusion problem
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../src/")

from fvMesh import fvMesh
from fvFields import fvField, surfaceField
from fvEquations import fvEqn
from runTime import runTime


# create time control
time = runTime(
    {"startTime":0.,
     "endTime":10.,
     "dt":0.1}
)

# create mesh
cellFaces = np.linspace(0, 1, 100)
mesh = fvMesh(cellFaces, time)


# create a field
T0 = np.zeros(mesh.nCells)  # initial field values
Tfield = fvField(
    "T", mesh, time,
    bc0={"type":"fixedValue", "value":0.},
    bcN={"type":"fixedGradient", "value":1.},
    values = T0
)

# diffusivity
diffCells = fvField("D", mesh, time, 0.01)
diffFaces = surfaceField("D", mesh, diffCells)

# prepare equations to solve
TEqn = fvEqn(mesh)
TEqn.addLaplacian(diffFaces, Tfield)
Tfield.update(TEqn.solve())


fig, ax = plt.subplots()
ax.plot(mesh.Xcells, Tfield.field)
ax.set_xlabel("X")
ax.set_ylabel("T")
ax.grid()
plt.show()
