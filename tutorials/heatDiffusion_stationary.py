"""
A simple stationnary heat diffusion problem
"""

import numpy as np
import matplotlib.pyplot as plt
from finVols1D import fv
from finVols1D.runTime import runTime

plt.rcParams["font.size"] = 15


# create time control
time = runTime(
    {"startTime":0.,
     "endTime":10.,
     "dt":0.1}
)

# create mesh
cellFaces = np.linspace(0, 1, 100)
mesh = fv.fvMesh(cellFaces, time)


# create a field
T0 = np.zeros(mesh.nCells)  # initial field values
Tfield = fv.fvField(
    "T", mesh, time,
    bc0={"type":"fixedValue", "value":0.},
    bcN={"type":"fixedGradient", "value":1.},
    values = T0
)

# diffusivity
diffCells = fv.fvField("D", mesh, time, 0.01)
diffFaces = fv.surfaceField("D", mesh, diffCells)

# prepare equations to solve
TEqn = fv.fvEqn(mesh)
TEqn.addLaplacian(diffFaces, Tfield)
Tfield.update(TEqn.solve())


fig, ax = plt.subplots()
ax.plot(mesh.Xcells, Tfield.field)
ax.set_xlabel("X")
ax.set_ylabel("T")
ax.grid()
plt.show()
