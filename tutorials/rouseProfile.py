"""
Solve advection diffusion transport equation
suspended sediment in a channel under equilibrium conditions
problem solution: Rouse profile
"""

import numpy as np
import matplotlib.pyplot as plt
from finVols1D import fv
from finVols1D.runTime import runTime

plt.rcParams["font.size"] = 15

# physical parameters
kappa = 0.41  # von karmann constant
Hwater = 0.1  # water depth
uf = 0.01  # friction velocity (m/s)
ws = -0.01  # settling velocity
csRef = 0.3  # reference concentration
aRef = 0.05 * Hwater  # reference height

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":10.,
     "dt":0.02}
)

# create mesh
cellFaces = np.linspace(aRef, Hwater, 200)
mesh = fv.fvMesh(cellFaces, time)


# create a field
Cs0 = np.zeros(mesh.nCells)
#Cs0 = 0.01 * (1 + np.tanh((mesh.Xcells - 0.5*Hwater)/0.002))  # initial field values
CsField = fv.fvField(
    "Cs", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedGradient", "value":0.},
    values = Cs0
)

# settling velocity field
WsField = fv.fvField(
    "Ws", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedValue", "value":0.},
    values=np.ones(mesh.nCells) * ws
)
phiWs = fv.surfaceField("Ws", mesh, WsField)

# diffusivity
nut = fv.fvField(
    "nut", mesh, time,
    bc0={"type":"fixedValue", "value":0.},
    bcN={"type":"fixedValue", "value":0.},
    values = kappa * uf * mesh.Xcells * (1 - mesh.Xcells/Hwater))
nutFaces = fv.surfaceField("nut", mesh, nut)

# prepare equations to solve
CsEqn = fv.fvEqn(mesh)

# reference concentration
csRefArr = np.zeros(mesh.nCells)
csRefArr[0] = csRef * np.abs(ws)
erosion = fv.fvField(
    "nut", mesh, time,
    values = csRefArr
)

CsEqn.addDiv(phiWs, CsField)
CsEqn.addLaplacian(nutFaces, CsField)
CsEqn.addSource(erosion)
CsField.update(CsEqn.solve())
CsEqn.reset()

Ro = np.abs(ws / (uf*kappa))  # Rouse number
print(f"Rouse number : {Ro}")
zbRef = mesh.Xcells[0]  # reference level for Rouse profile
CsRouse = csRef * (((Hwater-mesh.Xcells)/mesh.Xcells)*(zbRef/(Hwater-zbRef)))**Ro

fig, axCs = plt.subplots()

axCs.plot(CsRouse, mesh.Xcells, color="black", label="analytical")
axCs.plot(CsField.field, mesh.Xcells, label=r"$c_s$")
axCs.legend()
axCs.set_xlabel(r"$c_s$")
axCs.set_ylabel("Z")
axCs.set_xscale("log")
axCs.grid()

fig.tight_layout()
plt.show()
