"""
Solve velocity profile in straight channel
problem solution: Rouse profile
"""

import numpy as np
import matplotlib.pyplot as plt

from finVols1D import fv
from finVols1D.runTime import runTime
from finVols1D import turbulenceModels

plt.rcParams["font.size"] = 15


# physical parameters, for fluid
kappa = 0.41  # von karmann constant
nu = 1e-6  # kinematic water viscosity
Hwater = 0.1  # water depth
uf = 0.02  # friction velocity (m/s)
g = 9.81  # gravity acceleration
Uobj = 1.  # mean objective velocity m/s

# sediment properties
ws = -0.03  # settling velocity, m/s
csRef = 0.3  # reference concentration
aRef = 0.05 * Hwater  # reference height

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":5.,
     "dt":0.01}
)

# create mesh
cellFaces = np.linspace(0., Hwater, 100)
mesh = fv.fvMesh(cellFaces, time)

# create velocity field
U0 = np.zeros(mesh.nCells)
Ufield = fv.fvField(
    "U", mesh, time,
    bc0={"type":"fixedValue", "value":0.},
    bcN={"type":"fixedGradient", "value":0.},
    values = U0
)
# flux of U through mesh faces
phiU = fv.surfaceField("U", mesh, Ufield)

# create suspended concentration field
Cs0 = np.zeros(mesh.nCells)
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
# settling flux through mesh faces
phiWs = fv.surfaceField("Ws", mesh, WsField)

# instantiate turbulence model
turbulence = turbulenceModels.mixingLength(Ufield, length=Hwater, wall=0)

# turbulent eddy viscosity
nut = fv.fvField(
    "nut", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedGradient", "value":0.},
    values = turbulence.nut()
    )
nutFaces = fv.surfaceField("nut", mesh, nut)

# mean velocity force
Usource = fv.fvField(
    "F_Umean", mesh, time,
    values = np.ones(mesh.nCells) * (Uobj-Ufield.field)
)

# turbulent diffusivity for sediment
diffSed = fv.fvField(
    "diffSed", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedGradient", "value":0.},
    values = turbulence.nut()
    )
diffSedFaces = fv.surfaceField("diffSed", mesh, diffSed)

# reference concentration and erosion condition
csRefArr = np.zeros(mesh.nCells)
erosion = fv.fvField(
    "nut", mesh, time,
    values = csRefArr
)

# instantiate equations to solve
UEqn = fv.fvEqn(mesh)
CsEqn = fv.fvEqn(mesh)

while time.loop():
    print("\n", time)
    # solve equation for velocity U
    print("solve equation for U")
    Usource.update(np.ones(mesh.nCells) * (Uobj-Ufield.field))
    UEqn.addDdt(Ufield)
    UEqn.addLaplacian(nutFaces, Ufield)
    UEqn.addSource(Usource)
    Ufield.update(UEqn.solve())
    UEqn.reset()
    nut.update(turbulence.nut())
    nutFaces.update(nut)
    # solve suspended load
    print("solve equation for Cs")
    diffSed.update(turbulence.nut())
    diffSedFaces.update(diffSed)
    gradU = fv.fvTools.getGradCells(Ufield)
    uf = np.sqrt(nu * gradU[0])
    print(f"friction velocity, ustar = {uf*100} cm/s")
    csRefArr[0] = csRef * np.abs(ws)  # reference concentration
    erosion.update(csRefArr)
    CsEqn.addDdt(CsField)
    CsEqn.addDiv(phiWs, CsField, scheme="linearUpwind")
    CsEqn.addLaplacian(diffSedFaces, CsField)
    CsEqn.addSource(erosion)
    CsField.update(CsEqn.solve())
    CsEqn.reset()


fig, (axU, axGradU, axML, axCs) = plt.subplots(ncols=4)

axU.plot(Ufield.field, mesh.Xcells/Hwater)
axU.set_xlabel(r"$u$")
axU.set_ylabel("z/H")
axU.grid()

axGradU.plot(fv.fvTools.getGradCells(Ufield), mesh.Xcells/Hwater)
axGradU.set_xlabel(r"$\partial_zu$")
axGradU.grid()

axML.plot(nut.field, mesh.Xcells/Hwater)
axML.set_xlabel(r"$\nu_t$")
axML.grid()

axCs.plot(CsField.field, mesh.Xcells/Hwater)
axCs.set_xlabel(r"$c_s$")
axCs.set_xscale("log")
axCs.grid()

fig.tight_layout()
plt.show()
