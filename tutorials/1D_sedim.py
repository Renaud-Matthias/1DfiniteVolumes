import numpy as np
import matplotlib.pyplot as plt

from finVols1D import fv
from finVols1D.runTime import runTime

plt.rcParams["font.size"] = 15


# physical parameters
ws = -0.01  # settling velocity
Hwater = 0.1  # water column height

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":15.,
     "dt":0.2}
)

# create mesh
ncells = 100  # number of cells
cellFaces = np.linspace(0., Hwater, ncells)
mesh = fv.fvMesh(cellFaces, time)

# initial concentration profile
Cs0 = np.zeros(mesh.nCells) + 0.05
CsField = fv.fvField(
    "Cs", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedGradient", "value":0},
    values = Cs0
)

# settling velocity field
wsVals = ws * np.ones(mesh.nCells)
wsField = fv.fvField(
    "ws", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedValue", "value":0.},
    values=wsVals
)
phiWs = fv.surfaceField("ws", mesh, wsField)


# list to store mass of sediment in domain
CsMass = [np.sum(mesh.dX * CsField.field)]
massOut = [0.]

# prepare equations to solve
CsEqn = fv.fvEqn(mesh)

while time.loop():
    print(time)
    CsEqn.addDdt(CsField)
    CsEqn.addDiv(phiWs, CsField)
    #CsEqn.addLaplacian(nutFaces, CsField)
    CsField.update(CsEqn.solve())
    # store Cs mass in domain
    CsMass.append(np.sum(mesh.dX * CsField.field))
    # get flux of Cs through bottom boundary
    massOut.append(massOut[-1] - wsField[0] * CsField[0] * time._dt)
    # clean matrix before next iteration
    CsEqn.reset()

CsMass = np.array(CsMass)
massOut = np.array(massOut)
massError = 100. * (CsMass + massOut - CsMass[0]) / CsMass[0]

fig, (axCs, axMass, axErr) = plt.subplots(3)

axCs.plot(Cs0, mesh.Xcells/Hwater, color="black")
axCs.plot(CsField.field, mesh.Xcells/Hwater)
axCs.set_xlabel(r"$C_s$")
axCs.set_ylabel("z/H")
axCs.grid()

axMass.plot(CsMass)
axMass.plot(massOut, ls="dashed")
axMass.grid()

axErr.plot(massError)
axErr.grid()

fig.tight_layout()
plt.show()
