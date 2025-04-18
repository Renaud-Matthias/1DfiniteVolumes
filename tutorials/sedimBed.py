import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../src/")

from fvMesh import dynamicFvMesh
from fvFields import fvField, surfaceField
from fvEquations import fvEqn
from runTime import runTime

# physical parameters
ws = -0.01  # settling velocity
Hwater = 0.1  # water column height

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":10.,
     "dt":0.1}
)

# create mesh
ncells = 20  # number of cells
cellFaces = np.linspace(0., Hwater, ncells)
mesh = dynamicFvMesh(cellFaces, time)

# initial concentration profile
Cs0 = np.zeros(mesh.nCells) + 0.05
CsField = fvField(
    "Cs", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedGradient", "value":0},
    values = Cs0
)

# settling velocity field
wsVals = ws * np.ones(mesh.nCells)
wsField = fvField(
    "ws", mesh, time,
    bc0={"type":"fixedGradient", "value":0.},
    bcN={"type":"fixedValue", "value":0.},
    values=wsVals
)
phiWs = surfaceField("ws", mesh, wsField)


# list to store mass of sediment in domain
CsMass = [np.sum(mesh.dX * CsField.field)]
massOut = [0.]
massBed = [0.]

# prepare equations to solve
CsEqn = fvEqn(mesh)

dzBed = 0.  # initial bed elevation increment

while time.loop():
    print("\n", time)
    mesh.meshMotion(dzBed, 0.)
    print(mesh.Xfaces)
    CsEqn.addDdt(CsField)
    phiWs.makeRelative()  # make flux relative to mesh motion
    CsEqn.addDiv(phiWs, CsField)
    #CsEqn.addLaplacian(nutFaces, CsField)
    CsField.update(CsEqn.solve())
    dzBed = -phiWs[0] * CsField[0] * time._dt
    # store Cs mass in domain
    CsMass.append(np.sum(mesh.dX * CsField.field))
    # get flux of Cs through bottom boundary
    massOut.append(massOut[-1] - phiWs[0] * CsField[0] * time._dt)
    # get bed mass from boundary position
    massBed.append(mesh.Xfaces[0])
    # clean matrix before next iteration
    CsEqn.reset()

CsMass = np.array(CsMass)
massOut = np.array(massOut)
massBed = np.array(massBed)
errorMassOut = 100. * (CsMass + massOut - CsMass[0]) / CsMass[0]
errorMassBed = 100. * (CsMass + massBed - CsMass[0]) / CsMass[0]

fig, (axMass, axErr) = plt.subplots(2)

axMass.plot(CsMass, label="suspension")
axMass.plot(massOut, ls="dashed", label="massOut")
axMass.plot(massBed, ls="dashdot", label="bed")
axMass.grid()
axMass.set_ylabel("sediment mass")
axMass.legend()

axErr.plot(errorMassOut, label="massOut")
axErr.plot(errorMassBed, label="bed")
axErr.grid()
axErr.set_ylabel("relative error")

fig.tight_layout()
plt.show()
