import numpy as np
import matplotlib.pyplot as plt

from finVols1D import fv
from finVols1D.runTime import runTime

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
mesh = fv.dynamicFvMesh(cellFaces, time)

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
massBed = [0.]

# prepare equations to solve
CsEqn = fv.fvEqn(mesh)

dzBed = 0.  # initial bed elevation increment

while time.loop():
    print("\n", time)
    print("bed position: ", mesh.Xfaces[0])
    mesh.meshMotion(dzBed, 0.)
    CsEqn.addDdt(CsField)
    phiWs.makeRelative()  # make flux relative to mesh motion
    CsEqn.addDiv(phiWs, CsField, scheme="linearUpwind")
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

fig, (axMass, axErr) = plt.subplots(2, figsize=(8.3, 8))

axMass.plot(CsMass, lw=2, ls="solid", color="#D55E00", label="suspension")
axMass.plot(massOut, lw=2, ls="dashed", color="#0072B2", label="massOut")
axMass.plot(massBed, lw=2, ls="dashdot", color="#009E73", label="bed")
axMass.grid()
axMass.set_ylabel("sediment mass", fontsize=15)
axMass.legend(fontsize=15)
axMass.tick_params(
    axis="both", which="major", labelsize=15)

axErr.plot(errorMassOut, lw=2, color="#0072B2", label="massOut")
axErr.plot(errorMassBed, lw=2, color="#D55E00", label="bed")
axErr.grid()
axErr.set_ylabel("relative error in %", fontsize=15)
axErr.tick_params(
    axis="both", which="major", labelsize=15)

fig.tight_layout()
plt.show()
