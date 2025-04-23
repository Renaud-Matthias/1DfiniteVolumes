import numpy as np
import matplotlib.pyplot as plt

from finVols1D import fv
from finVols1D.runTime import runTime

plt.rcParams["font.size"] = 15

# physical parameters
Lx = 1.  # domain length

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":0.1,
     "dt":0.001}
)

# create mesh
cellFaces = np.linspace(0., Lx, 100)
mesh = fv.fvMesh(cellFaces, time)

# create a field, U0: initial condition
U0 = -1 + 0.2 * np.cos(2*np.pi * mesh.Xcells/Lx)
Ufield = fv.fvField(
    "u", mesh, time,
    bc0={"type":"cyclic"},
    bcN={"type":"cyclic"},
    values = U0
)

phiU = fv.surfaceField("u", mesh, Ufield)

# prepare equations to solve
UEqn = fv.fvEqn(mesh)

while time.loop():
    print(time)
    UEqn.addDdt(Ufield)
    UEqn.addDiv(phiU, Ufield)
    Ufield.update(UEqn.solve())
    phiU.update(Ufield)
    UEqn.reset()


fig, ax = plt.subplots()
ax.plot(mesh.Xcells, U0, color="black")
ax.plot(mesh.Xcells, Ufield.field)
ax.set_xlabel(r"$u$")
ax.set_ylabel("X")
ax.grid()
plt.show()
