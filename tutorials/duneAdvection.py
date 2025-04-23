import numpy as np
import matplotlib.pyplot as plt
from finVols1D import fv
from finVols1D.runTime import runTime

plt.rcParams["font.size"] = 15

# physical parameters
Qwater = 1.  # water discharge
Hwater = 1.  # water level
xdune = 0.5  # dune crest position
hdune = 0.1  # initial dune height
duneWidth = 0.2

# parameter for bedload flux
alpha = 0.05
beta = 1.5

def bedload(zb):
    qb = alpha * Qwater**beta / (Hwater-zb)**beta
    return qb

def celerityDune(zb):
    cdune = alpha * beta * Qwater**beta / (H-zb)**(beta+1)
    return cdune

Lx = 1.  # domain length

# create time control
time = runTime(
    {"startTime":0.,
     "endTime":1.,
     "dt":0.01}
)

# create mesh
cellFaces = np.linspace(0., 1., 100)
mesh = fv.fvMesh(cellFaces, time)

# create a field
zb0 = hdune * np.exp(-(mesh.Xcells - xdune)**2/duneWidth**2)
zbField = fv.fvField(
    "zb", mesh, time,
    bc0={"type":"zeroGradient"},
    bcN={"type":"zeroGradient"},
    values = zb0
)

# settling velocity field
qbField = fv.fvField(
    "qb", mesh, time,
    bc0={"type":"zeroGradient"},
    bcN={"type":"zeroGradient"},
    values=bedload(zbField[:])
)
phib = fv.surfaceField("qb", mesh, qbField)

# prepare equations to solve
zbEqn = fv.fvEqn(mesh)

while time.loop():
    print(f"time: {time.time}")
    zbEqn.addDdt(zbField)
    zbEqn.addDiv(phib, zbField)
    #CsEqn.addLaplacian(nutFaces, CsField)
    zbField.update(zbEqn.solve())
    zbEqn.reset()

fig, ax = plt.subplots()
ax.plot(mesh.Xcells, zb0, color="black")
ax.plot(mesh.Xcells, zbField.field)
ax.set_xlabel(r"$z_b$")
ax.set_ylabel("x")
ax.grid()
plt.show()
