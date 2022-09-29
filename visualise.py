import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

U = np.loadtxt("solution.dat")
assert U.shape[0]==2, "The solution file should contain exactly 2 rows"

n = U.shape[1]
m = int(np.sqrt(n)) # Reshape U into spatial coordinates
U = U.reshape(2, m, m)


fig = plt.figure(figsize=plt.figaspect(2.0))

ax = fig.add_subplot(2, 1, 1)
surf1 = ax.imshow(U[0,:,:], cmap=cm.viridis, aspect="equal", origin="lower", extent=[0,1,0,1], interpolation="antialiased")
plt.xlabel("x")
plt.ylabel("y")
plt.title("$u$")
fig.colorbar(surf1, shrink=0.6, aspect=16)

ax = fig.add_subplot(2, 1, 2)
surf2 = ax.imshow(U[1,:,:], cmap=cm.viridis, aspect="equal", origin="lower", extent=[0,1,0,1], interpolation="antialiased")
plt.xlabel("x")
plt.ylabel("y")
plt.title("$v$")
fig.colorbar(surf2, shrink=0.6, aspect=16)

plt.savefig("solution.pdf")

