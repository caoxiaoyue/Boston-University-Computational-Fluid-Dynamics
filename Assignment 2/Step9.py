from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from initial_conditions5_8 import *
from differences import *

# Initial problem parameters
xDomain = (0.0, 2.0)  # x domain
yDomain = (0.0, 1.0)  # y domain
nx = 60  # number of x-grid points
ny = 60  # number of y-grid points
nit = 1000  # number of iterations
X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
dx = max(X) / (nx - 1)  # delta x
dy = max(Y) / (ny - 1)  # delta y

# Create an empty array for the pressure across the domain
p = np.zeros((nx, ny), dtype=np.float64)
pn = np.zeros((nx, ny), dtype=np.float64)

# Set the boundary conditions
p[0, :] = 0.0
p[nx - 1, :] = Y
p[:, 0] = p[:, 1]  # Neumann conditions
p[:, ny - 1] = p[:, ny - 2]  # ... same as above

# Explicit scheme with Central Difference in space (5-point difference)
for iit in range(nit):
    pn = p
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            p[i, j] = ((dy ** 2) * (pn[i + 1, j] + pn[i - 1, j]) +
                       (dx ** 2) * (pn[i, j + 1] + pn[i, j - 1])) / (2 * (dx ** 2 + dy ** 2))
    p[0, :] = 0.0
    p[nx - 1, :] = Y
#    p[1:nx - 2, 0] = p[1:nx - 2, 1]
#    p[1:nx - 2, ny - 1] = p[1:nx - 2, ny - 2]
    p[:, 0] = p[:, 1]               # Neumann conditions
    p[:, ny - 1] = p[:, ny - 2]     # ... same as abovels



# Plot the pressure surface at the initial and final conditions
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(X, Y)

# ax.plot_surface(X, Y, Z0, rstride=1, cstride=1, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=True)

surf = ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
