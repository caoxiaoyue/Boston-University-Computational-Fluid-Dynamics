from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from initial_conditions5_8 import *
from differences import *

# Initial problem parameters
xDomain = (0.0, 2.0)    # x domain
yDomain = (0.0, 2.0)    # y domain
nx = 100     # number of x-grid points
ny = 100     # number of y-grid points
nt = 50     # number of time steps
sigma = 0.25  # CFL number (sigma = c*dt/dx or sigma = c*dt/dy)
c = 1.0     # transport velocity
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = float( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y
dt = np.minimum(sigma * dx / c, sigma * dy / c)  # time step size

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, ny, nt), dtype=np.float64)

# Set the initial and boundary conditions
u[:, :, 0] = set_init_conditions(nx, dx, ny, dy)
u[0, :, :] = 1.0
u[nx-1, :, :] = 1.0
u[:, 0, :] = 1.0
u[:, ny-1, :] = 1.0

# March along time steps
for n in range(nt-1):
    un = u[ :, :, n]
    for i in range(1,nx):
        for j in range(1, ny):
            u[i, j, n+1] = un[i, j] - \
                           c*dt*firstDerBD(un[:, j], i, dx) - \
                           c*dt*firstDerBD(un[i, :], j, dy)

# Plot the velocity surface at the initial and final conditions
fig = plt.figure()
ax = fig.gca(projection='3d')

X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
X, Y = np.meshgrid(X, Y)
Z0 = u[:,:,0]
Z = u[:,:,nt-1]

ax.plot_surface(X, Y, Z0, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()