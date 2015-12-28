from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from initial_conditions5_8 import *

# Initial problem parameters
xDomain = (0.0, 2.0)    # x domain
yDomain = (0.0, 2.0)    # y domain
nx = 20     # number of x-grid points
ny = 20     # number of y-grid points
nt = 50     # number of time steps
dt = 0.01   # time step size
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = float( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, ny, nt), dtype=np.float64)
v = np.zeros((nx, ny, nt), dtype=np.float64)

# Set the initial and boundary conditions
u[:, :, 0] = set_init_conditions(nx, dx, ny, dy)
v[:, :, 0] = set_init_conditions(nx, dx, ny, dy)
u[0, :, :] = 1.0
u[nx-1, :, :] = 1.0
u[:, 0, :] = 1.0
u[:, ny-1, :] = 1.0
v[0, :, :] = 1.0
v[nx-1, :, :] = 1.0
v[:, 0, :] = 1.0
v[:, ny-1, :] = 1.0

# March along time steps
for n in range(nt-1):
    un = u[:, :, n]
    vn = v[:, :, n]
    for i in range(1,nx):
        for j in range(1, ny):
            u[i, j, n+1] = un[i, j] - \
                           un[i, j]*dt/dx*(un[i, j] - un[i-1, j]) - \
                           vn[i, j]*dt/dy*(un[i, j] - un[i, j-1])
            v[i, j, n+1] = vn[i, j] - \
                           un[i, j]*dt/dx*(vn[i, j] - vn[i-1, j]) - \
                           vn[i, j]*dt/dy*(vn[i, j] - vn[i, j-1])

# Plot the vector plot of the velocity at the initial and final conditions
X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
X, Y = np.meshgrid(X, Y)
U0 = u[:, :, 0]
V0 = v[:, :, 0]

U = u[:, :, nt-1]
V = v[:, :, nt-1]

fig0, ax0 = plt.subplots()
strm0 = ax0.streamplot(X, Y, U0, V0, color=U, linewidth=2, cmap=plt.cm.autumn)
fig0.colorbar(strm0.lines)

fig1, ax1 = plt.subplots()
strm1 = ax1.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
fig1.colorbar(strm1.lines)

plt.show()