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
nx = 20     # number of x-grid points
ny = 20     # number of y-grid points
nt = 50     # number of time steps
beta = np.float64(0.1)	    # beta = nu*dt/dx**2 or nu*dt/dy**2
nu = 0.1    # viscosity
dx = np.float64( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = np.float64( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y
dt = np.minimum(beta*(dx**2)/nu, beta*(dy**2)/nu)   # time step size

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, ny, nt), dtype=np.float64)
v = np.zeros((nx, ny, nt), dtype=np.float64)

# Set the initial and boundary conditions
u[:, :, 0] = set_init_conditions(nx, dx, ny, dy)
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
    un = u[ :, :, n]
    vn = v[:, :, n]
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            u[i, j, n+1] = un[i, j] - \
                           un[i, j] * dt * firstDerBD(un[:, j], i, dx) - \
                           vn[i, j] * dt * firstDerBD(un[i, :], j, dy) +\
                           nu * dt * secDerCD(un[:, j], i, dx) + \
                           nu * dt * secDerCD(un[i, :], j, dy)
            v[i, j, n+1] = vn[i, j] - \
                           un[i, j] * dt * firstDerBD(vn[:, j], i, dx) - \
                           vn[i, j] * dt * firstDerBD(vn[i, :], j, dy) + \
                           nu * dt * secDerCD(vn[:, j], i, dx) + \
                           nu * dt * secDerCD(vn[i, :], j, dy)


# Plot the vector plot of the velocity at the initial and final conditions
X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
X, Y = np.meshgrid(X, Y)
U0 = u[:, :, 0]
V0 = v[:, :, 0]

U = u[:, :, nt-1]
V = v[:, :, nt-1]

print X.shape, Y.shape
print U0.shape, V0.shape
print U.shape, V.shape

fig0, ax0 = plt.subplots()
strm0 = ax0.streamplot(X, Y, U0, V0, color=U, linewidth=2, cmap=plt.cm.coolwarm)
fig0.colorbar(strm0.lines)

fig1, ax1 = plt.subplots()
strm1 = ax1.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.coolwarm)
fig1.colorbar(strm1.lines)

plt.show()