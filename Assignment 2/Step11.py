from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from differences import *

# Initial problem parameters
xDomain = (0.0, 2.0)  # x domain
yDomain = (0.0, 2.0)  # y domain
nx = 20  # number of x-grid points
ny = 20  # number of y-grid points
nt = 50    # number of time steps
nit = 100   # number of iterations
nu = 0.1   # viscosity
rho = 1.0   # density
beta = np.float64(0.1)	    # beta = nu*dt/dx**2 or nu*dt/dy**2
dx = np.float64( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = np.float64( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y
#dt = min(beta * (dx**2) / nu, beta * (dy ** 2) / nu)   # time step size
dt = 0.01

# Initial and Boundary Conditions
u = np.zeros((nx, ny, nt), dtype=np.float64)
v = np.zeros((nx, ny, nt), dtype=np.float64)
p = np.zeros((nx, ny, nt), dtype=np.float64)
pn = np.zeros((nx, ny), dtype=np.float64)
un = np.zeros((nx, ny), dtype=np.float64)
vn = np.zeros((nx, ny), dtype=np.float64)
b = np.zeros((nx, ny), dtype=np.float64)
u[:, ny-1, :] = 1.0

# Loop over time
for n in range(nt-1):

    for i in range(nx-1):
        for j in range(ny-1):
            b[i, j] = rho * ((firstDerCD(u[:, j, n], i, dx) + firstDerCD(v[i, :, n], j, dy)) / dt +
                             (firstDerCD(u[:, j, n], i, dx)) ** 2 +
                             2 * firstDerCD(u[i, :, n], j, dy) * firstDerCD(v[:, j, n], i, dx) +
                             (firstDerCD(v[i, :, n], j, dy)) ** 2)

    for iit in range(nit):
        pn = p
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                p[i, j] = (dy ** 2 * (pn[i + 1, j] + pn[i - 1, j]) +
                       dx ** 2 * (pn[i, j + 1] + pn[i, j - 1]) -
                       b[i,j] * dx**2 * dy**2) / (2 * (dx ** 2 + dy ** 2))
        p[0, :, n] = p[1, :, n]
        p[nx-1, :, n] = p[nx-2, :, n]
        p[:, 0, n] = p[:, 1, n]
        p[:, ny-1, n] = p[:, ny-2, n]

    pn = p[:, :, n]
    un = u[:, :, n]
    vn = v[:, :, n]
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            u[i, j, n + 1] = un[i, j] - \
                            un[i, j] * dt * firstDerBD(un[:, j], i, dx) - \
                            vn[i, j] * dt * firstDerBD(un[i, :], j, dy) - \
                            1.0/rho * dt * firstDerCD(pn[:, j], i, dx) + \
                             nu * dt * secDerCD(un[:, j], 1, dx) + \
                             nu * dt * secDerCD(un[i, :], j, dy)
            v[i, j, n + 1] = vn[i, j] - \
                            un[i, j] * dt * firstDerBD(vn[:, j], i, dx) - \
                            vn[i, j] * dt * firstDerBD(vn[i, :], j, dy) - \
                            1.0/rho * dt * firstDerCD(pn[i, :], j, dy) + \
                             nu * dt * secDerCD(vn[:, j], 1, dx) + \
                             nu * dt * secDerCD(vn[i, :], j, dy)
    u[0, :, n] = 0.0
    u[nx-1, :, n] = 0.0
    u[:, 0, n] = 0.0
    u[:, ny-1, n] = 1.0
    v[0, :, n] = 0.0
    v[nx-1, :, n] = 0.0
    v[:, 0, n] = 0.0
    v[:, ny-1, n] = 0.0

# Plot the vector plot of the velocity at the initial and final conditions
# Plot the vector plot of the velocity at the initial and final conditions
X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
X, Y = np.meshgrid(X, Y)
U0 = u[:, :, 0]
V0 = v[:, :, 0]

U = u[:, :, nt-1]
V = v[:, :, nt-1]

print U

fig0, ax0 = plt.subplots()
strm0 = ax0.streamplot(X, Y, U0, V0, color=U, linewidth=2, cmap=plt.cm.coolwarm)
fig0.colorbar(strm0.lines)

fig1, ax1 = plt.subplots()
strm1 = ax1.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.coolwarm)
fig1.colorbar(strm1.lines)

# Plot the pressure surface at the initial and final conditions
fig2 = plt.figure()
ax = fig2.gca(projection='3d')

P0 = p[:, :, 0]
P = p[:, :, nt-1]

ax.plot_surface(X, Y, P0, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=True)

surf = ax.plot_surface(X, Y, P, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=True)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig2.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
