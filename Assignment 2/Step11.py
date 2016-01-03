from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from differences import *

# Initial problem parameters
xDomain = (0.0, 2.0)  # x domain
yDomain = (0.0, 2.0)  # y domain
nx = 50  # number of x-grid points
ny = 50  # number of y-grid points
nt = 50    # number of time steps
nit = 100   # number of iterations
nu = 0.1   # viscosity
rho = 1.0   # density
beta = np.float64(0.1)	    # beta = nu*dt/dx**2 or nu*dt/dy**2
dx = np.float64( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = np.float64( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y
dt = min(beta * (dx**2) / nu, beta * (dy ** 2) / nu)   # time step size

# Initial and Boundary Conditions
u = np.zeros((nx, ny, nt), dtype=np.float64)
v = np.zeros((nx, ny, nt), dtype=np.float64)
p = np.zeros((nx, ny, nt), dtype=np.float64)
pn = np.zeros((nx, ny), dtype=np.float64)
un = np.zeros((nx, ny), dtype=np.float64)
vn = np.zeros((nx, ny), dtype=np.float64)
b = np.zeros((nx, ny), dtype=np.float64)
speed = np.zeros((nx, ny, nt), dtype=np.float64)
u[:, ny-1, :] = 1.0

# Loop over time
for n in range(nt-1):

    for i in range(nx-1):
        for j in range(ny-1):
            b[i, j] = rho * ((firstDerCD(u[:, j, n], i, dx) + firstDerCD(v[i, :, n], j, dy)) / dt -
                             (firstDerCD(u[:, j, n], i, dx)) ** 2 -
                             2 * firstDerCD(u[i, :, n], j, dy) * firstDerCD(v[:, j, n], i, dx) -
                             (firstDerCD(v[i, :, n], j, dy)) ** 2)

    for iit in range(nit):
        pn = p
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                p[i, j] = (dy ** 2 * (pn[i + 1, j] + pn[i - 1, j]) +
                       dx ** 2 * (pn[i, j + 1] + pn[i, j - 1]) -
                       b[i,j] * dx**2 * dy**2) / (2 * (dx ** 2 + dy ** 2))
        p[:, ny-1, n] = 0.0
        p[0, :, n] = p[1, :, n]
        p[nx-1, :, n] = p[nx-2, :, n]
        p[:, 0, n] = p[:, 1, n]
#        p[:, ny-1, n] = p[:, ny-2, n]

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
            speed[i, j, n + 1] = math.sqrt(u[i, j, n + 1]**2 + v[i, j, n + 1]**2)
    u[0, :, n] = 0.0
    u[nx-1, :, n] = 0.0
    u[:, 0, n] = 0.0
    u[:, ny-1, n] = 1.0
    v[0, :, n] = 0.0
    v[nx-1, :, n] = 0.0
    v[:, 0, n] = 0.0
    v[:, ny-1, n] = 0.0

# Plots
X = np.linspace(xDomain[0], xDomain[1], nx)
Y = np.linspace(yDomain[0], yDomain[1], ny)
extents = (min(X), max(X), min(Y), max(Y))
extents2 = (min(X)-dx, max(X)+dx, min(Y)-dy, max(Y)+dy)
X, Y = np.meshgrid(X, Y)
P0 = p[:, :, 0]
P = p[:, :, nt-1]
U0 = u[:, :, 0]
V0 = v[:, :, 0]
U = u[:, :, nt-1]
V = v[:, :, nt-1]

# Plot pressure contours at the final time step
plt.figure()
plt.hot()
im = plt.imshow(P, interpolation='bilinear', origin='lower', cmap=cm.coolwarm, extent=extents2)
CS = plt.contour(X, Y, P, 9, origin='lower', antialiased = True, linewidths=2, extent=extents2)
CB = plt.colorbar(im, shrink=0.8, extend='both')
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Pressure Contour Plot')

plt.show()

# Plot velocity contours at the final time step
S = speed[:, :, nt-1]
plt.figure()
plt.hot()
im1 = plt.imshow(S, interpolation='bilinear', origin='lower', cmap=cm.coolwarm, extent=extents2)
CS1 = plt.contour(X, Y, S, 9, origin='lower', antialiased = True, linewidths=2, extent=extents2)
CB1 = plt.colorbar(im1, shrink=0.8, extend='both')
plt.clabel(CS1, inline=1, fontsize=10)
plt.title('Speed Contour Plot')

plt.show()

# Plot the vector plot of the velocity at the initial and final conditions
# fig0, ax0 = plt.subplots()
# strm0 = ax0.streamplot(X, Y, U0, V0, color=U, linewidth=2, cmap=plt.cm.coolwarm)
# fig0.colorbar(strm0.lines)
#
# fig1, ax1 = plt.subplots()
# strm1 = ax1.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.coolwarm)
# fig1.colorbar(strm1.lines)

plt.show()
