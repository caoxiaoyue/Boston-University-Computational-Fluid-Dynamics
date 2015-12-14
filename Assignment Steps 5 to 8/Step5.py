import numpy as np
from initial_conditions58 import *
import matplotlib.pyplot as plt

# Initial problem parameters
xDomain = (0.0, 2.0)    # x domain
yDomain = (0.0, 2.0)    # y domain
nx = 20     # number of x-grid points
ny = 20     # number of y-grid points
nt = 50     # number of time steps
dt = 0.01   # time step size
c = 1.0     # transport velocity
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dy = float( (yDomain[1]-yDomain[0])/(ny - 1) )  # delta y

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, ny, nt))

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
                           c*dt/dx*(un[i, j] - un[i-1, j]) - \
                           c*dt/dy*(un[i, j] - un[i, j-1])

# Plot the velocities at the end of the computation
x = np.arange(xDomain[0], xDomain[1]+dx, dx)
y = np.arange(yDomain[0], yDomain[1]+dy, dy)

# fig = plt.figure()
# ax = plt.subplot(111)
# plt.ylabel('velocity')
#
# for it in range(0,nt,5):
#     t = it * dt
#     ax.plot(x, u[:,it], label="t="+str(t))
#
# # Shrink current axis' height by 10% on the bottom
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
#
# ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
# plt.show()