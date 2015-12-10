import numpy as np
from initial_conditions import *
import matplotlib.pyplot as plt

# Initial problem parameters
xDomain = (0.0,2.0)     # x domain
nx = 20     # number of x-grid points
nt = 50     # number of time steps
dt = 0.01   # time step size
nu = 0.1    # diffusion coefficient
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, nt))

# Create the initial conditions
u[:,0] = set_init_conditions(nx, dx)

# March along time steps
for n in range(1,nt):
    for i in range(2,nx):
        u[i,n] = u[i,n-1] + nu*dt/(dx**2)*(u[i+1,n-1] - 2.0*u[i,n-1] + u[i-1,n-1])

# Plot the velocities at the end of the computation
x = np.arange(xDomain[0], xDomain[1]+dx, dx)

fig = plt.figure()
ax = plt.subplot(111)
plt.ylabel('velocity')

for it in range(0,nt,5):
    t = it * dt
    ax.plot(x, u[:,it], label="t="+str(t))

# Shrink current axis' height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
plt.show()