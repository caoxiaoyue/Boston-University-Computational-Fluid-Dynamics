import numpy as np
from boundary_conditions import *
import matplotlib.pyplot as plt

# Initial problem parameters
xDomain = (0.0,2.0*np.pi)     # x domain
nx = 70     # number of x-grid points
nt = 50     # number of time steps
dt = 0.01   # time step size
nu = 0.1    # diffusion coefficient
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, nt))

# Create an array of x values based on the discretization
x = np.arange(xDomain[0], xDomain[1]+dx, dx)

# Create the initial conditions.  The initial boundary conitions are taken care of
def phi(x,nu):
	return np.exp(-x**2/(4*nu)) + np.exp(-(x-2*np.pi)**2/(4*nu))

def dphi_dx(x,nu):
	return -x/(2*nu)*np.exp(-x**2/(4*nu)) - (2*x-4*np.pi)/(4*nu)*np.exp(-(x-2*np.pi)**2/(4*nu))

u[:,0] = -2*nu*dphi_dx(x,nu)/phi(x,nu) + 4
bc = (u[0,0], u[nx-1,0])

# March along time steps
for n in range(nt-1):
    u[0,n+1], u[nx-1,n+1] = set_boundary_conditions(bc)
    for i in range(1,nx-1):
        u[i, n+1] = u[i,n] - u[i,n]*dt/dx*(u[i,n] - u[i-1,n]) + nu*dt/(dx**2)*(u[i+1,n] - 2.0*u[i,n] + u[i-1,n])

# Plot the velocities at the end of the computation
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