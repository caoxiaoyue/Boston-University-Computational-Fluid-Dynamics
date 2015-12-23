import numpy as np
from boundary_conditions import *
import matplotlib.pyplot as plt

# Initial problem parameters
xDomain = (0.0,2.0*np.pi)     # x domain
nx = 70     # number of x-grid points
nt = 50     # number of time steps
beta = 0.25	# beta = nu*dt/dx**2
nu = 0.1    # diffusion coefficient
dx = float( (xDomain[1]-xDomain[0])/(nx - 1) )  # delta x
dt = beta*(dx**2)/nu   # time step size

# Create an empty array for all velocity time steps including t=0
u = np.zeros((nx, nt))

# Create an array of x values based on the discretization
x = np.arange(xDomain[0], xDomain[1]+dx, dx)

# Create auxiliary function for periodic boundary conditions
ip1 = np.zeros(nx)
im1 = np.zeros(nx)
for i in range(nx):
    ip1[i] = i+1
    im1[i] = i-1
ip1[nx-1] = 1
im1[0] = nx-1

#print ip1
#print im1

# Create the initial conditions.  The initial boundary conitions are taken care of
def phi(x):
	return np.exp(-x**2/(4*nu)) + np.exp(-(x-2*np.pi)**2/(4*nu))

def dphi_dx(x):
	return -x/(2*nu)*np.exp(-x**2/(4*nu)) - (2*x-4*np.pi)/(4*nu)*np.exp(-(x-2*np.pi)**2/(4*nu))

u[:,0] = -2*nu*dphi_dx(x)/phi(x) + 4
bc = (u[0,0], u[nx-1,0])
u[0,:] = u[0,0]
u[nx-1,:] = u[0,:]

# March along time steps
for n in range(nt-1):
    un = u[:,n]
    #u[0,n+1], u[nx-1,n+1] = set_boundary_conditions(bc)
    for i in range(0,nx):
        u[i, n+1] = un[i] - un[i]*dt/dx*(un[i] - un[im1[i]]) +\
                    beta*(un[ip1[i]] - 2.0*un[i] + un[im1[i]])

# Analytical solution
def phiAnalytical(x,t):
    return np.exp(-(x-4*t)**2/(4*nu*(t+1))) + np.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1)))

def dphi_dx_Analytical(x,t):
    return -1.0/(2*nu*(t+1))*((x-4*t)*np.exp(-(x-4*t)**2/(4*nu*(t+1))) +
                              (x-4*t-2*np.pi)*np.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1))))

def uAnalytical(x,t):
    return -2*nu*dphi_dx_Analytical(x,t)/phiAnalytical(x,t) + 4.0

# Plot the velocities at the end of the computation
fig = plt.figure()
ax = plt.subplot(111)
plt.ylabel('velocity')

for it in range(0,nt,5):
    time = it * dt
    ax.plot(x, u[:,it], label="t="+str(time))

# Plot the analytical solution
ax.plot(x,uAnalytical(x, nt*dt), label="Analytical")

# Shrink current axis' height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
plt.show()