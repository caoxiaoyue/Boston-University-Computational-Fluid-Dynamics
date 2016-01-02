from pylab import *
from boundary_conditions1_4 import *
from initial_conditions1_3 import *
from differences import *

# Initial problem parameters
xDomain = (0.0, 2.0)  # x domain
nx = 100  # number of x-grid points
nt = 50  # number of time steps
sigma = 0.25  # CFL number (sigma = c*dt/dx)
c = 1.0  # transport velocity
bc = (1, 1)  # Tuple of boundary conditions
dx = float((xDomain[1] - xDomain[0]) / (nx - 1))  # delta x
dt = sigma * dx / c  # time step size

# Create an empty array for all velocity time steps including t=0
u = zeros((nx, nt))

# Set the initial and boundary conditions
u[:, 0] = set_init_conditions(nx, dx)
u[0, :], u[nx - 1, :] = set_boundary_conditions(bc)

# March along time steps
for n in range(nt - 1):
    un = u[:, n]
    for i in range(1, nx):
        u[i, n + 1] = un[i] - c * dt * firstDerBD(un, i, dx)

# Plot the velocities at the end of the computation
x = arange(xDomain[0], xDomain[1] + dx, dx)

fig = figure()
ax = subplot(111)
ylabel('velocity')

for it in range(0, nt, 5):
    t = it * dt
    ax.plot(x, u[:, it], label="t=" + str(t))

# Shrink current axis' height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
show()
