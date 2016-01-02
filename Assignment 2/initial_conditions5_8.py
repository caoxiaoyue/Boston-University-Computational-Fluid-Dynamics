import numpy as np


def set_init_conditions(nx, dx, ny, dy):
    """
    Sets the initial conditions of the velocity at time t=0,
    :param nx: the number of grid points on the x-axis
    :param dx: the x grid spacing, delta x
    :param ny: the number of grid points on the y-axis
    :param dy: the y grid spacing, delta y
    :return: returns an numpy array with the initial velocity
    """

    # Create an array of zeros
    u = np.zeros((nx, ny), dtype=np.float64)

    # Loop over all x grid points and set the initial conditions
    for i in range(nx):
        for j in range(ny):
            x = i * dx
            y = j * dy
            if 0.5 <= x <= 1.0 and 0.5 <= y <= 1.0:
                u[i, j] = 2.0
            else:
                u[i, j] = 1.0

    return u