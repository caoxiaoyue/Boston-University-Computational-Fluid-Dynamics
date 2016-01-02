import numpy as np


def set_init_conditions(nx, dx):
    """
    Sets the initial conditions of the velocity at time t=0,
    :param nx: the number of grid points on the x-axis
    :param dx: the x grid spacing, delta x
    :return: returns an numpy array with the initial velocity
    """

    # Create an array of zeros
    u = np.zeros(nx, dtype = np.float64)

    # Loop over all x grid points and set the initial conditions
    for i in range(nx):
        x = i*dx
        if 0.5 <= x <= 1.0 :
            u[i] = 2.0
        else:
            u[i] = 1.0

    return u