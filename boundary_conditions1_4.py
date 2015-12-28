import numpy as np


def set_boundary_conditions(bc):
    """
    Sets the boundary conditions of the domain,
    :param bc: tuple of values at each point of the domain
    :return: returns an numpy array with the boundary conditions
    """

    # Create an array of zeros
    val = np.zeros(2)

    val[0] = bc[0]
    val[1] = bc[1]

    return val