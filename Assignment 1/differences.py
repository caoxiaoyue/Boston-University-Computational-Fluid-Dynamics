# This file has all of the differences that discretize derivates of first and second order

from pylab import *


def firstDerBD(vec, index, delta):
    """ Returns the backward difference approximation to the first derivative
    :param vec: data vector
    :param index: index about which the difference calculations pivots
    :param delta: grid spacing
    :return: returns a number
    """
    return (vec[index] - vec[index - 1]) / delta

def firstDerFD(vec, index, delta):
    """ Returns the forward difference approximation to the first derivative
    :param vec: data vector
    :param index: index about which the difference calculations pivots
    :param delta: grid spacing
    :return: returns a number
    """
    return (vec[index + 1] - vec[index]) / delta

def firstDerCD(vec, index, delta):
    """ Returns the central difference approximation to the first derivative
    :param vec: data vector
    :param index: index about which the difference calculations pivots
    :param delta: grid spacing
    :return: returns a number
    """
    return (vec[index + 1] - vec[index - 1]) / (2 * delta)

def secDerCD(vec, index, delta):
    """ Returns the central difference approximation to the second derivative
    :param vec: data vector
    :param index: index about which the difference calculations pivots
    :param delta: grid spacing
    :return: returns a number
    """
    return (vec[index + 1] - 2 * vec[index] + vec[index - 1]) / (delta**2)
