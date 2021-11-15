'''
math_utils.py

Miscellaneous math utility functions
'''

import copy
import numpy as np
from py_platypus.utils import constants


def points_to_area(p1, p2):
    '''
    calculates the area of a rectangle bounded by two points
    inputs: p1 - (x, y) tuple of coordinates for first point
            p2 - (x, y) tuple of coordinates for second point
    '''

    return abs((p1[0] - p2[0]) * (p1[1] - p2[1]))


def points_to_volume(p1, p2):
    '''
    calculates the volume of a rectangular prism bounded by two points
    inputs: p1 - (x, y) tuple of coordinates for first point
            p2 - (x, y) tuple of coordinates for second point
    '''

    return abs((p1[0] - p2[0]) * (p1[1] - p2[1]) * (p1[2] - p2[2]))


def union_lists(l1, l2):
    '''
    Returns the union of two lists without modifying the original lists.
    '''
    union = copy.deepcopy(l2)

    # iterate through elements in l1, and add them to the union set
    # if they don't appear in l2
    for element in l1:
        if element not in l2:
            union.append(element)

    return union


def wrap_idx_2d(array, i, j):
    '''
    Wrap indices to fit along a 2d array and return value at that index
    inputs: array - 2D numpy array
            i - row index
            j - column index
    '''
    i_wrapped = i % array.shape[0]
    j_wrapped = j % array.shape[1]
    return array[i_wrapped, j_wrapped]


def ev_to_vth(ev, mass=constants.ELECTRON_MASS):
    '''
    Convert from temperature of particles (in electron volts) to thermal
    velocity (in meters per second)
    '''

    joules = ev * constants.JOULES_PER_EV
    return np.sqrt(2 * joules / mass)
