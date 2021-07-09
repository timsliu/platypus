# utils.py
# miscellaneous utility functions

import pickle

def save_pickle(name, data):
    '''save some data to a pickle file'''
    file_name = name + ".p" 
    pickle.dump(data, open(file_name, "wb"))
    return

def points_to_area(p1, p2):
    '''calculates the area of a rectangle bounded by two points
    inputs: p1 - (x, y) tuple of coordinates for first point
            p2 - (x, y) tuple of coordinates for second point'''

    return abs((p1[0] - p2[0]) * (p1[1] - p2[1]))

def points_to_volume(p1, p2):
    '''calculates the volume of a rectangular prism bounded by two points
    inputs: p1 - (x, y) tuple of coordinates for first point
            p2 - (x, y) tuple of coordinates for second point'''

    return abs((p1[0] - p2[0]) * (p1[1] - p2[1]) * (p1[2] - p2[2]))
