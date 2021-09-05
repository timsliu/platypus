'''
math_utils.py

Miscellaneous math utility functions
'''

import copy

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
