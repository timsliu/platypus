import numpy as np


def get_subplot_config(subplots):
    '''return an arrangement of subplots given the total number of subplots
    to plot. The function returns the number of rows and columns needed to 
    accomodate all subplots, and tries to arrange them in a square matrix, 
    defaulting to more columns than rows.
    inputs: subplots - number of subplots
    outputs: rows - rows of subplots
             cols - columns of subplots'''

    # start with square array large enough to fit subplots
    rows = int(np.ceil(np.sqrt(subplots)))
    cols = rows

    # reduce excess rows
    while (rows - 1) * cols >= subplots:
        rows -= 1

    # try to reduce rows further for a fully filled rectangle
    if (rows - 1) * (cols + 1) == subplots:
        rows -= 1
        cols += 1

    return rows, cols

