"""
Helper classes for updating current density in 2D EM simulations
"""
import numpy as np
import py_platypus as plat

####################
# Helper functions
####################


def delta_midpoint(x0, dx, delta_x):
    '''
    calculates displacement to get from x0 to the nearest midpoint of
    the cell in the direction the particle is traveling. The arguments are
    named for the x direction, but this method
    can be used for motion in direction in any direction.
    inputs: x0 - starting coordinate
            dx - size of cell along this dimension
            delta_x - total displacement for this step
    '''

    # wrap the coordinate so that it is inside a single cell
    rel_x0 = x0 % dx

    # particle motion in the positive direction
    if np.sign(delta_x) == 1:
        # particle starts left half of the cells - nearest midpoint in cell
        if rel_x0 < dx / 2:
            delta_x0 = 0.5 * dx - rel_x0
        # particle starts right half of cell - nearest midpoint is next cell
        else:
            delta_x0 = 1.5 * dx - rel_x0
    # particle motion in the negative direction
    else:
        # particle starts in lower half of cell - nearest midpoint in
        # previous cell
        if rel_x0 < dx / 2:
            delta_x0 = -0.5 * dx - rel_x0
        # particle starts in upper half of cell - nearest midpoint is in
        # this cell
        else:
            delta_x0 = 0.5 * dx - rel_x0

    return delta_x0




####################
# Class definitions
####################


class ChargeStep:
    """
    Class representing a charge moving
    """
    def __init__(self, x0, y0, dx, dy):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x0 + dx
        self.y1 = y0 + dy
        self.dx = dx
        self.dy = dy

    def __eq__(self, c1):
        '''
        Equality operator, mostly used for testing.
        '''
        # returns whether all corresponding elements are close to each other
        return np.logical_and.reduce(
            np.isclose([self.x0, self.y0, self.dx, self.dy],
                       [c1.x0, c1.y0, c1.dx, c1.dy]))

    def __repr__(self):
        '''
        String representing the charge step
        '''

        return "Charge step ({}, {})->({}, {})".format(self.x0, self.y0,
                                                       self.x1, self.y1)


class ChargeStepDivider:
    """
    Helper class that breaks the motion of a charge into a list of 
    ChargeSteps such that each ChargeStep only crosses 4 grid boundaries
    """
    def __init__(self, dx, cells):
        """
        initializer class
        inputs: dx - list of the sizes of each grid cell [dy, dx]
                cells - list of the number of cells in each dimension
                [rows, cols]
        """

        self.dx = dx[1]
        self.dy = dx[0]
        self.cells = cells
        self.xmax = self.dx * cells[1]
        self.ymax = self.dy * cells[0]

    def get_charge_steps(self, x0, y0, x1, y1):
        """
        Takes the initial and final positions of a particle and returns
        a list of ChargeSteps such that each ChargeStep only crosses four
        grid boundaries
        """
        # convert coordinates to extended, non-periodic grid
        x0_ext, x1_ext, y0_ext, y1_ext = self.extend_coords(x0, x1, y0, y1)

        # list of boundaries touched at the starting and ending positions
        hori_0, vert_0 = self.boundaries_touched(x0_ext, y0_ext)
        hori_1, vert_1 = self.boundaries_touched(x1_ext, y1_ext)

        # list of unique boundaries touched
        hori_all = plat.math_utils.union_lists(hori_0, hori_1)
        vert_all = plat.math_utils.union_lists(vert_0, vert_1)

        self.check_boundaries_valid(vert_all)
        self.check_boundaries_valid(hori_all)

        # count unique boundaries touched at start and stop positions
        boundaries_touched = len(hori_all) + len(vert_all)

        if boundaries_touched == 4:
            # get the list of four boundary crossing motions
            substeps = self.get_submotions_four(x0_ext, y0_ext, x1_ext, y1_ext)
        elif boundaries_touched == 7:
            # get the list of four boundary crossing motions
            substeps = self.get_submotions_seven(x0_ext, y0_ext, x1_ext,
                                                 y1_ext, len(hori_all))
        elif boundaries_touched == 8:
            # touching 8 unique boundaries at start and stop positions
            # corresponds to crossing 10 boundaries
            substeps = self.get_submotions_ten(x0_ext, y0_ext, x1_ext, y1_ext)
        else:
            raise ValueError(
                "Warning! {} boundaries touched, expected 4, 7, or 8")

        return substeps

    def get_submotions_four(self, x0, y0, x1, y1):
        '''
        returns a list of charge steps for the four boundary case
        '''
        delta_x = x1 - x0
        delta_y = y1 - y0

        return [ChargeStep(x0, y0, delta_x, delta_y)]

    def get_submotions_seven(self, x0, y0, x1, y1, num_hori):
        '''
        returns a list of charge steps for the seven boundary case
        '''
        # determine the direction leading to two steps
        delta_x = x1 - x0
        delta_y = y1 - y0

        if num_hori == 3:
            # four vertical boundaries, three horizontal boundaries
            # sub motions dependent on x

            # moving toward or away from the midpoint
            delta_x0 = delta_midpoint(x0, self.dx, delta_x)
            delta_y0 = delta_x0 * (delta_y / delta_x)

        else:
            # four horizontal boundaries, three vertical boundaries
            # sub motions depdendent on y
            delta_y0 = delta_midpoint(y0, self.dy, delta_y)
            delta_x0 = delta_y0 * (delta_x / delta_y)

        delta_x1 = delta_x - delta_x0
        delta_y1 = delta_y - delta_y0

        c0 = ChargeStep(x0, y0, delta_x0, delta_y0)
        c1 = ChargeStep(x0 + delta_x0, y0 + delta_y0, delta_x1, delta_y1)

        return [c0, c1]

    def get_submotions_ten(self, x0, y0, x1, y1):
        '''
        returns a list of charge steps for the ten boundary case
        '''
        delta_x = x1 - x0
        delta_y = y1 - y0

        # fraction of the move to execute to get to the midpoint of a cell
        # along either dimension
        x_frac = delta_midpoint(x0, self.dx, delta_x) / delta_x
        y_frac = delta_midpoint(y0, self.dy, delta_y) / delta_y

        # shorter to get to the x midpoint - movement to x midpoint is first
        if x_frac < y_frac:
            delta_x0 = delta_x * x_frac
            delta_y0 = delta_y * x_frac

            delta_x1 = delta_x * y_frac - delta_x0
            delta_y1 = delta_y * y_frac - delta_y0

        # shorter to get to the y midpoint - movement to y midpoint is first
        else:
            delta_x0 = delta_x * y_frac
            delta_y0 = delta_y * y_frac

            delta_x1 = delta_x * x_frac - delta_x0
            delta_y1 = delta_y * x_frac - delta_y0

        delta_x2 = delta_x - delta_x1 - delta_x0
        delta_y2 = delta_y - delta_y1 - delta_y0

        c0 = ChargeStep(x0, y0, delta_x0, delta_y0)
        c1 = ChargeStep(c0.x1, c0.y1, delta_x1, delta_y1)
        c2 = ChargeStep(c1.x1, c1.y1, delta_x2, delta_y2)

        return [c0, c1, c2]

    def boundaries_touched(self, x, y):
        '''
        returns a list of boundaries touched
        inputs: x, y - coordinates of the particle
        returns: horizontal - list of indices of horizontal edges touched
                 vertical - list of indices of vertical edges touched
        '''
        # column index of the vertical boundary
        col_vert_idx = int(np.ceil((x - self.dx / 2) / self.dx))
        # row index of the horizontal boundary
        row_hori_idx = int(np.ceil((y - self.dy / 2) / self.dy))

        horizontal = [[row_hori_idx, col_vert_idx - 1],
                      [row_hori_idx, col_vert_idx]]
        vertical = [[row_hori_idx - 1, col_vert_idx],
                    [row_hori_idx, col_vert_idx]]

        return horizontal, vertical

    def extend_coords(self, x0, x1, y0, y1):
        """
        Take the coordinates from a wrapped grid space and convert them so the
        values are as if the grid was not wrapped from (0, 0) to (xmax, ymax).
        Only the final coordinates are modified
        inputs: x0 - starting x
                y0 - starting y
                x1 - starting x
                y1 - starting y
        """

        delta_x = x1 - x0
        delta_y = y1 - y0

        # particle transitioned from left side to right side
        if delta_x > 0.5 * self.xmax:
            x1 -= self.xmax
        # particle transitioned from right side to left side
        if delta_x < -0.5 * self.xmax:
            x1 += self.xmax
        # particle transitioned from low y to high y
        if delta_y > 0.5 * self.ymax:
            y1 -= self.ymax
        # particle transitioned from high y to low y
        if delta_y < -0.5 * self.ymax:
            y1 += self.ymax

        return x0, x1, y0, y1

    def check_boundaries_valid(self, boundaries):
        """
        Check that the collection of boundaries crossed is valid. If four or
        more unique boundaries are crossed along a single dimension, then
        the Courant-Friedrichs-Lewy condition has been violated.
        """

        x_coors = [boundary[0] for boundary in boundaries]
        y_coors = [boundary[1] for boundary in boundaries]

        if (max(x_coors) - min(x_coors) >= 3) or (max(y_coors) - min(y_coors)
                                                  >= 3):
            raise RuntimeError("{} - {}\n{} \nBoundaries: {}".format(
                "More than 10 boundaries touched",
                "Courant condition likely violated!",
                "Consider reducing timestep.", boundaries))
