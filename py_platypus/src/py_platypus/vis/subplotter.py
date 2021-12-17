"""
Classes for plotting a single subplot
"""
import numpy as np

class Subplotter:
    """
    Base class for a subplotter. Each subplotter has methods for plotting data
    and for determining consistent axes.
    """
    def __init__(self,
                 y_axis_zero,
                 x_axis_zero,
                 lim_factor=0.05):
        """
        inputs: data - array like data to be plotted; each element is the data
                for a single time series or line
                y_axis_zero - y axis should start at zero
                x_axis_zero - x axis should start at zero
                lim_factor - proportion the axis should exceed the data
        """
        self.data = None
        self.y_axis_zero = y_axis_zero
        self.x_axis_zero = x_axis_zero
        self.lim_factor = lim_factor
        self.lim_up = 1 + self.lim_factor
        self.lim_dn = 1 - self.lim_factor
        self.ylims = None
        self.xlims = None

    def set_data(self, data):
        self.data = data
        self.get_ylimits()

    def get_ylimits(self):
        """
        Default method for getting appropriate y limits for the plot given the
        data set
        """
        data_min = np.min(self.data)
        data_max = np.max(self.data)

        # rules for getting max limit
        if data_max > 0:
            lim_pos = self.lim_up * data_max
        else:
            lim_pos = 0.9 * data_max

        # rules for getting min limit
        if self.y_axis_zero:
            lim_neg = 0
        elif data_min < 0:
            lim_neg = self.lim_up * data_min
        else:
            lim_neg = self.lim_dn * data_min

        self.ylims = (lim_neg, lim_pos)
        return lim_neg, lim_pos

    def get_xlimits(self):
        """
        Get appropriate x limits for the plot given the data set
        """
        raise NotImplementedError

    def plot_axes(self, axs, data_idx):
        """
        Plot the data at index data_index of self.data on the passed matplotlib
        axes.
        """
        raise NotImplementedError


class SubplotLines(Subplotter):
    """
    Subplotter for when each subplot is a series of lines
    """
    def __init__(self,
                 y_axis_zero,
                 x_axis_zero,
                 lim_factor=0.05):

        super().__init__(y_axis_zero, x_axis_zero, lim_factor) 

    def get_xlimits(self):
       """
       Find the x limits of the data
       """

       # find maximum length of a series across all time steps
       lim_pos = max([max([len(series) for series in time_step]) for time_step in self.data])
       lim_neg = 0
       self.ylims = (lim_neg, lim_pos)
       return lim_neg, lim_pos


    def plot_axes(self, axs, data_idx):

        # iterate through the lines for this time step
        for i in range(len(self.data[data_idx])):
            plot_obj = axs.plot(self.data[data_idx][i])

        axs.set_ylim(self.ylims[0], self.ylims[1])

        return plot_obj[0]
