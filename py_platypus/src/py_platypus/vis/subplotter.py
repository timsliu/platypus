"""
Classes for plotting a single subplot
"""
import numpy as np


class Subplotter:
    """
    Base class for a subplotter. Each subplotter has methods for plotting data
    and for determining consistent axes.
    """
    def __init__(self, y_axis_zero, x_axis_zero, lim_factor=0.05):
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
        self.y_lims = None
        self.x_lims = None

    def set_data(self, data):
        """
        Set the data that will be used for the subplot
        """
        self.data = data

    def get_limits(self, data_min, data_max, axis_zero):
        """
        Get the axis limits from the minimum and maximum data values
        """
        # rules for getting max limit
        if data_max > 0:
            upper_lim = self.lim_up * data_max
        else:
            upper_lim = self.lim_dn * data_max

        # rules for getting min limit
        if axis_zero:
            lower_lim = 0
        elif data_min < 0:
            lower_lim = self.lim_up * data_min
        else:
            lower_lim = self.lim_dn * data_min

        return lower_lim, upper_lim

    def get_y_limits(self):
        """
        Default method for getting appropriate y limits for the plot given the
        data set
        """
        if self.y_lims is None:
            self.y_lims = self.get_limits(*self.get_y_min_max(),
                                          self.y_axis_zero)
        return self.y_lims

    def get_x_limits(self):
        """
        Get appropriate x limits for the plot given the data set
        """
        if self.x_lims is None:
            self.x_lims = self.get_limits(*self.get_x_min_max(),
                                          self.x_axis_zero)
        return self.x_lims

    def get_y_min_max(self):
        """
        Abstract method for getting the min and max y values in the data.
        """
        raise NotImplementedError

    def get_x_min_max(self):
        """
        Abstract method for getting the min and max x values in the data.
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
    Subplotter for a series of lines
    """
    def __init__(self, y_axis_zero, x_axis_zero, lim_factor=0.05):

        super().__init__(y_axis_zero, x_axis_zero, lim_factor)

    def get_y_min_max(self):
        """
        Find the minimum and maximum y values of the data
        """
        return np.min(self.data), np.max(self.data)

    def get_x_min_max(self):
        """
        Find the minimum and maximum x values of the data
        """

        x_max = max([
            max([len(series) for series in time_step])
            for time_step in self.data
        ])
        x_min = min([
            min([len(series) for series in time_step])
            for time_step in self.data
        ])
        return x_min, x_max

    def plot_axes(self, axs, data_idx):

        # iterate through the lines for this time step
        for i in range(len(self.data[data_idx])):
            plot_obj = axs.plot(self.data[data_idx][i])

        return plot_obj[0]


class SubplotScatter2D(Subplotter):
    """
    Subplotter for a series of 2d scatter plots
    """
    def __init__(self, y_axis_zero, x_axis_zero, lim_factor=0.05):

        super().__init__(y_axis_zero, x_axis_zero, lim_factor)

    def get_x_min_max(self):
        """
        Find the min and max x values of the data
        """

        x_max = max([np.max(x[0]) for x in self.data])
        x_min = min([np.min(x[0]) for x in self.data])
        return x_min, x_max

    def get_y_min_max(self):
        """
        Find the min and max y values of the data
        """

        y_max = max([np.max(x[1]) for x in self.data])
        y_min = min([np.min(x[1]) for x in self.data])
        return y_min, y_max

    def plot_axes(self, axs, data_idx):
        """
        Plot a 2D scatter plot
        """
        plot_obj = axs.scatter(self.data[data_idx][0], self.data[data_idx][1])
        return plot_obj


class SubplotHistogram(Subplotter):
    """
    Subplotter for a series of 1D histograms
    """
    def __init__(self, x_axis_zero=False, y_axis_zero=True, lim_factor=0.05):
        self.num_bins = 20
        self.bin_edges = None
        super().__init__(y_axis_zero, x_axis_zero, lim_factor)

    def get_x_min_max(self):
        """
        Find the min and max x values of the data
        """

        # compute the bins that would be used if all data were flattened
        # and put into a single histogram
        bins = np.histogram_bin_edges(self.data, bins=self.num_bins)
        x_min = min(bins)
        x_max = max(bins)
        self.bin_edges = bins
        return x_min, x_max

    def get_y_min_max(self):
        """
        Find the min and max y values of the data
        """
        y_max = max([
            np.max(np.histogram(time_step, bins=self.num_bins)[0])
            for time_step in self.data
        ])

        y_min = 0
        return y_min, y_max

    def plot_axes(self, axs, data_idx):
        """
        Plot a 1D histogram
        """
        if self.bin_edges is None:
            self.get_x_min_max()
        vals, bins, bar_container = axs.hist(self.data[data_idx],
                                             bins=self.bin_edges)
        return bar_container.patches[0]
