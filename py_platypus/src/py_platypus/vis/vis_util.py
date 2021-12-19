import matplotlib.pyplot as plt
import pickle
import numpy as np

BINS = 20
LIM_FACTOR = 0.05

LIM_UP = 1 + LIM_FACTOR
LIM_DN = 1 - LIM_FACTOR


def get_ylimits(data, zero, subplotter):
    '''generate reasonable upper and lower axis limits for the data
    inputs: data - array of all data
            zero - bool for setting the lower limit to zero'''

    data_min = np.min(data)
    data_max = np.max(data)

    # special case for 1D histogram
    if subplotter == subplot_histogram:
        lim_neg = 0
        lim_pos = LIM_UP * max(
            [np.max(np.histogram(x, bins=BINS)[0]) for x in data])
        return lim_neg, lim_pos

    # special case for 2D histogram
    if subplotter == subplot_histogram_2d:
        lim_neg = 0
        lim_pos = LIM_UP * max(
            [np.max(np.histogram2d(x[0], x[1], bins=BINS)[0]) for x in data])

        return lim_neg, lim_pos

    if subplotter == subplot_scatter_2d:
        data_max = max([np.max(x[1]) for x in data])
        data_min = min([np.min(x[1]) for x in data])

    if subplotter == subplot_scatter_3d:
        data_max = max([np.max(x[2]) for x in data])
        data_min = min([np.min(x[2]) for x in data])

    # rules for getting max limit
    if data_max > 0:
        lim_pos = LIM_UP * data_max
    else:
        lim_pos = 0.9 * data_max

    # rules for getting min limit
    if zero:
        lim_neg = 0
    elif data_min < 0:
        lim_neg = LIM_UP * data_min
    else:
        lim_neg = LIM_DN * data_min

    return lim_neg, lim_pos


def get_xlimits(data, zero, subplotter):
    """
    Generate reasonable upper and lower x-axis limits for a given
    data set.
    inputs: data - np array of data to plot
            zero - whether the lower value should be fixed at zero
            subplotter - function used for plotting
    """

    if subplotter == subplot_scatter_3d or subplotter == subplot_scatter_2d:
        data_max = max([np.max(x[0]) for x in data])
        data_min = min([np.min(x[0]) for x in data])

    elif subplotter == subplot_grid:
        data_max = np.array(data).shape[
            2]  # data shape dimensions are time, rows, cols
        data_min = 0

    else:
        raise ValueError(
            "Subplot type {} not supported by get_xlimits".format(subplotter))

    # rules for getting max limit
    if data_max > 0:
        lim_pos = LIM_UP * data_max
    else:
        lim_pos = 0.9 * data_max

    # rules for getting min limit
    if zero:
        lim_neg = 0
    elif data_min < 0:
        lim_neg = LIM_UP * data_min
    else:
        lim_neg = LIM_DN * data_min

    return lim_neg, lim_pos


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


def subplot_lines(axs, data, legend, lims):
    '''
    DEPRECATED 
    helper function that plots multiple lines on a subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - list of data to plot
            legend - list of labels in the same order as the lines
            lims - tuple specifying upper and lower bounds'''

    plot_lines = len(data)  # number of lines to plot

    # iterate through the lines
    for j in range(plot_lines):

        # plot with or without lines
        if legend is not None:
            plot_obj = axs.plot(data[j], label=legend[j])
        else:
            plot_obj = axs.plot(data[j])

    axs.set_ylim(lims[0], lims[1])

    return plot_obj


def subplot_grid(axs, data, legend, lims):
    '''helper function for plotting values that lie on a 2-D grid.
    inputs: axs - matplotlib axes object for a single subplot
            data - 2D array of data to plot
            legend - unused paramter'''

    if len(data) > 1:
        raise ValueError("subplot_grid mishapened input")
    plot_obj = axs.imshow(data[0],
                          interpolation='none',
                          vmin=lims[0],
                          vmax=lims[1])

    return plot_obj


def subplot_scatter_2d(axs, data, legend, lims):
    '''helper function for plotting a scatter plot of data with 2
    dimensions.
    inputs: axs - matplotlib axes object for a single subplot
            data - list with each element containing 2 arrays 
                   [[x_vals0, y_vals0], [x_vals1, y_vals1]]
            legend - unused paramter'''

    # plot with or without lines
    if legend is not None:
        plot_obj = axs.scatter(data[0], data[1], label=legend, s=0.1)
    else:
        plot_obj = axs.scatter(data[0], data[1], s=1)

    axs.set_ylim(lims[0], lims[1])

    return plot_obj


def subplot_scatter_3d(axs, data, legend, lims):
    '''helper function for plotting a scatter plot of data with 3
    dimensions. The third dimension is illustrated with color. Only a single
    data series can be plotted per subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - list with 3 arrays [x_vals0, y_vals0, z_vals0]
            legend - unused paramter'''

    plot_obj = axs.scatter(data[0],
                           data[1],
                           c=data[2],
                           s=0.1,
                           vmin=lims[0],
                           vmax=lims[1])

    return plot_obj


def subplot_histogram(axs, data, legend, lims):
    '''
    DEPRECATED
    helper function for plotting a histogram subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - 1D array of data
            legend - unused paramter'''

    plot_obj = axs.hist(data[0], bins=BINS)
    axs.set_ylim(0, lims[1])

    return plot_obj


def subplot_histogram_2d(axs, data, legend, lims):
    '''helper function for plotting 2D data as a histogram subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - array of data [x_values, y_values]
            legend - unused paramter'''

    plot_obj = axs.hist2d(data[0],
                          data[1],
                          bins=BINS,
                          vmin=lims[0],
                          vmax=lims[1])

    return plot_obj
