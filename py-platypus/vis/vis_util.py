import matplotlib.pyplot as plt
import pickle
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


def plot_lines(filename, data, x_axis, y_axis, title, 
               subplotter, log=False, legend=None, steps=None):
    '''plot a single line chart with several lines
    inputs: filename - full path to output filename
            data - 2d array of data; dimension 0 is for each subplot, 
                   dimension 1 is for multiple lines on a subplot
            x-axis - name of x-axis
            y-axis - name of y-axis
            titles - list of chart titles
            log - display y as log plot
            legend - list of strings labeling the data'''

    subplots = len(data)                       # number of subplots
    rows, cols, = get_subplot_config(subplots) # subplot arrangement
    fig, axs = plt.subplots(rows, cols, constrained_layout=True, squeeze=False)
    
    lim_neg = 1.1 * np.min(data) 
    lim_pos = 1.1 * np.max(data)

    if legend is None:
        legend = subplots * [None]

    # iterate through subplots
    for i in range(rows * cols):
        # location in the subplot grid
        col = i % cols
        row = int(np.floor(i/cols))
      
        # remove axis for unused subplots 
        if i >= subplots:
            axs[row, col].set_axis_off()
            continue

        # call function to plot the data
        axs[row, col] = subplotter(axs[row, col], data[i], legend[i], (lim_neg, lim_pos))
        
        # add subplot title if available
        if steps is not None: 
            axs[row, col].set_title("Step {}".format(steps[i]), fontsize=10)
        
        # add optional features
        if legend[i] is not None:
            axs[row, col].legend()
        
        if log: 
            axs[row, col].yscale("log")
        
        axs[row, col].grid(True)

    # set the axis labels
    for ax in axs.flat:
        ax.set_xlabel(x_axis, fontsize = 10)
        ax.set_ylabel(y_axis, fontsize = 10)



    # only have axis titles on the outer edge
    for ax in axs.flat:
        ax.label_outer()


    # add figure title
    fig.suptitle(title)
   
    # save file
    plt.savefig(filename, dpi=800)

    return
       

def subplot_lines(axs, data, legend, lims):
    '''helper function that plots multiple lines on a subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - list of data to plot
            legend - list of labels in the same order as the lines
            lims - tuple specifying upper and lower bounds''' 
    
    plot_lines = len(data)    # number of lines to plot
  
    # iterate through the lines
    for j in range(plot_lines):
        
        # plot with or without lines 
        if legend is not None:
            axs.plot(data[j], label=legend[j])
        else:
            axs.plot(data[j])
            
    axs[row, col].set_ylim(lim[0], lim[1])

    return axs

def subplot_grid(axs, data, legend, lims):
    '''helper function for plotting values that lie on a 2-D grid.
    inputs: axs - matplotlib axes object for a single subplot
            data - 2D array of data to plot
            legend - unused paramter'''

    if len(data) > 1:
        raise ValueError("subplot_grid mishapened input")
    axs.imshow(data[0], interpolation = 'none', vmin=lims[0], vmax=lims[1])

    return axs

def subplot_scatter_2d(axs, data, legend, lims):
    '''helper function for plotting a scatter plot of data with 2
    dimensions.
    inputs: axs - matplotlib axes object for a single subplot
            data - list with each element containing 2 arrays 
                   [[x_vals0, y_vals0], [x_vals1, y_vals1]]
            legend - unused paramter'''

    plot_series = len(data)    # number of series to plot
  
    # iterate through the lines
    for j in range(plot_lines):
        
        # plot with or without lines 
        if legend is not None:
            axs.scatter(data[j][0], data[j][1], label=legend[j])
        else:
            axs.scatter(data[j][0], data[j][1])
    
    axs[row, col].set_ylim(lim[0], lim[1])

    return axs


def subplot_scatter_3d(axs, data, legend, lims):
    '''helper function for plotting a scatter plot of data with 3
    dimensions. The third dimension is illustrated with color. Only a single
    data series can be plotted per subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - list with 3 arrays [x_vals0, y_vals0, z_vals0]
            legend - unused paramter'''

    axs.scatter(data[0], data[1], c=data[2])

    return axs

def subplot_histogram(axs, data, legend, lims):
    '''helper function for plotting a histogram subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - 1D array of data
            legend - unused paramter'''

    axs.hist(data)

    return axs

def plot_histograms(filenames, x_axis, y_axis, title, fig_name):
    '''plot a histogram of some data'''
    rows = 2
    cols = 3
    # load data from pickle
    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    timesteps = [0, 200, 400, 600, 800, 999] 

    for i, filename in enumerate(filenames):
        print(filename)
        data = pickle.load(open(filename, "rb"))
        col = i % cols
        row = int(np.floor(i/cols))
        print(row, col)
        
        axs[row, col].hist(data, bins=20)
        axs[row, col].set_title("Timestep: {}".format(timesteps[i]))
        axs[row, col].set_ylim(0, 850)
        axs[row, col].grid(True)
    
    # set the axis labels
    for ax in axs.flat:
        ax.set(xlabel = x_axis, ylabel = y_axis)

    # only have axis titles on the outer edge
    for ax in axs.flat:
        ax.label_outer()
    
    plt.savefig(fig_name, dpi=800)


    return

