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


def plot_lines(filename, data, x_axis, y_axis, titles, 
               subplotter, log=False, legend=None):
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
    
    if legend is None:
        legend = subplots * [None]

    # iterate through subplots
    for i in range(subplots):
        # location in the subplot grid
        col = i % cols
        row = int(np.floor(i/cols))

        # call function to plot the data
        axs[row, col] = subplotter(axs[row, col], data[i], legend[i])
        
        # add axes and title
        axs[row, col].set_title(titles[i])
        axs[row, col].grid(True)
        
        # add optional features
        if legend[i] is not None:
            axs[row, col].legend()
        
        if log: 
            axs[row, col].yscale("log")

    # set the axis labels
    for ax in axs.flat:
        ax.set(xlabel = x_axis, ylabel = y_axis)

    # only have axis titles on the outer edge
    for ax in axs.flat:
        ax.label_outer()
   
    # save file
    plt.savefig(filename, dpi=800)

    return
       

def subplot_lines(axs, data, legend):
    '''helper function that plots multiple lines on a subplot
    inputs: axs - matplotlib axes object for a single subplot
            data - list of data to plot
            legend - list of labels in the same order as the lines''' 
    
    plot_lines = len(data)    # number of lines to plot
  
    # iterate through the lines
    for j in range(plot_lines):
        
        # plot with or without lines 
        if legend is not None:
            axs.plot(data[j], label=legend[j])
        else:
            axs.plot(data[j])

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

def plot_scatters_2d(x_file, y_file, x_axis, y_axis, title):
    pass

def plot_scatters_2d(x_file, y_file, x_axis, y_axis, title):
    '''plot a histogram of some data'''

    # load data from pickle
    x_data = pickle.load(open(x_file, "rb"))
    y_data = pickle.load(open(y_file, "rb"))
    
    plt.figure()
    plt.scatter(x_data, y_data, s=1)
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.grid(True)
    plt.savefig(title.replace(" ", "_") + ".png", dpi=800)

    return
