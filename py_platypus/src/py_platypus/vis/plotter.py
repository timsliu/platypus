# plotter class for plotting different data

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pickle

import py_platypus as plat
from py_platypus.vis import vis_util as vis_util


class Plotter:
    def __init__(self, name, params):
        self.params = params

        self.out_dir = os.path.join(
            plat.PLATYPUS_HOME,
            "py_platypus/out/{}/graphs".format(self.params["name"]))

        self.data_dir = os.path.join(
            plat.PLATYPUS_HOME,
            "py_platypus/out/{}/data".format(self.params["name"]))

        self.subplots = []

    def get_files(self, id_str):
        '''get all files with a certain id string'''
        files = os.listdir(self.data_dir)

        # filter for files with the ID type
        target_files = list(filter(lambda x: id_str in x, files))

        # order alphabetically, which also sorts by step
        target_files.sort(key=lambda x: int(x.split("_")[1]))

        # parse for the step count
        steps = [int(f.split("_")[1]) for f in target_files]

        return target_files, steps

    def plot_series(self,
                    suffixes,
                    out_name,
                    x_axis,
                    y_axis,
                    title,
                    subplotter,
                    animate=False,
                    plot_all=False):
        '''helper function to plot a data series across several time steps
        inputs: suffixes - identifier at end of pickle files to parse
                out_name - name of output file with .png extension
                x_axis - x_axis for charts
                y_axis - y_axis for charts
                title - title for each subplot
                subplotter - function for plotting a subplot'''

        self.subplots = []  # list of axs objects for each subplot

        values = []  # array holding data for each subplot
        all_files = []  # 2D array with suffixes by subplots elements

        # iterate through the data types we need to open
        for suffix in suffixes:
            files, steps = self.get_files(suffix)
            all_files.append(files)

        # iterate through the data for each subplot
        for i in range(len(steps)):
            subplot_values = []
            # add the data for each data type at each timestep
            for j in range(len(suffixes)):
                subplot_values.append(
                    pickle.load(
                        open(os.path.join(self.data_dir, all_files[j][i]),
                             "rb")))
            values.append(subplot_values)

        # name of output file name
        graph_file_name = os.path.join(self.out_dir, out_name)

        # plot the subplots
        self.plot_subplots(graph_file_name,
                           values,
                           x_axis,
                           y_axis,
                           title,
                           subplotter,
                           steps=steps)

        # optionally create an animation
        if animate:
            self.plot_animation(graph_file_name,
                                values,
                                x_axis,
                                y_axis,
                                title,
                                subplotter,
                                steps=steps)

        # optionally plot separate plots
        if plot_all:
            self.plot_all_plots(graph_file_name,
                                values,
                                x_axis,
                                y_axis,
                                title,
                                subplotter,
                                steps=steps)

    def plot_all_plots(self,
                       filename,
                       data,
                       x_axis,
                       y_axis,
                       title,
                       subplotter,
                       log=False,
                       legend=None,
                       steps=None,
                       zero=False):
        '''plot a single chart for each timestep
        inputs: filename - full path to output filename
                data - 2d array of data; dimension 0 is for each plot, 
                       dimension 1 is for multiple data series on a plot
                x-axis - name of x-axis
                y-axis - name of y-axis
                title - figure title
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from'''

        plots = len(data)
        # call function to get upper and lower figure limits
        lim_neg, lim_pos = vis_util.get_limits(data, zero, subplotter)

        # create filename_{}.png from filename.png so that step can be added
        # to the name of each tile
        dot_index = filename.find(".")
        filename_format = filename[:dot_index] + "_{}" + filename[dot_index:]

        if legend is None:
            legend = plots * [None]

        for i in range(plots):
            fig, axs = plt.subplots(1,
                                    1,
                                    constrained_layout=True,
                                    squeeze=False)
            subplotter(axs[0, 0], data[i], legend[i], (lim_neg, lim_pos))

            # add optional features
            if legend is not None:
                axs[0, 0].legend()

            if log:
                axs[0, 0].yscale("log")

            axs[0, 0].grid(True)

            # set the axis labels
            for ax in axs.flat:
                ax.set_xlabel(x_axis, fontsize=10)
                ax.set_ylabel(y_axis, fontsize=10)

            # add figure title
            fig.suptitle(title)

            # save file
            plt.savefig(filename_format.format(i), dpi=800)

        return

    def plot_subplots(self,
                      filename,
                      data,
                      x_axis,
                      y_axis,
                      title,
                      subplotter,
                      log=False,
                      legend=None,
                      steps=None,
                      zero=False):
        '''plot a single chart with multiple subplots
        inputs: filename - full path to output filename
                data - 2d array of data; dimension 0 is for each subplot, 
                       dimension 1 is for multiple lines on a subplot
                x-axis - name of x-axis
                y-axis - name of y-axis
                title - figure title
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from'''

        subplots = len(data)  # number of subplots
        rows, cols, = vis_util.get_subplot_config(
            subplots)  # subplot arrangement
        fig, axs = plt.subplots(rows,
                                cols,
                                constrained_layout=True,
                                squeeze=False)

        # call function to get upper and lower figure limits
        lim_neg, lim_pos = vis_util.get_limits(data, zero, subplotter)

        if legend is None:
            legend = subplots * [None]

        # iterate through subplots
        for i in range(rows * cols):
            # location in the subplot grid
            col = i % cols
            row = int(np.floor(i / cols))

            # remove axis for unused subplots
            if i >= subplots:
                axs[row, col].set_axis_off()
                continue

            # call function to plot the data
            subplotter(axs[row, col], data[i], legend[i], (lim_neg, lim_pos))

            # add subplot title if available
            if steps is not None:
                axs[row, col].set_title("Step {}".format(steps[i]),
                                        fontsize=10)

            # add optional features
            if legend[i] is not None:
                axs[row, col].legend()

            if log:
                axs[row, col].yscale("log")

            axs[row, col].grid(True)

        # set the axis labels
        for ax in axs.flat:
            ax.set_xlabel(x_axis, fontsize=10)
            ax.set_ylabel(y_axis, fontsize=10)

        # only have axis titles on the outer edge
        for ax in axs.flat:
            ax.label_outer()

        # add figure title
        fig.suptitle(title)

        # save file
        plt.savefig(filename, dpi=800)

        return

    def plot_electric_field(self):
        '''plot the electric field at multiple timesteps'''

        print("Plotting electric field")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            self.plot_series(["ef"], "ef.png", "Position (x)",
                             "Electric field", "E field",
                             plat.vis_util.subplot_lines)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            self.plot_series(["efx"], "efx.png", "Position (x)",
                             "Position (y)", "E field",
                             plat.vis_util.subplot_grid)

            self.plot_series(["efy"], "efy.png", "Position (x)",
                             "Position (y)", "E field",
                             plat.vis_util.subplot_grid)

        return

    def plot_energy(self):
        '''plot the kinetic and electrostatic energy as a function of time
        step'''

        print("Plotting energy")

        ee = pickle.load(open(os.path.join(self.data_dir, "ee.p"), "rb"))
        ke = pickle.load(open(os.path.join(self.data_dir, "ke.p"), "rb"))
        ee_file = os.path.join(self.out_dir, "ee.png")
        ke_file = os.path.join(self.out_dir, "ke.png")
        combined_file = os.path.join(self.out_dir, "energy.png")

        # plot kinetic and electrostatic energy on one chart
        plat.vis_util.plot_subplots(
            combined_file, [[ke, ee]],
            "Time step",
            "Normalized energy",
            "Energy",
            plat.vis_util.subplot_lines,
            legend=[["Kinetic energy", "Electrostatic energy"]],
            zero=True)

        # plot kinetic energy
        plat.vis_util.plot_subplots(ee_file, [[ee]],
                                    "Time step",
                                    "Normalized energy",
                                    "Electrostatic energy",
                                    plat.vis_util.subplot_lines,
                                    zero=True)

        # plot potential energy
        plat.vis_util.plot_subplots(ke_file, [[ke]],
                                    "Time step",
                                    "Normalized energy",
                                    "Kinetic energy",
                                    plat.vis_util.subplot_lines,
                                    zero=True)

        return

    def plot_density(self):
        '''plot the electron density at multiple timesteps'''

        print("Plotting density")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            self.plot_series(["ne"], "ne.png", "Position (x)", "Density",
                             "Electron density", plat.vis_util.subplot_lines)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            self.plot_series(["ne"], "ne.png", "Position (x)", "Position (y)",
                             "Electron density", plat.vis_util.subplot_grid)

        return

    def plot_phase(self):
        '''plot the phase space illustrating velocity as a function of
        position at several time steps'''

        print("Plotting phase space")

        if self.params["dimensions"] == 1:
            self.plot_series(["x", "v"], "phase.png", "Position (x)",
                             "Velocity (v)", "Phase plot",
                             plat.vis_util.subplot_scatter_2d)

        if self.params["dimensions"] == 2:
            self.plot_series(["ex", "ey", "evx"], "phase_vx.png",
                             "Position (x)", "Position (y)", "Phase plot (Vx)",
                             plat.vis_util.subplot_scatter_3d)

            self.plot_series(["ex", "ey", "evy"], "phase_vy.png",
                             "Position (x)", "Position (y)", "Phase plot (Vy)",
                             plat.vis_util.subplot_scatter_3d)

        return

    def plot_velocity(self):
        '''plot histogram showing velocity distribution at several different
        time steps'''

        print("Plotting velocity")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            self.plot_series("v", "v.png", "Position (x)", "Velocity",
                             "Velocity", plat.vis_util.subplot_histogram)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            self.plot_series(["vx", "vy"], "v.png", "Vx", "Vy", "Velocity",
                             plat.vis_util.subplot_histogram_2d)

    def plot_position(self):
        '''
       plot the positions of the particles
       '''

        print("Plotting particle positions")

        if self.params["dimensions"] == 2:
            self.plot_series(["ex", "ey"],
                             "position.png",
                             "Position (x)",
                             "Position (y)",
                             "Position plot",
                             plat.vis_util.subplot_scatter_2d,
                             plot_all=True)

    def plot_all(self):
        '''plot all default graphs'''

        self.plot_electric_field()
        self.plot_energy()
        self.plot_density()
        self.plot_phase()
        self.plot_velocity()

        return
