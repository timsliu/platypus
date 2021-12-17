# plotter class for plotting different data

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pickle
import enum

import py_platypus as plat
from py_platypus.vis import vis_util as vis_util
from py_platypus.vis import animator as animator


class PlotOutputs(enum.Enum):
    """
    Enumerate the different possible ways to plot a data series
    """
    SUBPLOTS = 0  # plot multiple subplots
    ANIMATION = 1  # create an animation
    ALL_PLOTS = 2  # create individual plots for each step


class Plotter:
    """
    Plotting class with several methods for plotting physical quantities
    from the simulation.
    """
    def __init__(self, name, params):
        self.params = params

        # directory to output graphics to
        self.out_dir = os.path.join(
            plat.PLATYPUS_HOME,
            "py_platypus/out/{}/graphs".format(self.params["name"]))

        # directory with the data
        self.data_dir = os.path.join(
            plat.PLATYPUS_HOME,
            "py_platypus/out/{}/data".format(self.params["name"]))

        # output types to plot
        self.output_types = set([])

    def add_animation(self):
        """
        add animation to the list of plot types
        """
        self.add_type(PlotOutputs.ANIMATION)

    def add_subplots(self):
        """
        add subplots to the list of plot types
        """
        self.add_type(PlotOutputs.SUBPLOTS)

    def add_all_plots(self):
        """
        add all plots to the list of plot output types
        """
        self.add_type(PlotOutputs.ALL_PLOTS)

    def remove_animation(self):
        """
        remove animation to the list of plot types
        """
        self.remove_type(PlotOutputs.ANIMATION)

    def remove_subplots(self):
        """
        remove subplots to the list of plot types
        """
        self.remove_type(PlotOutputs.SUBPLOTS)

    def remove_all_plots(self):
        """
        remove all plots to the list of plot output types
        """
        self.remove_type(PlotOutputs.ALL_PLOTS)

    def add_type(self, output_type):
        """
        Add a type to the set of output types
        """
        self.output_types.add(output_type)

    def remove_type(self, output_type):
        """
        Remove a type from the set of output types
        """
        self.output_types.discard(output_type)

    def get_files(self, id_str):
        """
        get all files with a certain id string
        """
        files = os.listdir(self.data_dir)

        # filter for files with the ID type
        target_files = list(filter(lambda x: id_str in x, files))

        # order alphabetically, which also sorts by step
        target_files.sort(key=lambda x: int(x.split("_")[1]))

        # parse for the step count
        steps = [int(f.split("_")[1]) for f in target_files]

        return target_files, steps

    def plot_series(self, suffixes, out_name, x_label, y_label, title,
                    subplotter):
        """
        helper function to plot a data series across several time steps
        inputs: suffixes - list of identifier at end of pickle files to parse
                out_name - name of output file(s) with no extension
                x_label - x_label for charts
                y_label - y_label for charts
                title - title for each subplot
                subplotter - class for plotting a subplot
        """

        if len(self.output_types) == 0:
            print("No plot types listed")
            return

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

        subplotter.set_data(values)

        output_to_method = {
            PlotOutputs.SUBPLOTS: self.plot_subplots,
            PlotOutputs.ANIMATION: self.plot_animation,
            PlotOutputs.ALL_PLOTS: self.plot_all_plots
        }

        # loop through the output types and call methods to plot each output
        for output in self.output_types:
            plot_method = output_to_method[output]
            plot_method(out_name,
                        x_label,
                        y_label,
                        title,
                        subplotter,
                        steps=steps)

    def plot_all_plots(self,
                       filename,
                       data,
                       x_label,
                       y_label,
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
        lim_neg, lim_pos = vis_util.get_ylimits(data, zero, subplotter)

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
                ax.set_xlabel(x_label, fontsize=10)
                ax.set_ylabel(y_label, fontsize=10)

            # add figure title
            fig.suptitle(title)

            # save file
            plt.savefig(os.path.join(self.out_dir,
                                     "{}_{}.png".format(filename, i)),
                        dpi=800)
            plt.close("all")

        return

    def plot_subplots(self,
                      filename,
                      x_label,
                      y_label,
                      title,
                      subplotter,
                      log=False,
                      legend=None,
                      steps=None,
                      zero=False):
        """
        TODO make this a separate class like animator
        plot a single chart with multiple subplots
        inputs: filename - full path to output filename
                data - 2d array of data; dimension 0 is for each subplot, 
                       dimension 1 is for multiple lines on a subplot
                x-axis - name of x-axis
                y-axis - name of y-axis
                title - figure title
                subplotter - class for plotting a subplot
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from
        """

        subplots = len(subplotter.data)  # number of subplots
        rows, cols, = vis_util.get_subplot_config(
            subplots)  # subplot arrangement
        fig, axs = plt.subplots(rows,
                                cols,
                                constrained_layout=True,
                                squeeze=False)

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
            subplotter.plot_axes(axs[row, col], i)

            # add subplot title if available
            if steps is not None:
                axs[row, col].set_title("Step {}".format(steps[i]),
                                        fontsize=10)
            # add optional features
            if legend is not None:
                axs[row, col].legend()

            if log:
                axs[row, col].yscale("log")

            axs[row, col].grid(True)

        # set the axis labels
        for ax in axs.flat:
            ax.set_xlabel(x_label, fontsize=10)
            ax.set_ylabel(y_label, fontsize=10)

        # only have axis titles on the outer edge
        for ax in axs.flat:
            ax.label_outer()

        # add figure title
        fig.suptitle(title)

        # save file
        plt.savefig(os.path.join(self.out_dir, "{}.png".format(filename)),
                    dpi=800)
        plt.close("all")

        return

    def plot_animation(self,
                       filename,
                       x_label,
                       y_label,
                       title,
                       subplotter,
                       log=False,
                       legend=None,
                       steps=None,
                       zero=False):
        """
        Create an animation from the passed data. This method is a wrapper
        that uses the Animator class to create the animation.
        inputs: filename - full path to output filename
                data - 2d array of data; dimension 0 is for each subplot, 
                       dimension 1 is for multiple lines on a subplot
                x-axis - name of x-axis
                y-axis - name of y-axis
                title - figure title
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from
        """

        save_path = os.path.join(self.out_dir, "{}.mp4".format(filename))
        anim = animator.Animator(save_path,
                                 subplotter,
                                 x_label=x_label,
                                 y_label=y_label,
                                 title=title)

        anim.create_animation()

        return

    def plot_electric_field(self):
        '''plot the electric field at multiple timesteps'''

        print("Plotting electric field")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            self.plot_series(["ef"], "ef", "Position (x)", "Electric field",
                             "E field", plat.vis_util.subplot_lines)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            self.plot_series(["efx"], "efx", "Position (x)", "Position (y)",
                             "E field", plat.vis_util.subplot_grid)

            self.plot_series(["efy"], "efy", "Position (x)", "Position (y)",
                             "E field", plat.vis_util.subplot_grid)

        return

    def plot_energy(self):
        '''plot the kinetic and electrostatic energy as a function of time
        step'''

        print("Plotting energy")

        ee = pickle.load(open(os.path.join(self.data_dir, "ee.p"), "rb"))
        ke = pickle.load(open(os.path.join(self.data_dir, "ke.p"), "rb"))
        ee_file = os.path.join(self.out_dir, "ee")
        ke_file = os.path.join(self.out_dir, "ke")
        combined_file = os.path.join(self.out_dir, "energy")

        # plot kinetic and electrostatic energy on one chart
        self.plot_subplots(combined_file, [[ke, ee]],
                           "Time step",
                           "Normalized energy",
                           "Energy",
                           plat.vis_util.subplot_lines,
                           legend=[["Kinetic energy", "Electrostatic energy"]],
                           zero=True)

        # plot kinetic energy
        self.plot_subplots(ee_file, [[ee]],
                           "Time step",
                           "Normalized energy",
                           "Electrostatic energy",
                           plat.vis_util.subplot_lines,
                           zero=True)

        # plot potential energy
        self.plot_subplots(ke_file, [[ke]],
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
            subplotter = plat.subplotter.SubplotLines(y_axis_zero=True, x_axis_zero=True)
            self.plot_series(["ne"], "ne", "Position (x)", "Density",
                             "Electron density", subplotter)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            # TODO update for subplotter classes
            self.plot_series(["ne"], "ne", "Position (x)", "Position (y)",
                             "Electron density", plat.vis_util.subplot_grid)

        return

    def plot_phase(self):
        '''plot the phase space illustrating velocity as a function of
        position at several time steps'''

        print("Plotting phase space")

        if self.params["dimensions"] == 1:
            # TODO update for subplotter classes
            self.plot_series(["x", "v"], "phase", "Position (x)",
                             "Velocity (v)", "Phase plot",
                             plat.vis_util.subplot_scatter_2d)

        if self.params["dimensions"] == 2:
            # TODO update for subplotter classes
            self.plot_series(["ex", "ey", "evx"], "phase_vx", "Position (x)",
                             "Position (y)", "Phase plot (Vx)",
                             plat.vis_util.subplot_scatter_3d)

            self.plot_series(["ex", "ey", "evy"], "phase_vy", "Position (x)",
                             "Position (y)", "Phase plot (Vy)",
                             plat.vis_util.subplot_scatter_3d)

        return

    def plot_velocity(self):
        '''plot histogram showing velocity distribution at several different
        time steps'''

        print("Plotting velocity")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            # TODO update for subplotter classes
            self.plot_series("v", "v", "Position (x)", "Velocity", "Velocity",
                             plat.vis_util.subplot_histogram)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            # TODO update for subplotter classes
            self.plot_series(["vx", "vy"], "v", "Vx", "Vy", "Velocity",
                             plat.vis_util.subplot_histogram_2d)

    def plot_position(self):
        '''
       plot the positions of the particles
       '''

        print("Plotting particle positions")

        if self.params["dimensions"] == 2:
            # TODO update for subplotter classes
            self.plot_series(["ex", "ey"], "position", "Position (x)",
                             "Position (y)", "Position plot",
                             plat.vis_util.subplot_scatter_2d)

    def plot_all(self):
        '''plot all default graphs'''

        self.plot_electric_field()
        self.plot_energy()
        self.plot_density()
        self.plot_phase()
        self.plot_velocity()

        return
