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
                    subplotter, log=False, legend=None):
        """
        helper function to plot a data series across several time steps
        inputs: suffixes - list of identifier at end of pickle files to parse
                out_name - name of output file(s) with no extension
                x_label - x_label for charts
                y_label - y_label for charts
                title - title for each subplot
                subplotter - class for plotting a subplot
                log - plot the y_label as a log plot
                legend - plot legends 
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
                        log=log,
                        legend=legend,
                        steps=steps)

    def plot_all_plots(self,
                       filename,
                       x_label,
                       y_label,
                       title,
                       subplotter,
                       log=False,
                       legend=None,
                       steps=None):
        """
        Plot a single chart for each timestep.
        inputs: filename - full path to output filename
                x_label - name of x-axis
                y_label - name of y-axis
                title - figure title
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from
        """
        plots = len(subplotter.data)
        all_data = subplotter.data  # hold data across all subplots

        for i in range(plots):
            # set data to the data from a single time step
            # but don't recalculate the limits
            subplotter.set_data([all_data[i]], reset_lims=False)
            self.plot_subplots("{}_{}".format(filename, i),
                               x_label,
                               y_label,
                               title,
                               subplotter,
                               log=log,
                               legend=legend,
                               steps=[steps[i]])
        # reset the subplotter to its original state
        subplotter.set_data(all_data)

    def plot_subplots(self,
                      filename,
                      x_label,
                      y_label,
                      title,
                      subplotter,
                      log=False,
                      legend=None,
                      steps=None):
        """
        plot a single chart with multiple subplots
        inputs: filename - full path to output filename
                x_label - name of x-axis
                y_label - name of y-axis
                title - figure title
                log - display y as log plot
                legend - list of strings labeling the data
                steps - list of time steps the data is from
        """
        subplots = len(subplotter.data)  # number of subplots
        # subplot arrangement
        rows, cols, = vis_util.get_subplot_config(subplots)
        fig, axs = plt.subplots(rows,
                                cols,
                                constrained_layout=True,
                                squeeze=False)

        # iterate through subplots
        for i in range(rows * cols):
            # location in the subplot grid
            col = i % cols
            row = int(np.floor(i / cols))
            subplot_axes = axs[row, col]

            # remove axis for unused subplots
            if i >= subplots:
                subplot_axes.set_axis_off()
                continue

            # call function to plot the data
            plot_objs = subplotter.plot_axes(subplot_axes, i)

            # add subplot title if available
            if steps is not None:
                subplot_axes.set_title("Step {}".format(steps[i]),
                                        fontsize=10)
            # add optional features
            if legend is not None:
                subplot_axes.legend(plot_objs, legend)

            if log:
                subplot_axes.yscale("log")

            subplot_axes.grid(True)
            subplot_axes.set_xlim(subplotter.get_x_limits())
            subplot_axes.set_ylim(subplotter.get_y_limits())

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
                       steps=None):
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
            subplotter = plat.subplotter.SubplotLines(x_axis_zero=True,
                                                      y_axis_zero=False)
            self.plot_series(["ef"], "ef", "Position (x)", "Electric field",
                             "E field", subplotter)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            subplotter = plat.subplotter.Subplot2DGrid()
            self.plot_series(["efx"], "efx", "Position (x)", "Position (y)",
                             "Ex field", subplotter)

            self.plot_series(["efy"], "efy", "Position (x)", "Position (y)",
                             "Ey field", subplotter)

        return

    def plot_energy(self):
        '''plot the kinetic and electrostatic energy as a function of time
        step'''

        print("Plotting energy")

        ee = pickle.load(open(os.path.join(self.data_dir, "ee.p"), "rb"))
        ke = pickle.load(open(os.path.join(self.data_dir, "ke.p"), "rb"))

        # plot kinetic and electrostatic energy on one chart
        subplotter = plat.subplotter.SubplotLines(x_axis_zero=True,
                                                  y_axis_zero=True)
        subplotter.set_data([[ke, ee, np.array(ee) + np.array(ke)]])
        self.plot_subplots(
            "total_energy",
            "Time step",
            "Normalized energy",
            "Energy",
            subplotter,
            legend=["Kinetic", "Electrostatic", "Combined"]
        )

        # plot kinetic energy
        subplotter.set_data([[ee]])
        subplotter.y_axis_zero = False
        self.plot_subplots("electrostatic_energy",
                           "Time step",
                           "Normalized energy",
                           "Electrostatic energy",
                           subplotter)

        # plot potential energy
        subplotter.set_data([[ke]])
        self.plot_subplots("kinetic_energy",
                           "Time step",
                           "Normalized energy",
                           "Kinetic energy",
                           subplotter)


    def plot_density(self):
        '''plot the electron density at multiple timesteps'''

        print("Plotting density")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            subplotter = plat.subplotter.SubplotLines(y_axis_zero=True,
                                                      x_axis_zero=True)
            self.plot_series(["ne"], "ne", "Position (x)", "Density",
                             "Electron density", subplotter)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            subplotter = plat.subplotter.Subplot2DGrid()
            self.plot_series(["ne"], "ne", "Position (x)", "Position (y)",
                             "Electron density", subplotter)

        return

    def plot_phase(self):
        '''plot the phase space illustrating velocity as a function of
        position at several time steps'''

        print("Plotting phase space")

        if self.params["dimensions"] == 1:
            subplotter = plat.subplotter.SubplotScatter2D(y_axis_zero=False,
                                                          x_axis_zero=True)
            self.plot_series(["x", "v"], "phase", "Position (x)",
                             "Velocity (v)", "Phase plot", subplotter)

        if self.params["dimensions"] == 2:
            subplotter = plat.subplotter.SubplotScatter3D()
            self.plot_series(["ex", "ey", "evx"], "phase_vx", "Position (x)",
                             "Position (y)", "Phase plot (Vx)", subplotter)

            self.plot_series(["ex", "ey", "evy"], "phase_vy", "Position (x)",
                             "Position (y)", "Phase plot (Vy)", subplotter)

        return

    def plot_velocity(self):
        '''plot histogram showing velocity distribution at several different
        time steps'''

        print("Plotting velocity")

        # call helper for single dimension
        if self.params["dimensions"] == 1:
            subplotter = plat.subplotter.SubplotHistogram()
            self.plot_series("v", "v", "Velocity", "Count", "Velocity",
                             subplotter)

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            subplotter = plat.subplotter.SubplotHistogram2D()
            self.plot_series(["vx", "vy"], "v", "Vx", "Vy", "Velocity",
                             subplotter)

    def plot_position(self):
        '''
       plot the positions of the particles
       '''

        print("Plotting particle positions")

        if self.params["dimensions"] == 2:
            subplotter = plat.subplotter.SubplotScatter2D(y_axis_zero=True,
                                                          x_axis_zero=True)
            self.plot_series(["ex", "ey"], "position", "Position (x)",
                             "Position (y)", "Position plot", subplotter)

    def plot_all(self):
        '''plot all default graphs'''

        self.plot_electric_field()
        self.plot_energy()
        self.plot_density()
        self.plot_phase()
        self.plot_velocity()

        return
