"""
Animator class for generating matplotlib based animations.
"""

import matplotlib.pyplot as plt
import matplotlib.animation as ani

from py_platypus.vis import vis_util as vis_util


class Animator:
    def __init__(self,
                 save_path,
                 subplotter,
                 x_label=None,
                 y_label=None,
                 title=None,
                 interval=40,
                 blit=True,
                 dpi=800,
                 verbose=True):

        self.save_path = save_path  # path to save figure to
        self.x_label = x_label  # label for x_axis
        self.y_label = y_label  # label for y_axis
        self.title = title  # plot title
        self.subplotter = subplotter  # subplotting class
        self.interval = interval  # time between frames in miliseconds
        self.blit = blit  # only animate parts of frame that change
        self.verbose = verbose

        self.fig = plt.figure(dpi=dpi)
        self.ax = plt.axes()
        self.xlim = subplotter.get_xlimits()
        self.ylim = subplotter.get_ylimits()
        self.frames = len(subplotter.data)

    def setup_axes(self):
        """
        Add legend, labels, and resize the axes
        """
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        if self.x_label is not None:
            self.ax.set_xlabel(self.x_label)
        if self.y_label is not None:
            self.ax.set_ylabel(self.y_label)
        if self.title is not None:
            self.ax.set_title(self.title)

    def animate_frame(self, i):
        """
        Return a single frame
        """
        self.ax.clear()  # clear the data from the old frame
        self.setup_axes()
        plot_obj = self.subplotter.plot_axes(self.ax, i)

        if self.verbose and i % 50 == 0:
            print("Animating figure frame {} of {}".format(i, self.frames))
        return plot_obj,

    def create_animation(self):
        """
        Create and save the animation.
        """
        print("Starting animation...")
        anim = ani.FuncAnimation(self.fig,
                                 self.animate_frame,
                                 frames=self.frames,
                                 interval=self.interval,
                                 blit=self.blit)
        anim.save(self.save_path, writer='ffmpeg')
        print("Finished generating animation, saved to {}".format(
            self.save_path))
