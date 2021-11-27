"""
Test script for demonstrating creating animations using the animator class
and other py_platypus visual utilities
"""
import numpy as np
import py_platypus as plat

if __name__ == "__main__":
    frames = 200  # data frames
    points = 100  # data points per step
    data = np.random.rand(frames, 2, points)
    subplotter = plat.vis_util.subplot_scatter_2d

    animator = plat.animator.Animator("demo_animation.mp4",
                                      data,
                                      subplotter,
                                      x_label="x values",
                                      y_label="y values",
                                      title="random points")
    animator.create_animation()
