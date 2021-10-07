"""
Common functions and helper classes for unit tests
"""

import matplotlib.pyplot as plt
import numpy as np
import py_platypus as plat

def setup(pic_type):
    """
    Function for setting up different types of PIC instances
    """


    pic_map = {"2d_em": setup_2d_em_pic}
    
    if pic_type not in pic_map.keys():
        raise ValueError("Pic type {} not recognized by setup function".format(pic_type))

    return pic_map[pic_type]()


def setup_2d_em_pic():
    """
    Returns a 2D electromagnetic PIC for testing
    """
    
    params = {
        "length": [2 * np.pi, 2 * np.pi],
        "cells": [32, 32],
        "dimensions": 2,
        "nppc": 10,
        "single_stream": {       # defaults for single stream instability
            "stream_v": 3, 
            "stream_frac": 0.8, 
            "stream_width": 1
        },
        "landau": {              # defaults for Landau damping
            "amplitude": 0.8,
            "mode": 3
        },
        "two_stream": {          # defaults for two stream instability
            "vpos": 2,
            "vneg": -2,
            "stream_frac": 1,
            "stream_width": 0.8
        },
    }

    sim_params = plat.params.Parameters(2)
    sim_params.set_from_dict(params)
    pic = plat.pic_2d_em.PIC_2D_EM(sim_params)

    return pic


def check_graph(expect_str=""):
    """
    Helper function that displays a plot, asks for user feedback if the plot
    looks correct, then returns a boolean representing if the plot matches
    what's expected
    inputs - expect_str: optional string describing the expected image
    """

    plt.show(block=False)
   
    plot_good = "x"
    while plot_good not in ["y", "n"]:
        plot_good = input("Plot matches expectation: {} (y/n)? ".format(expect_str))

    plt.close("all")
    return plot_good == "y"
