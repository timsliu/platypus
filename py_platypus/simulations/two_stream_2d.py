# example of setting up, running, and graphing a two stream instability in
# two dimensions

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":

    # dictionary of override parameters
    params = {
        "print_every": 5,
        "nppc": 10,
        "cells": [32, 32],
        "runtime": 20,
        "save_every": 5,
        "two_stream": {
            "vpos": 1,
            "vneg": -1,
            "stream_frac": 1,
            "stream_width": 0.5,
        }
    }
    
    # set up a two stream simulation in 2 dimensions
    plat.run_sim.two_stream("two_stream_2d", 2, param_dict=params)
   
    # load the parameters from a json file
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/two_stream_2d/params.json")
    params = Parameters(2, load_file=param_json) 
   
    # create instance of the plotting class
    plotter = Plotter("two_stream_2d", params)
    plotter.add_animation()
    plotter.plot_all()  # plot all quantities


