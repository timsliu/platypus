# example of setting up, running, and graphing a two stream instability in
# two dimensions

import sys
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":

    params = {
        "print_every": 5,
        "nppc": 30,
        "cells": [32, 32],
        "runtime": 40,
        "save_every": 100,
        "two_stream": {
            "vpos": 1,
            "vneg": -1,
            "stream_frac": 1,
            "stream_width": 1,
        }
    }
    
    plat.run_sim.two_stream("two_stream_2d", 2, param_dict = params)
    
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/two_stream_2d/params.json")
    params = Parameters(2, load_file=param_json) 
    plotter = Plotter("two_stream_2d", params)
    plotter.plot_all()


