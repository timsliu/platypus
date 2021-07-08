# example of setting up, running, and graphing a single stream instability
# in two dimensions

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter
from params import Parameters 

if __name__ == "__main__":
    
    params = {
        "single_stream": {
            "stream_v": 1.5, 
            "stream_width": 1.5, 
            "stream_frac": 0.9
        },
        "runtime": 10,
        "print_every": 5,
        "nppc": 25
    }
    
    run_sim.single_stream("single_stream_2d", 2, param_dict=params)
    
    param_json = os.path.join(PLATYPUS_HOME, "py-platypus/out/single_stream_2d/params.json")
    params = Parameters(2, load_file=param_json)
    plotter = Plotter("single_stream_2d", params)
    plotter.plot_all() 
