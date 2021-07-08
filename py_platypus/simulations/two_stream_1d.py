# example of setting up, running, and graphing a two stream instability in
# one dimension

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter
from params import Parameters 


if __name__ == "__main__":
    run_sim.two_stream("two_stream_1d", 1, param_dict={"runtime": 40})
    
    param_json = os.path.join(PLATYPUS_HOME, "py-platypus/out/two_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
    plotter = Plotter("two_stream_1d", params)
    plotter.plot_all() 
