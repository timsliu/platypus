# example of setting up, running, and graphing a single stream instability
# in one dimension

import sys
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    plat.run_sim.single_stream("single_stream_1d", 1, param_dict={"runtime": 40})
    
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/single_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
    plotter = Plotter("single_stream_1d", params)
    plotter.plot_all() 