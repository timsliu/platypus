# example of setting up, running, and graphing a single stream instability
# in one dimension

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
   
    # setup a single stream, 1D simulation with default parameters
    plat.run_sim.single_stream("single_stream_1d", 1)

    # load the parameters from json
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/single_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
    
    plotter = Plotter("single_stream_1d", params)
    plotter.add_subplots()   # show quantities at different steps as subplots 
    plotter.plot_all() 
