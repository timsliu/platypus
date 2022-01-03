# example of setting up, running, and graphing a single stream instability
# in two dimensions

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":

    # override default parameters
    params = {
        "single_stream": {
            "stream_v": 1.5,     # velocity of the stream
            "stream_width": 1.5, # how wide the stream will be  
            "stream_frac": 0.9  # fraction of particles in the band in the stream
        },
        "runtime": 10,
        "print_every": 5,
        "nppc": 25
    }
    
    # set up and run a single stream simulation in 2 dimensions
    plat.run_sim.single_stream("single_stream_2d", 2, param_dict=params)
    
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/single_stream_2d/params.json")
    params = Parameters(2, load_file=param_json)
    
    # plot all quantities at different timesteps on the same subplot
    plotter = Plotter("single_stream_2d", params)
    plotter.add_subplots() 
    plotter.plot_all() 
