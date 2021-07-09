# example of setting up, running, and graphing landau damping

import sys
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    plat.run_sim.landau("landau_2d", 2, 
        param_dict={"runtime": 10, "save_every": 10}
    )
   
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/landau_2d/params.json")
    params = Parameters(2, load_file=param_json) 
    plotter = Plotter("landau_2d", params)
   
    plotter.plot_all()
