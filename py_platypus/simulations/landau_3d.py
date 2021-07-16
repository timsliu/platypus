# example of setting up, running, and graphing landau damping in 3D

import sys
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    params = {
        "dimensions": 3,
        "nppc": 10,
        "landau": {              # defaults for Landau damping
            "amplitude": 0.8,
            "mode": 3
        },
        "print_every": 1,
        "save_every": 1,
        "runtime": 1
    }
    
    plat.run_sim.landau("landau_3d", 3, 
        param_dict=params
    )
   
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/landau_3d/params.json")
    params = Parameters(3, load_file=param_json) 
    plotter = Plotter("landau_3d", params)
   
    plotter.plot_energy()
