# example of setting up, running, and graphing landau damping in 3D

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    
    # override default parameters
    params = {
        "dimensions": 3,
        "nppc": 4,
        "landau": {              # defaults for Landau damping
            "amplitude": 0.8,
            "mode": 3   # number of density peaks
        },
        "print_every": 1,   # print current step at every step
        "save_every": 1,  # save data at every time step
        "runtime": 1
    }
   
    # set up and run Landau damping simulation in 2D
    plat.run_sim.landau("landau_3d", 3, 
        param_dict=params
    )
   
    # load parameters from a json file
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/landau_3d/params.json")
    params = Parameters(3, load_file=param_json) 
   
    # create instance of plotting class and plot the energy
    # for the three dimension case, only the energy can be plotted
    plotter = Plotter("landau_3d", params)
    plotter.plot_energy()
