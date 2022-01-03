# example of setting up, running, and graphing landau damping

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    
    # set up and run Landau damping simulation in 2D
    plat.run_sim.landau("landau_2d", 2, 
        param_dict={"runtime": 5, "save_every": 1}
    )
  
    # load parameters from a json file
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/landau_2d/params.json")
    params = Parameters(2, load_file=param_json) 
   
    # plot all quantities and generate animations
    plotter = Plotter("landau_2d", params)
    plotter.add_animation() 
    plotter.plot_all()
