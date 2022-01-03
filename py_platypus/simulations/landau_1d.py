# example of setting up, running, and graphing landau damping
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":

    # set up and run a 1D simulation illustrating Landau damping
    plat.run_sim.landau("landau_1d", 1, param_dict={"runtime": 10, "save_every": 1})
 
    # load the parameters from the json file
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/landau_1d/params.json")
    params = Parameters(1, load_file=param_json)
    
    # create instance of the plotting class
    plotter = Plotter("landau_1d", params)
    plotter.add_animation()  
    plotter.plot_all()
