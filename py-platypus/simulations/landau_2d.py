# example of setting up, running, and graphing landau damping

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter
from params import Parameters

if __name__ == "__main__":
    #run_sim.landau("landau_2d", 2, 
    #    param_dict={"runtime": 5, "save_every": 10}
    #)
   
    param_json = os.path.join(PLATYPUS_HOME, "py-platypus/out/landau_2d/params.json")
    params = Parameters(2, load_file=param_json) 
    plotter = Plotter("landau_2d", params)
    #plotter.plot_electric_field()
    #plotter.plot_energy()
    #plotter.plot_density()
    plotter.plot_phase() 
    #plotter.plot_velocity()
