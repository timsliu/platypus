# example of setting up, running, and graphing landau damping

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter

if __name__ == "__main__":
    #run_sim.landau("landau_1d", 1, param_dict={"runtime": 10, "save_every": 20})
    
    plotter = Plotter("landau_1d")
    plotter.plot_electric_field()
