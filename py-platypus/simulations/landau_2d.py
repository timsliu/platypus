# example of setting up, running, and graphing landau damping

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter

if __name__ == "__main__":
    #params = run_sim.landau("landau_2d", 2, param_dict={"runtime": 5})
    
    plotter = Plotter("landau_2d")
    plotter.plot_energy()
