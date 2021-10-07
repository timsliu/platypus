# example of setting up, running, and graphing landau damping in 3D

import sys
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter
from py_platypus.models.pic_3d import PIC_3D as PIC_3D

if __name__ == "__main__":
    params = {
        "dimensions": 3,
        "nppc": 20,
        "landau": {              # defaults for Landau damping
            "amplitude": 0.8,
            "mode": 3
        },
        "print_every": 1,
        "save_every": 1,
        "runtime": 10
    }
    
    sim_params = Parameters(3)
    sim_params.set("name", "test")
    sim_params.set_from_dict(params)

    pic = PIC_3D(sim_params)
    
    #pic.init_x_random()           # initialize random x
    #pic.init_v_maxwellian()       # maxwellian velocity distribution
    #pic.density_perturbation()    # create charge density perturbation
  
    print(pic.particle_weight)
    print(pic.dx)
