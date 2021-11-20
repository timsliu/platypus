'''
Demonstration test for the 2D EM PIC with a uniform magnetic field. The
particles in the simulation should exhibit cyclotron motion.
'''

import py_platypus as plat

if __name__ == "__main__":

    # override default simulation values
    param_dict = {
        "name": "cyclotron",
        "print_every": 5,
        "nppc": 1,
        "cells": [32, 32],
        "runtime": 1,
        "save_every": 1
    }

    # create parameters object
    params = plat.params.Parameters(2)
    params.set_from_dict(param_dict)

    # instantiate the 2D electromagentic PIC
    pic = plat.pic_2d_em.PIC_2D_EM(params)
    pic.init_b_uniform(value=1.0)
    pic.init_e()
    pic.init_x_random()

    plat.run_sim.run_simulation(pic, params)
    plotter = Plotter("cyclotron")
    # eventually plot the particles instead
    plotter.plot_density()
