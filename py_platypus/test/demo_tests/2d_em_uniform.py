'''
Demonstration test for the 2D EM PIC with a uniform magnetic field. The
particles in the simulation should exhibit cyclotron motion.
'''

import py_platypus as plat

if __name__ == "__main__":

    # override default simulation values
    param_dict = {
        "name": "cyclotron",
        "print_every": 1,
        "nppc": 0.01,
        "cells": [32, 32],
        "runtime": 50,
        "timestep": 0.05,
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
    pic.init_v_maxwellian()

    plat.run_sim.run_simulation(pic, params)
    plotter = plat.plotter.Plotter("cyclotron", params)
    plotter.add_animation()
    plotter.plot_position()
    #plotter.plot_density()
