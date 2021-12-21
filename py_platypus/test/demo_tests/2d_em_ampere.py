'''
Demonstration test for the 2D EM PIC with a uniform magnetic field illustrating
Ampere's law. The charged particles move in a stream and generate a
magnetic field
'''

import py_platypus as plat

if __name__ == "__main__":

    # override default simulation values
    param_dict = {
        "name": "ampere",
        "print_every": 1,
        "nppc": 1,
        "cells": [32, 32],
        "runtime": 50,
        "timestep": 0.05,
        "save_every": 1,
        "single_stream": {
            "stream_v": 3,
            "stream_frac": 0.8,
            "stream_width": 0.5
        }
    }

    # create parameters object
    params = plat.params.Parameters(2)
    params.set_from_dict(param_dict)

    # instantiate the 2D electromagentic PIC
    pic = plat.pic_2d_em.PIC_2D_EM(params)
    pic.init_b_uniform()
    pic.init_e()
    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.init_v_single_stream()

    plat.run_sim.run_simulation(pic, params)
    plotter = plat.plotter.Plotter("ampere", params)
    #plotter.add_animation()
    #plotter.plot_position()
    #plotter.plot_b()
    #plotter.plot_density()
