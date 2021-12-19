'''
Demostrate plotting all values for a 2D simulation
'''

import py_platypus as plat

if __name__ == "__main__":

    # override default simulation values
    param_dict = {
        "name": "plot_2d_test",
        "print_every": 1,
        "nppc": 1,
        "cells": [32, 32],
        "runtime": 1,
        "timestep": 0.04,
        "save_every": 1
    }

    # create parameters object
    params = plat.params.Parameters(2)
    params.set_from_dict(param_dict)

    # instantiate the 2D PIC
    pic = plat.pic_2d.PIC_2D(params)
    pic.init_x_random()
    pic.init_v_maxwellian()

    plat.run_sim.run_simulation(pic, params)
    plotter = plat.plotter.Plotter("plot_2d_test", params)
    # eventually plot the particles instead
    plotter.add_animation()
    #plotter.add_all_plots()
    #plotter.add_subplots()
    #plotter.plot_position()
    #plotter.plot_phase()
    plotter.plot_electric_field()
    plotter.plot_density()
    #plotter.plot_energy()
