'''
Demostrate plotting all values for a 1D simulation
'''

import py_platypus as plat

if __name__ == "__main__":

    # override default simulation values
    param_dict = {
        "name": "plot_1d_test",
        "print_every": 1,
        "nppc": 32,
        "cells": [32],
        "runtime": 0.16,
        "timestep": 0.04,
        "save_every": 1
    }

    # create parameters object
    params = plat.params.Parameters(1)
    params.set_from_dict(param_dict)

    # instantiate the 1D PIC
    pic = plat.pic_1d.PIC_1D(params)
    pic.init_x_random()
    pic.init_v_two_stream()

    plat.run_sim.run_simulation(pic, params)
    plotter = plat.plotter.Plotter("plot_1d_test", params)
    # eventually plot the particles instead
    # plotter.add_animation()
    plotter.add_all_plots()
    plotter.add_subplots()
    #plotter.plot_all() 
    #plotter.plot_position()
    #plotter.plot_phase()
    #plotter.plot_density()
    #plotter.plot_velocity()
    plotter.plot_energy()
