# test instantiating a 2D electrostatic PIC

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import params
from pic_2d import *
import matplotlib.pyplot as plt


if __name__ == "__main__":
    sim_params = params.Parameters(2)
    # set up parameters 
    params = {
        "length": [2 * np.pi, 2 * np.pi],
        "cells": [24, 32],
        "dimensions": 2,
        "nppc": 8
    }

    sim_params.set_from_dict(params)
    pic = PIC_2D(sim_params)
   
    # initialize x randomly and check distribution
    pic.init_x_random()
    plt.figure(1)
    plt.scatter(pic.electron_x, pic.electron_y, s=0.1)

    # initialize a maxwellian velocity distribution
    pic.init_v_maxwellian()
    plt.figure(2)
    bins = np.linspace(-3, 3, 40)
    plt.hist2d(pic.electron_vx, pic.electron_vy, bins = [bins, bins])
    plt.colorbar()
    
    # create a single stream and plot vx and vy
    pic.init_v_single_stream(1, 0.5, 2)
    plt.figure(3)
    plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vx, s = 1)
    plt.title("Vx")
    
    plt.figure(4)
    plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vy, s = 1)
    plt.title("Vy")
    
    # create a two stream setup and plot vx and vy
    pic.init_v_maxwellian()
    pic.init_v_two_beams(0.8, 0.5, 2, -2)
    plt.figure(5)
    plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vx, s = 2)
    plt.title("Two beams Vx")
    
    plt.figure(6)
    plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vy, s = 2)
    plt.title("Two beams Vy")

    # create a density perturbation
    pic.density_perturbation(0.8, 4)
    plt.figure(7)
    plt.scatter(pic.electron_x, pic.electron_y, s=1)

    pic.update_ne()
    plt.figure(8)
    ax = plt.imshow(pic.ne, interpolation = 'none')
    plt.colorbar()

    plt.show()
