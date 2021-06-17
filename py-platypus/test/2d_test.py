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
    params = {
        "length": [2 * np.pi, np.pi],
        "cells": [32, 16],
        "dimensions": 2,
        "nppc": 8
    }

    sim_params.set_from_dict(params)
    pic = PIC_2D(sim_params)
    pic.init_x_random()
    pic.init_v_maxwellian()

    plt.figure(1)
    plt.scatter(pic.electron_x, pic.electron_y, s=0.1)

    plt.figure(2)
    bins = np.linspace(-3, 3, 40)
    plt.hist2d(pic.electron_vx, pic.electron_vy, bins = [bins, bins])
      
    # Adding color bar
    plt.colorbar()
      
    plt.show()
