# simple test script that exercises methods in PIC_2D not tested
# by 2d_test.py

# 2d_test.py

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
        "length": [2 * np.pi, 4 * np.pi],
        "cells": [32, 64],
        "dimensions": 2,
        "nppc": 8
    }

    sim_params.set_from_dict(params)
    pic = PIC_2D(sim_params)

    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.batch = [0, 1, 2, 3, 4]

    for i in range(10):
        print("Simple 2d test step: {}".format(i))
        pic.step()
        pic.calc_electrostatic_energy()
        pic.calc_kinetic_energy()
        pic.calc_batch_kinetic_energy()

    plt.figure(1)
    plt.plot(pic.output["kinetic_energy"])
    plt.title("Kinetic energy")
    plt.xlabel("Time step")
    plt.ylabel("KE")
    
    plt.figure(2)
    plt.plot(pic.output["electrostatic_energy"])
    plt.title("Electrostatic energy")
    plt.xlabel("Time step")
    plt.ylabel("EE")

    plt.figure(3)
    plt.plot(pic.output["batch_ke"])
    plt.title("Batch kinetic energy")
    plt.xlabel("Time step")
    plt.ylabel("Batch KE")

    plt.show()
