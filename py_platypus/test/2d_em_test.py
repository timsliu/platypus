import py_platypus as plat
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # set up parameters 
    params = {
        "length": [2 * np.pi, 2 * np.pi, 2 * np.pi],
        "cells": [32, 32, 32],
        "dimensions": 2,
        "nppc": 10,
        "single_stream": {       # defaults for single stream instability
            "stream_v": 3, 
            "stream_frac": 0.8, 
            "stream_width": 1
        },
        "landau": {              # defaults for Landau damping
            "amplitude": 0.8,
            "mode": 3
        },
        "two_stream": {          # defaults for two stream instability
            "vpos": 2,
            "vneg": -2,
            "stream_frac": 1,
            "stream_width": 0.8
        },
    }

    sim_params = plat.params.Parameters(2)
    sim_params.set_from_dict(params)
    pic = plat.pic_2d_em.PIC_2D_EM(sim_params)



