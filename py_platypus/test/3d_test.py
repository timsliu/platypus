# test instantiating a 3D electrostatic PIC

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import py_platypus as plat


if __name__ == "__main__":
    sim_params = plat.params.Parameters(3)
    # set up parameters 
    params = {
        "length": [2 * np.pi, 4 * np.pi, 6 * np.pi],
        "cells": [32, 32, 32],
        "dimensions": 3,
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

    sim_params.set_from_dict(params)
    pic = plat.pic_3d.PIC_3D(sim_params)
   
    # initialize x randomly and check distribution
    pic.init_x_random()
    fig = plt.figure(1)
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pic.electron_x, pic.electron_y, pic.electron_z, s=0.05)

    # initialize a maxwellian velocity distribution
    pic.init_v_maxwellian()
    fig, axs = plt.subplots(1, 3, constrained_layout=True)
    axs[0].hist(pic.electron_vx, bins=20)
    axs[1].hist(pic.electron_vy, bins=20)
    axs[2].hist(pic.electron_vz, bins=20)
    
    # create a single stream and plot vx, vy, and vz
    # set up parameters 
    pic.init_v_single_stream()

    plt.figure(3)
    plt.scatter(pic.electron_y, pic.electron_z, c=pic.electron_vx, s = 1)
    plt.title("Vx")
    plt.xlabel("y")
    plt.ylabel("z")

    plt.figure(4)
    plt.scatter(pic.electron_y, pic.electron_z, c=pic.electron_vy, s = 1)
    plt.title("Vy")
    plt.xlabel("y")
    plt.ylabel("z")
    
    
    # create a two stream setup and plot vx and vy
    pic.init_v_maxwellian()
    pic.init_v_two_stream()
    
    plt.figure(5)
    plt.scatter(pic.electron_y, pic.electron_z, c=pic.electron_vx, s = 1)
    plt.title("Two stream Vx")
    plt.xlabel("y")
    plt.ylabel("z")

    plt.figure(6)
    plt.scatter(pic.electron_y, pic.electron_z, c=pic.electron_vy, s = 1)
    plt.title("Two stream Vy")
    plt.xlabel("y")
    plt.ylabel("z")

    # create a density perturbation
    pic.density_perturbation()
    plt.figure(7)
    plt.scatter(pic.electron_x, pic.electron_y, s=1)
    plt.title("Two stream xy projection")
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.figure(8)
    plt.scatter(pic.electron_x, pic.electron_y, s=1)
    plt.title("Two stream xz projection")
    plt.xlabel("x")
    plt.ylabel("z")
    
    fig = plt.figure(9)
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pic.electron_x, pic.electron_y, pic.electron_z, s=0.05)

    pic.update_ne()
    print("Average number density: ", np.mean(pic.ne))
    xy_number_density = np.sum(pic.ne, axis=2)   # collapse the z dimension 
    xz_number_density = np.sum(pic.ne, axis=1)   # collapse the y dimension 

    plt.figure(10)
    plt.title("Electron number density xy")
    plt.imshow(xy_number_density, interpolation = 'none')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar()
    
    plt.figure(11)
    plt.title("Electron number density xz")
    plt.imshow(xz_number_density, interpolation = 'none')
    plt.xlabel("x")
    plt.ylabel("z")
    plt.colorbar()

    ## calculate charge density
    #pic.update_ne()
    #pic.update_ni()
    #pic.update_rho()
    #plt.figure(9)
    #ax = plt.imshow(pic.rho, interpolation = 'none')
    #plt.title("Charge density")
    #plt.colorbar()

    # test calculating phi from rho
    #sin_2d = np.zeros(pic.cells)
    #for i in range(pic.cells[0]):
    #    for j in range(pic.cells[1]):
    #        sin_2d[i][j] = np.sin(pic.dx[0] * i ) + \
    #                       np.sin(pic.dx[1] * j )

    #pic.rho = sin_2d
    #pic.update_phi()

    #plt.figure(10)
    #ax = plt.imshow(pic.rho, interpolation = 'none')
    #plt.title("Sin Charge density")
    #plt.colorbar()
    #
    #plt.figure(11)
    #ax = plt.imshow(pic.phi, interpolation = 'none')
    #plt.title("Electric potential")
    #plt.colorbar()

    ## test calculating electric field at the nodes
    #pic.update_e()
    #plt.figure(12)
    #ax = plt.imshow(pic.ex, interpolation = 'none')
    #plt.title("Electric field Ex")
    #plt.colorbar()
    #
    #plt.figure(13)
    #ax = plt.imshow(pic.ey, interpolation = 'none')
    #plt.title("Electric field Ey")
    #plt.colorbar()

    ## test updating particle velocity
    #pic.electron_vx = np.zeros(pic.n_particles)   # zero out particle velocity
    #pic.electron_vy = np.zeros(pic.n_particles)

    #pic.update_v()        # update velocity based on E field

    #plt.figure(14)
    #plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vx, s = 2)
    #plt.title("Velocity vx")
    #plt.colorbar()
    #
    #plt.figure(15)
    #plt.scatter(pic.electron_x, pic.electron_y, c=pic.electron_vy, s = 2)
    #plt.title("Velocity vy")
    #plt.colorbar()

    plt.show()
