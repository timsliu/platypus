# test a 3D electrostatic PIC

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import py_platypus as plat

def test_random_x(pic):
    '''test the random distribution of position x'''
    # initialize x randomly and check distribution
    pic.init_x_random()
    fig = plt.figure(1)
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pic.electron_x, pic.electron_y, pic.electron_z, s=0.05)
    
    return

def test_maxwellian(pic):
    '''test the maxwellian velocity distribution'''
    # initialize a maxwellian velocity distribution
    pic.init_v_maxwellian()
    fig, axs = plt.subplots(1, 3, constrained_layout=True)
    axs[0].hist(pic.electron_vx, bins=20)
    axs[1].hist(pic.electron_vy, bins=20)
    axs[2].hist(pic.electron_vz, bins=20)
    
    return

def test_single_stream(pic):
    '''test creating a single stream velocity distribution'''
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
    
def test_two_stream(pic):
    '''test creating a two stream velocity distribution'''
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

def test_density_perturbation(pic):
    '''test creating a density perturbation'''
    # create a density perturbation
    pic.init_x_random()
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

    return

def test_update_n(pic):
    '''test updating the number density'''
    # test updating number density
    pic.init_x_random()
    pic.density_perturbation()
    pic.update_ne()
    
    print("Average number density: ", np.mean(pic.ne))
    xy_number_density = np.sum(pic.ne, axis=2)   # collapse the z dimension 
    xz_number_density = np.sum(pic.ne, axis=0)   # collapse the y dimension 

    # imshow treats input as (rows, cols) but density is in (x, y, z) coords
    # so swap the labels
    plt.figure(10)
    plt.title("Electron number density xy")
    plt.imshow(xy_number_density, interpolation = 'none')
    plt.xlabel("x") 
    plt.ylabel("y")
    plt.colorbar()
    
    plt.figure(11)
    plt.title("Electron number density xz")
    plt.imshow(xz_number_density, interpolation = 'none')
    plt.xlabel("z")
    plt.ylabel("x")
    plt.colorbar()

def test_rho(pic):
    '''test updating the charge density rho'''
    # calculate charge density
    pic.init_x_random()
    pic.density_perturbation()
    pic.update_ne()
    pic.update_rho()
    print("Average charge density: ", np.mean(pic.rho))
    xy_charge_density = np.sum(pic.rho, axis=2)   # collapse the z dimension 
    xz_charge_density = np.sum(pic.rho, axis=0)   # collapse the y dimension 

    plt.figure(12)
    plt.title("Charge density xy")
    plt.imshow(xy_charge_density, interpolation = 'none')
    plt.xlabel("x") 
    plt.ylabel("y")
    plt.colorbar()
    
    plt.figure(13)
    plt.title("Charge density xz")
    plt.imshow(xz_charge_density, interpolation = 'none')
    plt.xlabel("z")
    plt.ylabel("x")
    plt.colorbar()

def test_phi(pic):
    '''test updating the potential from the charge density'''
    # test calculating phi from rho
    sin_3d = np.zeros(pic.cells)
   
    # create a phi field sin(x) + sin(y) + sin(z)
    for i in range(pic.cells[0]):
        for j in range(pic.cells[1]):
            for k in range(pic.cells[2]):
                sin_3d[i][j][k] = np.sin(pic.dx[0] * i ) + \
                                  np.sin(pic.dx[1] * j ) + \
                                  np.sin(pic.dx[2] * k )

    pic.rho = sin_3d
    pic.update_phi()

    # compare against sin_3d; negative second derivative of sin is sin
    abs_error = np.abs((pic.phi - sin_3d))
    print("Average poisson solver error: ", np.mean(abs_error))
    print("Max poisson solver error: ", np.max(abs_error))
    print("Min poisson solver error: ", np.min(abs_error))

    return

def test_e_flat(pic):
    '''test calculating the electric field from the potential'''
    
    # create field with a uniform gradient
    flat_grad = np.zeros(pic.cells)
    for i in range(pic.cells[0]):
        for j in range(pic.cells[1]):
            for k in range(pic.cells[2]):
                flat_grad[i][j][k] = i + j + k
    
    pic.phi = flat_grad 
    pic.update_e()

    print("1 over dx: ", 1/np.array(pic.dx))
    print("Average E: ", np.mean(pic.ex), np.mean(pic.ey), np.mean(pic.ez))
    print("Max E:     ", np.max(pic.ex), np.max(pic.ey), np.max(pic.ez))
    print("Min E:     ", np.min(pic.ex), np.min(pic.ey), np.min(pic.ez))
 
    # plot the electric field in each direction; use random y to distribute points
    # for easier visualization

    plt.figure(14)
    plt.scatter(pic.ex.flatten(), np.random.rand(pic.ex.size), s=0.1)
    plt.title("Ex") 
    
    plt.figure(15)
    plt.scatter(pic.ex.flatten(), np.random.rand(pic.ey.size), s=0.1)
    plt.title("Ey") 
    
    plt.figure(16)
    plt.scatter(pic.ex.flatten(), np.random.rand(pic.ez.size), s=0.1)
    plt.title("Ez")

    plt.imshow(pic.ez[0], interpolation = 'none')
    #plt.imshow(pic.ex[1], interpolation = 'none')
    plt.colorbar()
    plt.title("Ez")

    return

def test_update_v(pic):
    '''test updating the velocity'''

    # test updating particle velocity
    pic.electron_vx = np.zeros(pic.n_particles)   # zero out particle velocity
    pic.electron_vy = np.zeros(pic.n_particles)
    pic.electron_vz = np.zeros(pic.n_particles)
    
    pic.init_x_random()

    # create sinusoidal electric field
    sin_3d = np.zeros(pic.nodes) 
    for i in range(pic.nodes[0]):
        for j in range(pic.nodes[1]):
            for k in range(pic.nodes[2]):
                sin_3d[i][j][k] = np.sin(pic.dx[0] * i ) + \
                                  np.sin(pic.dx[1] * j ) + \
                                  np.sin(pic.dx[2] * k )
    pic.ex = sin_3d
    pic.ey = sin_3d
    pic.ez = sin_3d

    pic.update_v()        # update velocity based on E field

    positions = [pic.electron_x, pic.electron_y, pic.electron_z]
    velocities = [pic.electron_vx, pic.electron_vy, pic.electron_vz]
    coords = ["x", "y", "z"]


    start_fig = 17
    dims = pic.dimensions

    # plot x, y, and z electric field along every pair of axes
    for i in range(dims):
        for j in range(dims): 
    
            plt.figure(start_fig + i * dims + j)
            plt.scatter(positions[i % dims], positions[(i + 1) % dims], c=velocities[j], s = 2)
            plt.title("Velocity {}".format(coords[j]))
            plt.xlabel(coords[i % dims])
            plt.ylabel(coords[(i + 1) % dims])
            plt.colorbar()

    return

def test_electrostatic_energy(pic):
    '''test calculating the electrostatic energy'''

    # set the electric field to all ones
    pic.ex = np.ones(pic.nodes)
    pic.ey = np.ones(pic.nodes)
    pic.ez = np.ones(pic.nodes)

    pic.calc_electrostatic_energy()
    print("Calculated electrostatic energy: ", pic.output["electrostatic_energy"][0])
    print("Expected: ", 0.5 * 3 * pic.xmax * pic.ymax * pic.zmax)

    return

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

    #test_e_flat(pic)
    #test_update_v(pic)
    test_electrostatic_energy(pic)
    #plt.show()
