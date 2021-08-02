import py_platypus as plat
import numpy as np
import matplotlib.pyplot as plt

NUM_TESTS = 5

def test_random_x(pic):
    pic.init_x_random()
    plt.figure(1)
    plt.scatter(pic.electron_x, pic.electron_y, s=0.05)

    return

def test_init_E(pic):
    '''test calculating the inital electric field'''

    pic.init_v_maxwellian()
    pic.density_perturbation()
    pic.init_E()

    plt.figure(2)
    plt.imshow(pic.ex)
    plt.title("Ex field") 
    plt.colorbar()

    plt.figure(3) 
    plt.imshow(pic.ey)
    plt.title("Ey field") 
    plt.colorbar()
    
    plt.figure(4)
    plt.imshow(pic.Ex_edges)
    plt.title("Ex field edges") 
    plt.colorbar()

    plt.figure(5) 
    plt.imshow(pic.Ey_edges)
    plt.title("Ey field edges") 
    plt.colorbar()

    return

def test_update_B(pic):
    '''test calculating the curl of a 2D vector field'''

    # set up x vectors
    rows, cols = pic.Ex_edges.shape
    for i in range(rows):
        for j in range(cols):
            x = i - rows/2
            y = j - cols/2
            
            pic.Ex_edges[i][j] = y 
   
    # set up y vectors
    rows, cols = pic.Ey_edges.shape
    for i in range(rows):
        for j in range(cols):
            x = i - rows/2
            y = j - cols/2
            
            pic.Ey_edges[i][j] = -x 

    pic.calc_B_update()
    plt.figure(6)
    plt.imshow(pic.delta_Bz)
    plt.title("B update for uniform curl")

def test_interpolate(pic):
    '''test helper function for interpolating field properties at particles'''

    print("==== Testing Interpolate ====")
    corners = [0, 1, 0, 1]     # [x0, x1, y0, y1] coordinates of corners
    x_ns = [0.5, 0.75]         # list of x_n particle positions
    y_ns = [0.5, 0.75]         # list of y_n particle positions
    values = [1, 2, 3, 4]      # field values at the four corners

    # expected values
    expected = [0.25 * 1 + 0.25 * 2 + 0.25 * 3 + 0.25 * 4, 
                1/16 * 1 + 3/16 * 2 + 3/16 * 3 + 9/16 * 4]
    
    # scale expected values by 1/ area of cells 
    expected = 1 / np.prod(pic.dx) * np.array(expected)

    # iterate through test cases
    for i in range(len(x_ns)):
        int_value = pic.interpolate(x_ns[i], y_ns[i], corners, values) 
        print("Expected: {} Actual: {}".format(expected[i], int_value))

    return

def test_interpolate_ex(pic):
    '''test interpolating the electric field at each particle from the field
    values'''

    print("==== Testing Interpolate Ex field ====")
    pic.Ex_edges = np.random.rand(pic.Ex_edges.shape[0], pic.Ex_edges.shape[1])
    pic.init_x_random() 
    pic.interpolate_e()

    print("Average Ex field: {} Average Ex particle: {}".format(
        np.mean(pic.Ex_edges), np.mean(pic.e_particle[:, 0])))
  
    expected = np.zeros(NUM_TESTS)
    for i in range(NUM_TESTS):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[i] = (1 + offset) * dx 
        pic.electron_y[i] = (1 - offset) * dy 

        short = 0.5 - offset
        expected[i] = (pic.Ex_edges[0][0] * offset * short +\
                       pic.Ex_edges[0][1] * (1 - short) * offset +\
                       pic.Ex_edges[1][0] * (1 - offset) * short +\
                       pic.Ex_edges[1][1] * (1 - offset) * (1 - short)) 
    pic.interpolate_e()
   
    for i in range(NUM_TESTS):
        print("Expected: {} Actual: {}".format(expected[i], pic.e_particle[i][0]))

    return

def test_interpolate_ey(pic):
    '''test interpolating the electric field at each particle from the field
    values'''

    print("==== Testing Interpolate Ey field ====")
    pic.Ey_edges = np.indices(pic.Ey_edges.shape)[0]
    pic.init_x_random() 
    pic.interpolate_e()

    print("Average Ey field: {} Average Ey particle: {}".format(
        np.mean(pic.Ey_edges), np.mean(pic.e_particle[:, 1])))
    
    expected = np.zeros(NUM_TESTS)
    for i in range(NUM_TESTS):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[i] = (1 + offset) * dx 
        pic.electron_y[i] = (1 - offset) * dy 

        short = 0.5 - offset
        expected[i] = (pic.Ex_edges[0][0] * offset * short +\
                       pic.Ex_edges[0][1] * (1 - short) * offset +\
                       pic.Ex_edges[1][0] * (1 - offset) * short +\
                       pic.Ex_edges[1][1] * (1 - offset) * (1 - short)) 
    pic.interpolate_e()
   
    for i in range(NUM_TESTS):
        print("Expected: {} Actual: {}".format(expected[i], pic.e_particle[i][0]))

    return

def test_interpolate_b(pic):
    '''test interpolate the magnetic field at each particle from the field
    values'''

    print("==== Testing Interpolate B field ====")
    pic.Bz = np.indices(pic.Bz.shape)[0]
    pic.init_x_random() 
    pic.interpolate_b() 
    
    print("Average B field: {} Average B particle: {}".format(
        np.mean(pic.Bz), np.mean(pic.b_particle)))
   
    values = [[4, 1], [8, 2]]

    pic.Bz[0][0] = values[0][0]
    pic.Bz[0][1] = values[0][1]
    pic.Bz[1][0] = values[1][0]
    pic.Bz[1][1] = values[1][1]

    for i in range(5):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[0] = (1 + offset) * dx 
        pic.electron_y[0] = (1 - offset) * dy 

        short = 0.5 - offset
        expected = (pic.Bz[0][0] * (1 - short) * short +\
                    pic.Bz[0][1] * (1 - short) * (1 - short) +\
                    pic.Bz[1][0] * short * short +\
                    pic.Bz[1][1] * (1 - short) * short)
        pic.interpolate_b()
    
        print("Expected: {} Actual: {}".format(expected, pic.b_particle[0]))
    return

if __name__ == "__main__":
    # set up parameters 
    params = {
        "length": [2 * np.pi, 2 * np.pi],
        "cells": [32, 32],
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

    #test_random_x(pic)
    #test_init_E(pic)
    #test_update_B(pic)
    test_interpolate(pic)
    test_interpolate_b(pic) 
    test_interpolate_ex(pic) 
    test_interpolate_ey(pic) 
    plt.show()
