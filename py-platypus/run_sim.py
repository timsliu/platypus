import matplotlib.pyplot as plt
import os
import params
import utils

from pic_1d import *
from pic_2d import *

def run_simulation(pic, params):
    '''helper function for running a simulation and saving the outputs
    inputs: pic - initialized PIC_1D object instance
            params - parameter class''' 

    # initialize the pickle output directory
    params.init_output()
    steps = int(params["runtime"]/params["timestep"])
    dims = params["dimensions"]

    for step in range(steps):
        if step % params["print_every"] == 0:
            print("Step {}".format(step))
        pic.step()                       # step the simulatioin
        pic.calc_electrostatic_energy()  # calc the total electrostatic energy
        pic.calc_kinetic_energy()        # calc the total kinetic energy
        
        # save simulation information
        if step % params["save_every"] == 0 or step == steps - 1:
            if dims == 1: 
                utils.save_pickle("{}/step_{}_ev".format(
                    params.data_dir, step), pic.electron_v)
                utils.save_pickle("{}/step_{}_ex".format(
                    params.data_dir, step), pic.electron_x)
                utils.save_pickle("{}/step_{}_ef".format(
                    params.data_dir, step), pic.e)
            if dims == 2: 
                utils.save_pickle("{}/step_{}_evx".format(
                    params.data_dir, step), pic.electron_vx)
                utils.save_pickle("{}/step_{}_evy".format(
                    params.data_dir, step), pic.electron_vy)
                utils.save_pickle("{}/step_{}_ex".format(
                    params.data_dir, step), pic.electron_x)
                utils.save_pickle("{}/step_{}_ey".format(
                    params.data_dir, step), pic.electron_y)
                utils.save_pickle("{}/step_{}_efx".format(
                    params.data_dir, step), pic.ex)
                utils.save_pickle("{}/step_{}_efy".format(
                    params.data_dir, step), pic.ey)

            utils.save_pickle("{}/step_{}_ne".format(
                params.data_dir, step), pic.ne)
    
    # get the saved quantities and save to pickle
    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    
    utils.save_pickle("{}/ee".format(params.data_dir), ee)
    utils.save_pickle("{}/ke".format(params.data_dir), ke)

    return

def init_pic_by_dims(dims, params):
    '''helper function initializes a PIC with the specified number of
    dimensions'''
    if dims == 1:
        pic = PIC_1D(params)
    elif dims == 2:
        pic = PIC_2D(params)
    else:
        raise ValueError("Expected dimensions 1 or 2; passed: {}".\
                         format(dims))
    return pic

def two_stream(name, dims, param_dict={}):
    '''set up and run a two stream instability
    inputs: name - name of simulation
            dims - number of dimensions
            param_dict - optional dictionary of parameters'''

    # create a simulation object for holding simulation parameters
    sim_params = params.Parameters(dims)
    sim_params.set("name", name) 
    sim_params.set_from_dict(param_dict)

    pic = init_pic_by_dims(dims, sim_params)
    
    # initialize x randomly and create two streams
    pic.init_x_random()
    
    if dims == 2:
        pic.init_v_maxwellian()
    
    pic.init_v_two_stream()
    run_simulation(pic, sim_params)

    return sim_params

def single_stream(name, dims, param_dict={}):
    '''set up and run a single stream simulation demonstrating neutral beam
    injection heating.
    inputs: name - name of simulation
            dims - number of dimensions
            param_dict - optional dictionary of parameters'''
   
    sim_params = params.Parameters(dims)
    sim_params.set("name", name)
    sim_params.set_from_dict(param_dict)

    pic = init_pic_by_dims(dims, sim_params)
    
    # initialize x randomly and create a maxwellian velocity with one stream
    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.init_v_single_stream() 

    run_simulation(pic, sim_params)

    return sim_params

def landau(name, dims, param_dict={}):
    '''set up and run a landau damping simulation
    inputs: name - name of simulation
            dims - number of dimensions
            param_dict - optional dictionary of parameters'''
    
    sim_params = params.Parameters(dims)
    sim_params.set("name", name)
    sim_params.set_from_dict(param_dict)

    pic = init_pic_by_dims(dims, sim_params)
    
    pic.init_x_random()           # initialize random x
    pic.init_v_maxwellian()       # maxwellian velocity distribution
    pic.density_perturbation()    # create charge density perturbation
    
    run_simulation(pic, sim_params)
    
    return sim_params


