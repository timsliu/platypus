import matplotlib.pyplot as plt
import os

import py_platypus as plat

from py_platypus.models.pic_1d import PIC_1D as PIC_1D 
from py_platypus.models.pic_2d import PIC_2D as PIC_2D 
from py_platypus.models.pic_3d import PIC_3D as PIC_3D 
from py_platypus.utils.params import Parameters as Parameters

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
        if hasattr(pic, "bz"): 
            pic.calc_magnetic_energy()       # calc the total magnetic energy
        
        # save simulation information
        if step % params["save_every"] == 0 or step == steps - 1:
            if dims == 1: 
                plat.io_utils.save_pickle("{}/step_{}_ev".format(
                    params.data_dir, step), pic.electron_v)
                plat.io_utils.save_pickle("{}/step_{}_ex".format(
                    params.data_dir, step), pic.electron_x)
                plat.io_utils.save_pickle("{}/step_{}_ef".format(
                    params.data_dir, step), pic.e)
            if dims == 2: 
                plat.io_utils.save_pickle("{}/step_{}_evx".format(
                    params.data_dir, step), pic.electron_vx)
                plat.io_utils.save_pickle("{}/step_{}_evy".format(
                    params.data_dir, step), pic.electron_vy)
                plat.io_utils.save_pickle("{}/step_{}_ex".format(
                    params.data_dir, step), pic.electron_x)
                plat.io_utils.save_pickle("{}/step_{}_ey".format(
                    params.data_dir, step), pic.electron_y)
                if hasattr(pic, "ey_edges"):
                    plat.io_utils.save_pickle("{}/step_{}_efx".format(
                        params.data_dir, step), pic.ex_edges)
                    plat.io_utils.save_pickle("{}/step_{}_efy".format(
                        params.data_dir, step), pic.ey_edges)
                else:
                    plat.io_utils.save_pickle("{}/step_{}_efx".format(
                        params.data_dir, step), pic.ex)
                    plat.io_utils.save_pickle("{}/step_{}_efy".format(
                        params.data_dir, step), pic.ey)
            if hasattr(pic, "bz"):
                plat.io_utils.save_pickle("{}/step_{}_bz".format(
                    params.data_dir, step), pic.bz)

            plat.io_utils.save_pickle("{}/step_{}_ne".format(
                params.data_dir, step), pic.ne)
    
    # get the saved quantities and save to pickle
    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    
    plat.io_utils.save_pickle("{}/ee".format(params.data_dir), ee)
    plat.io_utils.save_pickle("{}/ke".format(params.data_dir), ke)
    
    if hasattr(pic, "bz"):
        me = pic.output["magnetic_energy"]
        plat.io_utils.save_pickle("{}/me".format(params.data_dir), me)

    return

def init_pic_by_dims(dims, params):
    '''helper function initializes a PIC with the specified number of
    dimensions'''
    if dims == 1:
        pic = PIC_1D(params)
    elif dims == 2:
        pic = PIC_2D(params)
    elif dims == 3:
        pic = PIC_3D(params)
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
    sim_params = Parameters(dims)
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
   
    sim_params = Parameters(dims)
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
    
    sim_params = Parameters(dims)
    sim_params.set("name", name)
    sim_params.set_from_dict(param_dict)

    pic = init_pic_by_dims(dims, sim_params)
    
    pic.init_x_random()           # initialize random x
    pic.init_v_maxwellian()       # maxwellian velocity distribution
    pic.density_perturbation()    # create charge density perturbation
    
    run_simulation(pic, sim_params)
    
    return sim_params


