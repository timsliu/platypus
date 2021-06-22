import matplotlib.pyplot as plt
import os
import params
import utils

from pic_1d import *

PRINT_EVERY = 100    # print step number
SAVE_EVERY = 100    # save simulation info

def run_simulation(pic, params):
    '''helper function for running a simulation and saving the outputs
    inputs: pic - initialized PIC_1D object instance
            params - parameter class''' 

    # initialize the pickle output directory
    params.init_output()
    steps = int(params["runtime"]/params["timestep"])

    for step in range(steps):
        if step % PRINT_EVERY == 0:
            print("Step {}".format(step))
        
        pic.step()                       # step the simulatioin
        pic.calc_electrostatic_energy()  # calc the total electrostatic energy
        pic.calc_kinetic_energy()        # calc the total kinetic energy

        # save simulation information
        if step % SAVE_EVERY == 0 or step == steps - 1:
            utils.save_pickle("{}/step_{}_ev".format(
                params.data_dir, step), pic.electron_v)
            utils.save_pickle("{}/step_{}_ex".format(
                params.data_dir, step), pic.electron_x)
            utils.save_pickle("{}/step_{}_ne".format(
                params.data_dir, step), pic.ne)
            utils.save_pickle("{}/step_{}_e".format(
                params.data_dir, step), pic.e)

    # get the saved quantities and save to pickle
    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    
    utils.save_pickle("{}/ee".format(params.data_dir), ee)
    utils.save_pickle("{}/ke".format(params.data_dir), ke)

    diff_ee = np.array(ee) - ee[0]
    max_diff_ee = max(abs(diff_ee))
    max_index = list(abs(diff_ee)).index(max_diff_ee)
    
    delta_ee = diff_ee[max_index]
    delta_ke = ke[max_index] - ke[0]

    ratio = delta_ee/delta_ke

    return params, delta_ee, delta_ke, ratio

def two_stream(name, vpos, vneg, param_dict={}):
    '''set up and run a two stream instability
    inputs: name - name of simulation
            vpos - velocity of positive stream
            vneg - velocity of negative stream
            params - optional dictionary of simulation parameters'''

    # create a simulation object for holding simulation parameters
    sim_params = params.Parameters(1)
    sim_params.set("name", name)
    sim_params.set("vpos", vpos)
    sim_params.set("vneg", vneg)
    sim_params.set_from_dict(param_dict)

    # initialize x randomly and two streams
    pic = PIC_1D(sim_params)
    pic.init_x_random()
    pic.init_v_two_beams(vpos, vneg)
    
    return run_simulation(pic, sim_params)


def single_stream(name, stream_v, stream_frac, param_dict={}):
    '''set up and run a single stream simulation demonstrating neutral beam
    injection heating'''
   
    sim_params = params.Parameters(1)
    sim_params.set("name", name)
    sim_params.set("stream_v", stream_v)
    sim_params.set("stream_frac", stream_frac)
    sim_params.set_from_dict(param_dict)

    # initialize x randomly and create a maxwellian velocity with one stream
    pic = PIC_1D(sim_params)
    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.init_v_single_stream(stream_frac, stream_v) 

    run_simulation(pic, sim_params)

def landau(name, amplitude, mode, param_dict={}):
    '''set up and run a landau damping simulation'''
    
    sim_params = params.Parameters(1)
    sim_params.set("name", name)
    sim_params.set("mode", mode)
    sim_params.set("amplitude", amplitude)
    sim_params.set_from_dict(param_dict)

    pic = PIC_1D(sim_params)
    pic.init_x_random()              # initialize random x
    pic.init_v_maxwellian()          # maxwellian velocity distribution
    pic.density_perturbation(amplitude, mode) # create charge density perturbation
    
    run_simulation(pic, sim_params)
    
    return


