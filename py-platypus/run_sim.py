import matplotlib.pyplot as plt
import os
import params
import utils

from pic_1d import *

PRINT_EVERY = 20    # print step number
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

    return

def two_stream(name, vpos, vneg, param_dict={}):
    '''set up and run a two stream instability
    inputs: name - name of simulation
            vpos - velocity of positive stream
            vneg - velocity of negative stream
            params - optional dictionary of simulation parameters'''

    # create a simulation object for holding simulation parameters
    sim_params = params.Parameters()
    sim_params.set("name", name)
    sim_params.set("vpos", vpos)
    sim_params.set("vneg", vneg)
    sim_params.set_from_dict(param_dict)

    # initialize x randomly and two streams
    pic = PIC_1D(sim_params)
    pic.init_x_random()
    pic.init_v_two_beams(vpos, vneg)
    
    run_simulation(pic, sim_params)

    return

def single_stream(nppc):
    '''set up and run a single stream simulation demonstrating neutral beam
    injection heating'''
    length = 2 * np.pi
    cells = 32
    timestep = 0.04
    runtime = 80
    steps = int(runtime/timestep)
    
    # initialize x randomly and two streams
    pic = PIC_1D(cells, length/cells, timestep, cells * nppc, steps)
    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.init_v_single_stream(0.1, 1.5) 

    pickle_names = []
    for step in range(steps):
        print("Step {}".format(step))
        pic.step()
        pic.calc_electrostatic_energy()
        pic.calc_kinetic_energy()
        pic.calc_batch_kinetic_energy()

        # periodically save the electron velocity
        if step % 200 == 0 or step == steps - 1:
            name = "single_stream_npcc_{}_step_{}_ev".format(nppc, step)
            save_pickle(name, pic.electron_v)
            pickle_names.append(name + ".p")

    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    batch_ke = pic.output["batch_ke"]
    save_pickle("single_stream_npcc_{}_extended_ee".format(nppc), ee)
    save_pickle("single_stream_npcc_{}_extended_ke".format(nppc), ke)
    save_pickle("single_stream_npcc_{}_extended_batch_ke".format(nppc), batch_ke)

    return pickle_names


def landau(nppc, amplitude):
    '''set up and run a landau damping simulation'''
    length = 2 * np.pi
    cells = 32
    timestep = 0.04
    runtime = 30
    steps = int(runtime/timestep)

    # initialize x randomly and two streams
    pic = PIC_1D(cells, length/cells, timestep, cells * nppc, steps)
    pic.init_x_random()
    pic.init_v_maxwellian()
    pic.density_perturbation(amplitude, 1)
    
    # step through the simulation
    for step in range(steps):
        print("Step {}".format(step))
        pic.step()                       # update particle positions
        pic.calc_electrostatic_energy()  # save the total electrostatic energy
        pic.calc_kinetic_energy()        # save the total kinetic energy
        if step % 200 == 0 or step == steps - 1:
            save_pickle("landau_amplitude_{}_npcc_{}_step_{}_ev".format(amplitude, nppc, step), pic.electron_v)
            save_pickle("landau_amplitude_{}_npcc_{}_step_{}_ex".format(amplitude, nppc, step), pic.electron_x)
            save_pickle("landau_amplitude_{}_npcc_{}_step_{}_ne".format(amplitude, nppc, step), pic.ne)
            save_pickle("landau_amplitude_{}_npcc_{}_step_{}_e".format(amplitude, nppc, step), pic.e)

    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    save_pickle("landau_amplitude_{}_npcc_{}_ee".format(amplitude, nppc), ee)
    save_pickle("landau_amplitude_{}_npcc_{}_ke".format(amplitude, nppc), ke)

    return


