import matplotlib.pyplot as plt
import pickle
from pic_1d import *
import os

PRINT_EVERY = 20    # print step number
SAVE_EVERY = 100    # save simulation info

def run_simulation(steps, pic, arg_str, sim_str):
    '''helper function for running a simulation and saving the outputs
    inputs: steps - number of steps to run simulation for
            pic - initialized PIC_1D object instance
            arg_str - string describing the arguments to be appended to pickle
            sim_str - string with name of simulation'''

    for step in range(steps):
        if step % PRINT_EVERY == 0:
            print("Step {}".format(step))
        
        pic.step()                       # step the simulatioin
        pic.calc_electrostatic_energy()  # calc the total electrostatic energy
        pic.calc_kinetic_energy()        # calc the total kinetic energy

        # save simulation information
        if step % SAVE_EVERY == 0 or step == steps - 1:
            save_pickle("{}_step_{}__ev".format(sim_str, step, arg_str), pic.electron_v)
            save_pickle("{}_step_{}__ex".format(sim_str, step, arg_str), pic.electron_x)
            save_pickle("{}_step_{}__ne".format(sim_str, step, arg_str), pic.ne)
            save_pickle( "{}_step_{}__e".format(sim_str, step, arg_str), pic.e)

    # get the saved quantities and save to pickle
    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    save_pickle("{}_{}_ee".format(sim_str, arg_str), ee)
    save_pickle("{}_{}_ke".format(sim_str, arg_str), ke)

    return

def two_stream(name, vpos, vneg, params={}):
    '''set up and run a two stream instability
    inputs: name - name of simulation
            vpos - velocity of positive stream
            vneg - velocity of negative stream
            params - optional dictionary of simulation parameters'''

    # create a simulation object for holding simulation parameters
    sim_params = Simulation()
    sim_params.set("name", name)
    sim_params.set("vpos", vpos)
    sim_params.set("vneg", vneg)
    sim_params.set_from_dict(params)

    # initialize x randomly and two streams
    pic = PIC_1D(sim_params)
    pic.init_x_random()
    pic.init_v_two_beams(vpos, vneg)
    
    vpos_str = str(vpos).replace(".", "dot")
    vneg_str = str(vneg).replace(".", "dot")
    sim_str = "two_stream"
    arg_str = 


    run_simulation(steps, pic, arg_str, sim_str)

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


def plot_many():
    vel_files = [
        "landau_npcc_1000_step_0_ev.p",
        "landau_npcc_1000_step_200_ev.p",
        "landau_npcc_1000_step_400_ev.p",
        "landau_npcc_1000_step_600_ev.p",
        "landau_npcc_1000_step_749_ev.p",
    ]
    
    pos_files = [
        "landau_npcc_1000_step_0_ex.p",
        "landau_npcc_1000_step_200_ex.p",
        "landau_npcc_1000_step_400_ex.p",
        "landau_npcc_1000_step_600_ex.p",
        "landau_npcc_1000_step_749_ex.p",
    ]
    
    #for f in vel_files:
    #    plot_histogram(f, "Velocity", "Count", f[0:-2])

    #for f in pos_files:
    #    plot_histogram(f, "Position", "Count", f[0:-2])

    #for i in range(len(pos_files)):
    #    title = "Pos vs Velocity " + pos_files[i][:-5]
    #    plot_scatter(pos_files[i], vel_files[i], "Position", "Velocity", title)

    files = os.listdir()

    for f in files:
        #if "ee.p" in f and "landau" in f and ".png" not in f:
        #    print(f)
        #    plot_line(f, "Timestep", "Electrostatic Energy", "Electrostatic energy vs. time")
        #if "ke.p" in f and "landau" in f and ".png" not in f:
        #    print(f)
        #    plot_line(f, "Timestep", "Kinetic Energy", "Kinetic energy vs. time")
        #if "_ne.p" in f and ".png" not in f and "_50_" in f:
        #    print(f)
        #    plot_line(f, "Timestep", "Electron density", "Electron density vs time")
        #if "_e.p" in f and ".png" not in f and "_50_" in f:
        #    print(f)
        #    plot_line(f, "Timestep", "Electric Field", "Electric field vs time")
        if "amplitude" in f and "_ee.p" in f and ".png" not in f:
            print(f)
            plot_line(f, "Timestep", "Electrostatic Energy", "Electrostatic energy vs. time")

    return
def many_two_stream():
    two_stream(0.5, -0.5, 25)
    two_stream(0.5, -0.5, 100)
    two_stream(0.5, -0.5, 200)

def many_landau():
    landau(100, 0.5)
    #landau(250)
    #landau(500)

if __name__== "__main__":
    pickle_names = single_stream(200)
    print(len(pickle_names))
    
    #pickle_names = [
    #    "single_stream_npcc_200_step_0_ev.p", 
    #    "single_stream_npcc_200_step_200_ev.p", 
    #    "single_stream_npcc_200_step_400_ev.p", 
    #    "single_stream_npcc_200_step_600_ev.p", 
    #    "single_stream_npcc_200_step_800_ev.p", 
    #    "single_stream_npcc_200_step_999_ev.p", 
    #]

    #plot_histograms(
    #    pickle_names, 
    #    "Velocity", 
    #    "Number of particles", 
    #    "Velocity distribution", 
    #    "single_stream_ev.png")
   
    b = "single_stream_npcc_200_extended_batch_ke.p"
    k = "single_stream_npcc_200_extended_ke.p"
    e = "single_stream_npcc_200_extended_ee.p"
    plot_line(b, "Timestep", "Batch Kinetic Energy", "Kinetic energy vs. time")
    plot_line(k, "Timestep", "Kinetic Energy", "Kinetic energy vs. time")
    plot_line(e, "Timestep", "Electrostatic Energy", "Electrostatic energy vs. time")


