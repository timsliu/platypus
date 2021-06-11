import matplotlib.pyplot as plt
import pickle
from pic_1d import *
import os

def save_pickle(name, data):
    '''save some data to a pickle file'''
    file_name = name + ".p" 
    pickle.dump(data, open(file_name, "wb"))
    return

def plot_line(filename, x_axis, y_axis, title):
    '''plot a histogram of some data'''

    # load data from pickle
    data = pickle.load(open(filename, "rb"))
    
    plt.figure()
    plt.plot(data)
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.grid(True)
    #plt.yscale("log")
    plt.savefig(filename.replace(".p", ".png"), dpi=800)

    return

def plot_histograms(filenames, x_axis, y_axis, title, fig_name):
    '''plot a histogram of some data'''
    rows = 2
    cols = 3
    # load data from pickle
    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    timesteps = [0, 200, 400, 600, 800, 999] 

    for i, filename in enumerate(filenames):
        print(filename)
        data = pickle.load(open(filename, "rb"))
        col = i % cols
        row = int(np.floor(i/cols))
        print(row, col)
        
        axs[row, col].hist(data, bins=20)
        axs[row, col].set_title("Timestep: {}".format(timesteps[i]))
        axs[row, col].set_ylim(0, 850)
        axs[row, col].grid(True)
    
    # set the axis labels
    for ax in axs.flat:
        ax.set(xlabel = x_axis, ylabel = y_axis)

    # only have axis titles on the outer edge
    for ax in axs.flat:
        ax.label_outer()
    
    plt.savefig(fig_name, dpi=800)


    return

def plot_histogram(filename, x_axis, y_axis, title):
    '''plot a histogram of some data'''

    # load data from pickle
    data = pickle.load(open(filename, "rb"))
    
    plt.figure()
    plt.hist(data, bins=20)
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.yscale("log")
    plt.grid(True)
    plt.savefig(filename.replace(".p", ".png"), dpi=800)

    return

def plot_scatter(x_file, y_file, x_axis, y_axis, title):
    '''plot a histogram of some data'''

    # load data from pickle
    x_data = pickle.load(open(x_file, "rb"))
    y_data = pickle.load(open(y_file, "rb"))
    
    plt.figure()
    plt.scatter(x_data, y_data, s=1)
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.grid(True)
    plt.savefig(title.replace(" ", "_") + ".png", dpi=800)

    return

####### Functions for running simulation #######

def two_stream(vpos, vneg, nppc):
    '''set up and run a two stream instability'''
    length = 2 * np.pi
    cells = 32
    timestep = 0.04
    runtime = 30
    steps = int(runtime/timestep)

    # initialize x randomly and two streams
    pic = PIC(cells, length/cells, timestep, cells * nppc, steps)
    pic.init_x_random()
    pic.init_v_two_beams(vpos, vneg)
    
    vpos_str = str(vpos).replace(".", "dot")
    vneg_str = str(vneg).replace(".", "dot")

    # step through the simulation
    for step in range(steps):
        print("Step {}".format(step))
        pic.step()                       # update particle positions
        pic.calc_electrostatic_energy()  # save the total electrostatic energy
        pic.calc_kinetic_energy()        # save the total kinetic energy

        if step % 200 == 0 or step == steps - 1:
            save_pickle("two_stream_npcc_{}_step_{}_v_pos_{}_v_neg_{}_ev".format(nppc, step, vpos_str, vneg_str), pic.electron_v)
            save_pickle("two_stream_npcc_{}_step_{}_v_pos_{}_v_neg_{}_ex".format(nppc, step, vpos_str, vneg_str), pic.electron_x)
            save_pickle("two_stream_npcc_{}_step_{}_v_pos_{}_v_neg_{}_ne".format(nppc, step, vpos_str, vneg_str), pic.ne)
            save_pickle("two_stream_npcc_{}_step_{}_v_pos_{}_v_neg_{}_e".format(nppc, step, vpos_str, vneg_str), pic.e)

    ee = pic.output["electrostatic_energy"]
    ke = pic.output["kinetic_energy"]
    save_pickle("two_stream_npcc_{}_v_pos_{}_v_neg_{}_ee".format(nppc, vpos_str, vneg_str), ee)
    save_pickle("two_stream_npcc_{}_v_pos_{}_v_neg_{}_ke".format(nppc, vpos_str, vneg_str), ke)
    
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
    pic = PIC(cells, length/cells, timestep, cells * nppc, steps)
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
    pic = PIC(cells, length/cells, timestep, cells * nppc, steps)
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

def test_n():
    length = 2 * np.pi
    cells = 32
    timestep = 0.04
    runtime = 30
    nppc = 50
    steps = int(runtime/timestep)

    # initialize x randomly and two streams
    pic = PIC(cells, length/cells, timestep, cells * nppc, steps)
    pic.test_n()

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


