# example of setting up, running, and graphing a two stream instability in
# one dimension

import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter

if __name__ == "__main__":
    
    # override default simulation parameters
    param_dict = {
        "runtime": 40,    # time to run simulation for
        "timestep": 0.04, # size of each timestep
        "save_every": 4,  # save outputs at every fourth step
    }

    # setup and run a 1D two stream simulation using the parameters
    # specified in param_dict. This method takes care of initialization
    # and runs the simulation
    plat.run_sim.two_stream("two_stream_1d", 1, param_dict=param_dict)
   
    # load the full parameters object from a json
    param_json = os.path.join(plat.PLATYPUS_HOME, "py_platypus/out/two_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
   
    # create instance of the Plotter class for generating charts
    plotter = Plotter("two_stream_1d", params)
    plotter.add_animation()   # outputs should be animations

    plotter.plot_velocity()   # plot the velocity distribution
    plotter.plot_phase()      # plot the phase chart (velocity and position)
    plotter.plot_density()    # plot electron density
    plotter.plot_energy()     # plot energy
