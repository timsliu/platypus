# plotter class for plotting different data

import matplotlib.pyplot as plt
import numpy as np
import vis_util

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim
import params
import utils
import pickle

class Plotter:
    def __init__(self, name):
        #self.params = params
        params = {"name": name}
        self.dims = 1
        self.dpi = 800
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/graphs".format(params["name"]))
        
        self.data_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/data".format(params["name"]))

    def get_files(self, id_str):
        '''get all files with a certain id string'''
        files = os.listdir(self.data_dir)
        
        # filter for files with the ID type
        target_files = list(filter(lambda x: id_str in x, files))
        
        # order alphabetically, which also sorts by step 
        target_files.sort(key=lambda x: int(x.split("_")[1]))

        # parse for the step count
        steps = [int(f.split("_")[1]) for f in target_files]
        print(target_files)
        print(steps)
        
        return target_files, steps

    def plot_electric_field(self):
        '''plot the electric field at multiple timesteps'''
       
        if self.dims == 1:
            E_fields = []
            E_files, steps = self.get_files("E")
           
            # open all the electric field files
            for f in E_files:
                E_fields.append(
                    [pickle.load(open(os.path.join(self.data_dir, f), "rb"))])
            
            # name of output file name 
            graph_file_name = os.path.join(self.out_dir, "E.png")
            
            # plot the subplots
            vis_util.plot_lines(
                graph_file_name,
                E_fields,
                "Position (x)",
                "Electric field",
                ["E field step {}".format(x) for x in steps],
                vis_util.subplot_lines
            )

        #if self.dims == 2
        #    Ex_files, steps = get_files("Ex")
        #    Ey_files, _ = get_files("Ey")



        return

    def plot_energy(self):
        '''plot the kinetic and electrostatic energy as a function of time
        step'''
        
        ee = pickle.load(open(os.path.join(self.data_dir, "ee.p"), "rb"))
        ke = pickle.load(open(os.path.join(self.data_dir, "ke.p"), "rb"))
        ee_file = os.path.join(self.out_dir, "ee.png") 
        ke_file = os.path.join(self.out_dir, "ke.png") 
        combined_file = os.path.join(self.out_dir, "energy.png")

        # plot kinetic and electrostatic energy on one chart
        vis_util.plot_lines(
            combined_file,
            [[ke, ee]], 
            "Time step", 
            "Normalized energy", 
            ["Energy"],
            vis_util.subplot_lines, 
            legend=[["Kinetic energy", "Electrostatic energy"]]
        )

        # plot kinetic energy
        vis_util.plot_lines(
            ee_file, 
            [[ee]], 
            "Time step", 
            "Normalized energy", 
            ["Electrostatic energy"],
            vis_util.subplot_lines, 
        )
        
        # plot potential energy
        vis_util.plot_lines(
            ke_file, 
            [[ke]], 
            "Time step", 
            "Normalized energy", 
            ["Kinetic energy"],
            vis_util.subplot_lines, 
        )


        return

    def plot_density(self):
        return

    def plot_phase(self):
        return

    def plot_all(self):
        return

    def plot_velocity(self):
        return

    def plot_position(self):
        return
