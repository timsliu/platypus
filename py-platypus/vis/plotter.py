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
    def __init__(self, name, params):
        self.params = params
        
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/graphs".format(self.params["name"]))
        
        self.data_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/data".format(self.params["name"]))

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

    def plot_series(self, suffix, out_name, x_axis, y_axis, title, subplotter):
        '''helper function to plot a data series across several time steps
        inputs: suffix - identifier at end of pickle files to parse
                out_name - name of output file with .png extension
                x_axis - x_axis for charts
                y_axis - y_axis for charts
                title - title for each subplot
                subplotter - function for plotting a subplot'''

        values = []
        files, steps = self.get_files(suffix)
        
        # open all the files
        for f in files:
            values.append(
                [pickle.load(open(os.path.join(self.data_dir, f), "rb"))])
        
        # name of output file name 
        graph_file_name = os.path.join(self.out_dir, out_name)
        
        # plot the subplots
        vis_util.plot_lines(
            graph_file_name,
            values, 
            x_axis,
            y_axis,
            ["{} step {}".format(title, x) for x in steps],
            subplotter 
        )

    def plot_electric_field(self):
        '''plot the electric field at multiple timesteps'''
      
        # call helper for single dimension
        if self.params["dimensions"] == 1:
            self.plot_series("ef", "ef.png", "Position (x)", 
                "Electric field", "E field", vis_util.subplot_lines) 

        # call helper for two dimensions
        if self.params["dimensions"] == 2:
            self.plot_series("efx", "efx.png", "Position (x)", 
                "Position (y)", "E field", vis_util.subplot_grid) 
            
            self.plot_series("efy", "efy.png", "Position (x)", 
                "Position (y)", "E field", vis_util.subplot_grid) 

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
