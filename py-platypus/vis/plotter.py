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
        self.dpi = 800
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/graphs".format(params["name"]))
        
        self.data_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}/data".format(params["name"]))

    def plot_electric_field(self):
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
            combined_file, [ke, ee], "Time step", 
            "Normalized energy", "Energy", 
            legend = ["Kinetic energy", "Electrostatic energy"]
        )

        # plot kinetic energy
        vis_util.plot_lines(
            ee_file, [ee], "Time step", 
            "Normalized energy", "Electrostatic energy"
        )
        
        # plot potential energy
        vis_util.plot_lines(
            ke_file, [ke], "Time step", 
            "Normalized energy", "Kinetic energy"
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
