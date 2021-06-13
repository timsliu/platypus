# class for describing a simulation
import numpy as np
import json
import os
import shutil

# full path to platypus home directory
PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")

class Parameters:

    def __init__(self):
        # dictionary with default parameters
        self.params = {
            "name": "default-simulation",
            "seed": 0,              # random seed
            "version": "1.0",
            "length": 2 * np.pi,
            "cells": 32,
            "timestep": 0.04, 
            "runtime": 30,
            "dimensions": 1,
            "nppc": 50,
            "vpos": None,
            "vneg": None,
            "dx": None,              # cell size (derived)
            "steps": None,           # total steps (derived)
            "n_particles": None,     # total particles (derived)
        }
       
        # output directories
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}".format(self["name"])
        )
        self.data_dir = os.path.join(self.out_dir, "data")
        self.graph_dir = os.path.join(self.out_dir, "graph")

    def set_from_dict(self, params):
        '''set parameters from a passed dictionary. Parameters not
        specified in params will be left as the default value'''

        # set the object's parameters to the values in the passed dict
        for param in params.keys():
            self.set(param, params[param])

    def set(self, param, value, refresh=True):
        '''set a parameter to a specificied value'''

        # check valid parameter
        if param not in self.params.keys():
            raise ValueError(
                "Parameter {} not a simulation parameters".format(param))

        # TODO additional checks on parameters
        self.params[param] = value   # set the value

        # make sure derived fields are up to date
        if refresh:
            self.refresh()
   
        return

    def refresh(self):
        '''update the derived self.params'''
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}".format(self["name"])
        )
        self.data_dir = os.path.join(self.out_dir, "data")
        self.graph_dir = os.path.join(self.out_dir, "graph")
   
        self.set("dx", self["length"]/self["cells"], refresh=False)
        self.set("steps", int(self["runtime"]/self["timestep"]), refresh=False)
        self.set("n_particles", self["nppc"] * self["cells"], refresh=False)

        return

    def save_json(self):
        '''save the parameters as a json file'''
        
        json_data = json.dumps(self.params, indent = 4)
        with open(os.path.join(self.out_dir, "params.json"), "w") as f:
            f.write(json_data)
        
        return

    def init_output(self):
        '''create a folder for output simulation information'''
        
        # create output folder for the simulation
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)
        
        os.mkdir(self.out_dir)                           # create output directory
        os.mkdir(os.path.join(self.out_dir, "data"))     # subfolder for pickles
        os.mkdir(os.path.join(self.out_dir, "graphs"))   # subfolder for data

        self.save_json()      # dump the parameters to file

    def __getitem__(self, key):
        '''overload the braces [] operator to look up fields in the params
        attribute '''
        
        return self.params[key]
