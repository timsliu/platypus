# class for describing a simulation
import numpy as np
import json
import os
import shutil

# full path to platypus home directory
PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")

class Parameters:

    def __init__(self, dims):
        # dictionary with default parameters
        # TODO have different defaults based on dimension
        self.params = {
            "name": "default-simulation",
            "seed": 0,              # random seed
            "version": "1.0",
            "length": None,
            "cells": None,
            "timestep": 0.04, 
            "runtime": 30,
            "dimensions": 1,
            "nppc": 50,
            "vpos": None,            # velocity of positive stream (2 stream)
            "vneg": None,            # velocity of negative stream (2 stream)
            "mode": None,            # number of density waves
            "amplitude": None,       # amplitude of charge perturbation
            "dx": None,              # cell size (derived)
            "steps": None,           # total steps (derived)
            "n_particles": None,     # total particles (derived)
        }

        # initialize defaults based on number of dimensions
        if dims == 1:
            self.init_1d()

        if dims == 1:
            self.init_2d()
       
        # output directories
        self.out_dir = os.path.join(
            PLATYPUS_HOME, "py-platypus/out/{}".format(self["name"])
        )
        self.data_dir = os.path.join(self.out_dir, "data")
        self.graph_dir = os.path.join(self.out_dir, "graph")

    def init_1d(self):
        '''initialize default parameters for 1d simulation'''
        params = {
            "length": [2 * np.pi],
            "cells": [32],
            "dimensions": 1
        }

        self.set_from_dict(params)
        return

    def init_2d(self):
        '''initialize default parameters for 2d simulation'''
        params = {
            "length": [2 * np.pi, 2 * np.pi],
            "cells": [32, 32],     # matrix notation row, col
            "dimensions": 2
        }

        self.set_from_dict(params)
        return

    def set_from_dict(self, params):
        '''set parameters from a passed dictionary. Parameters not
        specified in params will be left as the default value'''

        # set the object's parameters to the values in the passed dict
        print(params) 
        for param in params.keys():
            print(param) 
            self.set(param, params[param], refresh=False)

        self.refresh()
        return

    def set(self, param, value, refresh=True):
        '''set a parameter to a specificied value'''

        # check valid parameter
        if param not in self.params.keys():
            raise ValueError(
                "Parameter {} not a simulation parameters".format(param))

        # TODO additional checks on parameters
        self.params[param] = value   # set the value
        print(param, value)
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
        print(self["length"])
        print(self["cells"])
        self.set("dx", list(map(lambda a, b: a/b, self["length"], self["cells"])), refresh=False)
        self.set("steps", int(self["runtime"]/self["timestep"]), refresh=False)
        # TODO this is dependent on the dimension
        self.set("n_particles", self["nppc"] * np.prod(self["cells"]), refresh=False)

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
