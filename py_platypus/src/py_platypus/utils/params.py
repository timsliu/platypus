# class for describing a simulation
import numpy as np
import json
import os
import shutil

# full path to platypus home directory
PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")

class Parameters:

    def __init__(self, dims, load_file=None):
        '''inputs: dims - dimensions of the simulation
                   load_file - full path to json file to load parameters from'''
        # dictionary with default parameters
        
        self.params = {
            "name": "default_simulation",
            "seed": 0,              # random seed
            "version": "1.0",
            "length": [2 * np.pi for x in range(dims)],
            "cells": [32 for x in range(dims)],
            "timestep": 0.04, 
            "runtime": 30,
            "dimensions": dims,
            "nppc": 100,
            "dx": None,              # cell size (derived)
            "steps": None,           # total steps (derived)
            "n_particles": None,     # total particles (derived)
            "print_every": 20,       # time steps between printing current step
            "save_every": 100,       # interval between saving data
            "single_stream": {       # defaults for single stream instability
                "stream_v": 1, 
                "stream_frac": 0.5, 
                "stream_width": 0.5
            },
            "landau": {              # defaults for Landau damping
                "amplitude": 0.5,
                "mode": 1
            },
            "two_stream": {          # defaults for two stream instability
                "vpos": 0.5,
                "vneg": -0.5,
                "stream_frac": 1,
                "stream_width": 0.5
            },
        }

        # initialize the parameters from a json file
        if load_file is not None:
            self.params = json.load(open(load_file))

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
        self.set("dx", list(map(lambda a, b: a/b, self["length"], self["cells"])), refresh=False)
        self.set("steps", int(self["runtime"]/self["timestep"]), refresh=False)
        self.set("n_particles", self["nppc"] * np.prod(self["cells"]), refresh=False)

        return

    def save_json(self):
        '''save the parameters as a json file'''
        params = {}
        for key in self.params.keys():
            # convert numpy type to python types 
            if type(self[key]) == type(np.array([])):
                params[key] = list(self[key])
            if type(self[key]) == type(np.int64(0)):
                params[key] = int(self[key])
            else:
                params[key] = self[key]

        json_data = json.dumps(params, indent = 4)
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
