# class for describing a simulation
import numpy
import json
import os

# full path to platypus home directory
PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")

class Parameters:

    def __init__(self):
        # dictionary with default parameters
        self.params = {
            "name": "default-simulation",
            "length": 2 * np.pi,
            "cells": 32,
            "timestep": 0.04, 
            "runtime": 30,
            "dimensions": 1,
            "nppc": 50,
            "vpos": None,
            "vneg": None,
        }
        
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

    def set(self, param, value):
        '''set a parameter to a specificied value'''

        # check valid parameter
        if param not in self.params.keys():
            raise ValueError("Parameter {} not a simulation parameters".format(param)

        # TODO additional checks on parameters
        self.params[param] = value   # set the value

    def save_json(self):
        '''save the parameters as a json file'''
        
        json_data = json.dumps(self.params, indent = 4)
        with open(os.path.join(self.out_dir, "params.json")) as f:
            f.write(json_data)

    def init_output():
        '''create a folder for output simulation information'''
        # create output folder for the simulation 
        
        os.mkdir(self.out_dir)                           # create output directory
        os.mkdir(os.path.join(self.out_dir, "data"))     # subfolder for pickles
        os.mkdir(os.path.join(self.out_dir, "graphs"))   # subfolder for data

        self.save_json(out_dir)      # dump the parameters to file

    def __getitem__(self, key):
        '''overload the braces [] operator to look up fields in the params
        attribute '''
        
        return self.params[key]
