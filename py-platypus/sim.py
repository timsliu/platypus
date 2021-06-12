# class for describing a simulation
import numpy
import json

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

    def save_json(self, path = "."):
        '''save the parameters as a json file'''

