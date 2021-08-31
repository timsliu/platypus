'''
io_utils.py

miscellaneous io utility functions related to reading/writing to files and
other I/O operations
'''

import pickle

def save_pickle(name, data):
    '''save some data to a pickle file'''
    file_name = name + ".p" 
    pickle.dump(data, open(file_name, "wb"))
    return

