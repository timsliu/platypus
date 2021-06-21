# test the dependency of the ratio of KE to PE on various simulation
# parameters

import sys
import os
import numpy as np
import pickle

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim
RUNTIME = 20


if __name__ == "__main__":
    # global lists 
    delta_ee_list = []
    delta_ke_list = []
  
    # vary nppc
    nppc_list = [10, 20, 40, 80, 160]
    nppc_ratio = []
    nppc_result = []
    
    for nppc in nppc_list:
        print("\n=== Nppc: {} ===".format(nppc)) 
        for i in range(5):
            
            params = {"runtime": RUNTIME,
                      "nppc": nppc}
            
            # run simulation 
            full_params, delta_ee, delta_ke, ratio = run_sim.two_stream("two-stream_many", 0.5, -0.5, param_dict=params)
  
            delta_ee_list.append(delta_ee)
            delta_ke_list.append(delta_ke) 
            nppc_ratio.append(ratio)
            nppc_result.append(nppc)
    
    pickle.dump(nppc_ratio, open("nppc_ratio.p", "wb"))
    pickle.dump(nppc_result, open("nppc_result.p", "wb"))

    # vary number of cells and length, keeping dx constant
    cell_list = [8, 16, 32, 64]
    len_list = [np.pi, 2 * np.pi, 4 * np.pi, 8 * np.pi]

    cells_ratio = []
    cells_result = []
    
    for j in range(len(cell_list)):
        print("\n=== Cells: {} ===".format(cell_list[j])) 
        for i in range(5):
            
            params = {"runtime": RUNTIME,
                      "cells": [cell_list[j]],
                      "length": [len_list[j]]}
            
            # run simulation 
            full_params, delta_ee, delta_ke, ratio = run_sim.two_stream(
                "two-stream_many", 0.5, -0.5, param_dict=params)
   
            delta_ee_list.append(delta_ee)
            delta_ke_list.append(delta_ke) 
            cells_ratio.append(ratio)
            cells_result.append(cell_list[j])
    
    pickle.dump(cells_ratio,  open("cell_ratio.p", "wb"))
    pickle.dump(cells_result, open("cell_result.p", "wb"))

    # vary number of cells and fix length, moving dx
    cell_list = [8, 16, 32, 64, 128]
    
    dx_ratio = []
    dx_result = []
    
    for cells in cell_list:
        print("\n=== Dx: {} ===".format(4 * np.pi / cells)) 
        for i in range(5):
            
            params = {"runtime": RUNTIME,
                      "cells": [cells],
                      "length": [4 * np.pi]}
            
            # run simulation 
            full_params, delta_ee, delta_ke, ratio = run_sim.two_stream(
                "two-stream_many", 0.5, -0.5, param_dict=params)
   
            delta_ee_list.append(delta_ee)
            delta_ke_list.append(delta_ke) 
            dx_ratio.append(ratio)
            dx_result.append(full_params["dx"][0])
    
    pickle.dump(dx_ratio,  open("dx_ratio.p", "wb"))
    pickle.dump(dx_result, open("dx_result.p", "wb"))

    print("\n ===== NPPC =======")
    print(nppc_ratio)
    print(nppc_result)
    print("\n ===== Cells =======")
    print(cells_ratio)
    print(cells_result)
    print("\n ===== Dx =======")
    print(dx_ratio)
    print(dx_result)

