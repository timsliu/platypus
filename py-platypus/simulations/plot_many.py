import numpy as np
import matplotlib.pyplot as plt
import pickle


if __name__ == "__main__":
    #nppc_ratio = pickle.load(open("nppc_ratio.p", "rb"))
    #nppc_result = pickle.load(open("nppc_result.p", "rb"))
    #
    #cells_ratio = pickle.load(open("cell_ratio.p", "rb"))
    #cells_result = pickle.load(open("cell_result.p", "rb"))
    #
    dx_ratio = pickle.load(open("dx_ratio.p", "rb"))
    dx_result = pickle.load(open("dx_result.p", "rb"))

    #plt.figure(1)
    #plt.scatter(nppc_result, nppc_ratio)
    #plt.xlabel("nppc")
    #
    #plt.figure(2)
    #plt.scatter(cells_result, cells_ratio)
    #plt.xlabel("cells")
    
    plt.figure(3)
    plt.scatter(dx_result, dx_ratio)
    plt.plot([0, 3.5], [0, -3.5], c="red")
    plt.xlabel("dx")

    plt.show()
