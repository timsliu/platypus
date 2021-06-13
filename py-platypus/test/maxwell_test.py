# test script demonstrating technique for creating a Maxwellian distribution

import numpy as np
import matplotlib.pyplot as plt


def maxwell_test():
    num_points = 100000
    v = np.zeros(num_points)
    
    # create array of points w/ Maxwellian distribution
    for i in range(num_points):
        r1 = np.random.rand()
        r2 = np.random.rand()
        v[i] = np.sqrt(-np.log(r1)) * np.cos(2 * np.pi * r2)

    # plot the distribution
    bins = np.linspace(-3, 3, 100)
    plt.hist(v, bins=bins)
    plt.show()

if __name__ == "__main__":
    maxwell_test()
