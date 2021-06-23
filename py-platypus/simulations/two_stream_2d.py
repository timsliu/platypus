# example of setting up, running, and graphing a two stream instability in
# two dimensions

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim

if __name__ == "__main__":

    params = {
        "print_every": 5,
        "nppc": 40,
        "runtime": 5
    }
    run_sim.two_stream("two-stream-2d", 2, param_dict = params)
