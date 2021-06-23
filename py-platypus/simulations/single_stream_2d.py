# example of setting up, running, and graphing a single stream instability
# in two dimensions

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim


if __name__ == "__main__":
    
    params = {
        "single_stream": {
            "stream_v": 1, 
            "stream_width": 0.5, 
            "stream_frac": 0.8
        },
        "runtime": 4,
        "print_every": 5,
        "nppc": 25
    }
    
    run_sim.single_stream("single-stream-2d", 2, param_dict=params)
