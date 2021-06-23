# example of setting up, running, and graphing landau damping

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim


if __name__ == "__main__":
    run_sim.landau("landau_2d", 2, param_dict={"runtime": 5})
