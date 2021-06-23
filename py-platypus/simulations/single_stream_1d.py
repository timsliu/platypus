# example of setting up, running, and graphing a single stream instability
# in one dimension

import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))

import run_sim


if __name__ == "__main__":
    run_sim.single_stream("single-stream-1d", 1)