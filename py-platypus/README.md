# py-platypus

py-platypus is a pure Python3 implementation of several particle in cell
(PIC) plasma simulators. Currently, py-platypus includes 1D and 2D electrostatic
simulators and a simple visualization suite. py-platypus is designed to be 
a proving ground for PIC simulators before they ``graduate`` and are implemented
in the main Platypus version. The performance of py-platypus is not optimized,
and is written to be easy for new users to understand.

## Table of contents
1. [Supported features](#supported-features)
2. [Installation](#installation)
3. [Running your first simulation](#running-first-sim)
4. [Directory contents](#directory-contents)
5. [Troubleshooting](#troubleshooting)


## [Supported features](#supported-features)
py-platypus currently includes two PIC simulators, ```pic_1d.py``` and
```pic_2d.py```, which are 1D and 2D simulators respectively. Both are
purely electrostatic, collisionless simulators with mobile electrons and
immobile ions. The boundary conditions are periodic along all axes.

There are several preset simulations that can be configured for either
PIC simulator. These simulations create an initial velocity or position
distribution that illustrate common plasma physics phenomena.

1. Landau damping - sinusoidal electron density. This illustrates the case
where an electrostatic wave in a plasma is damped.
2. Two-stream instability - two streams of electrons traveling in opposite
directions. This illustrates the two stream instability, where the electrostatic
energy suddenly increases.
3. Single stream - single stream of particles traveling through a plasma with
a Maxwellian velocity distribution. This illustrates neutral beam injection,
where a beam of particles is shot into a plasma to increasee its temperature.

All simulations have multiple parameters such as the simulation resolution and
run time that can be customized. For a description of how to customize a 
simulation, see the section [Customizable parameters](#parameters)

## [Installation](#installation)

## [Running your first simulation](#running-first-sim)
Simulations should be placed in the ```PLATYPUS_HOME/py-platypus/simulations```
folder. Here you will also find several example simulations. This section 
describes each part of the ```two_stream_1d.py``` simulation and how to write
your own simulation.

To run the simulation, simply run the command:

```bash
cd PLATYPUS_HOME/py-platypus/simulations
python3 two_stream_1d.py
```

The program will great a new folder in the ```py-platypus/out``` directory
containing output data, graphs, and a json file of the simulation parameters.
Below is a description of what each line of the program does.

At the start of the simulation are several import statements:

```python
import sys
import os

PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus"))
sys.path.append(os.path.join(PLATYPUS_HOME, "py-platypus/vis"))

import run_sim
from plotter import Plotter
from params import Parameters 

```

These lines tell the program where Platypus is installed and imports the
needed data structures and functions. These lines should be included in all 
simulations.

```python
if __name__ == "__main__":
    run_sim.two_stream("two_stream_1d", 1, param_dict={"runtime": 40})
  
```

When the program is run, the ```two_stream``` method is called. This function
creates an instance of the PIC simulator class and sets up a two-stream
instability. The arguments to this function are the name of the simulation
(which is used for naming the output directory), and the number of dimensions
for the simulation. These is also an optional argument ```param_dict``` which 
is a dictionary of parameters for customizing the simulation. If this argument
is not passed, then the default configuration for a two stream instability is
used. See the section [Customizable parameters](#parameters) for info on
how to customize the simulation. This single line instantiates and kicks off
the simulation.

Simulations are set up and configured by a ```Parameters``` class. This
class is created by the ```two_stream``` function in ```run_sim```. The full
list of parameters used in a simulation is saved to a json file when a
simulation starts running. The following two lines create the path to the
parameters json, and loads them into an instance of the ```Parameters``` class.


```python
    param_json = os.path.join(PLATYPUS_HOME, "py-platypus/out/two_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
```

py-platypus includes a visualization suite for generating graphs. The next lines
instantiate an instance of the ```Plotter``` class and generates several default
graphs.

```python
plotter = Plotter("two_stream_1d", params)
plotter.plot_all() 
```

These graphs are saved in the directory ```PLATYPUS_HOME/py-platypus/out/two_stream_1d/graphs```. The output of the simulation, including the density, electric
field, and velocity and particle distribution at each step, is saved as pickle
files. To make custom graphs, the user can directly load the data from the
pickle files.

## [Customizable parameters](#parameters)
The ```Parameters``` class specifies the parameters of a simulation. An instance
of the ```Parameters``` class is created by the functions in ```run_sim``` and
can be modified by the ```param_dict``` argument. The values of the class
are the defaults found in ```params.py``` until modified by the user. Below
is a full description of the parameters that can be modified along with their
default values. There are also several derived values that are not set by the
user but calculated based on other values.
            
* ```"name": "default_simulation"``` - (string) name of the simulation, used
for creating the output directory
* ```"seed": 0``` - (int) random seed used for generating initial x position and
velocity distribution
* ```"version": "1.0"``` - (string) version number of the simulation, allowing
for multiple versions of the same simulation to be saved.
* ```"length": [2 * np.pi for x in range(dims)]``` - (array) list with the
size of the simulation space along each axes.
* ```"cells": [32 for x in range(dims)]``` - (array) number of cells along each
axes
* ```"timestep": 0.04``` - (float) size of each timestep of the simulation
* ```"runtime": 30``` - (float) total amount of time to run the simulation for
* ```"dimensions": dims``` - (int) number of dimensions for the simulation
* ```"nppc": 100``` - (int) number of particles per cell
* ```"dx": None``` - (float array) - size of the cells along each dimension (derived)
* ```"steps": None``` - (float) total number of simulation steps (derived)  
* ```"n_particles"``` - (int) total number of particles (derived)
* ```"print_every": 20``` - (int) time steps between printing current step
* ```"save_every": 100``` - (int) interval between saving data
* ```"single_stream"``` - (dict) dictionary of parameters for the single stream
instability; all values MUST be set by users or left at default
    * ```"stream_v": 1``` - (int) velocity of particles in the stream 
    * ```"stream_frac": 0.5``` - (float 0 to 1) fraction of particles in the stream set to the stream velocity 
    * ```"stream_width": 0.5``` - (float) width of the stream
* ```"landau"``` - (dict) dictionary of parameters for Landau damping; all values
MUST be set by users or left at default
    * ```amplitude": 0.5``` - (float 0 to 1) amplitude of the density perturbation
    * ```mode": 1``` - (int) number of peaks in the density perturbation
* "two_stream" - (dict) dictionary of parameters for the two stream instability;
MUST be set by users or left at default
    * ```"vpos": 0.5``` - (float) velocity of positive stream
    * ```"vneg": -0.5``` - (float) velocity of negative stream
    * ```"stream_frac": 1 - (float 0 to 1) fraction of particles in the stream area that are in the stream
    * ```"stream_width": 0.5``` - (float) width of the stream in each direction; only for the 2D and 3D case.

## [Visualization](#visualization)

## [Troubleshooting](#troubleshooting)
