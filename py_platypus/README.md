# PyPlatypus

PyPlatypus is a pure Python3 implementation of several particle in cell
(PIC) plasma simulators. 

PyPlatypus currently includes the following simulators:
* 1-D electrostatic
* 2-D electrostatic
* 3-D electrostatic
* 2-D electromagnetic (in development)


along with example simulations, various tests, and a simple visualization suite.
PyPlatypus is designed to be a proving ground for PIC simulators before they
graduate and are implemented in the main Platypus version. PyPlatypus simulators
are not optimized for performance, but they were intentionally designed
to be easy to understand.


## Table of contents
1. [Directory contents](#directory-contents)
2. [Supported features](#supported-features)
3. [Installation](#installation)
4. [Running your first simulation](#running-your-first-simulation)
6. [Visualization](#visualization)
7. [Troubleshooting](#troubleshooting)

## [Directory contents](#directory-contents)
```
PyPlatypus
|- out: outputs from running simulations
|- simulations: example simulations using PyPlatypus
|- src: actual PyPlatypus Python library and Python egg
  |- py_platypus:
    |- models: PIC simulators
    |- utils: various utility functions used by the simulators
    |- vis: visualization tools
|- test: test scripts for PyPlatypus
  |- unit_tests: unit tests for PyPlatypus
  |- demo_tests: simulations for demonstrating and testing specific features
```

## [Supported features](#supported-features)
PyPlatypus currently includes three PIC simulators, ```pic_1d.py```,
```pic_2d.py```, ```pic_3d.py```, which are 1D, 2D, and 3D simulators respectively.
All three are purely electrostatic, collisionless simulators with mobile electrons and
immobile ions. The boundary conditions are periodic along all axes.

You will also find a ```pic_2d_em.py```, which is a 2D electromagnetic simulator.
This simulator can be used, but there are still some bugs being worked out. Use
with caution! Since the EM PIC is still being developed, the documention below
focuses on the more stable electrostatic simulators.

There are several preset simulations that can be configured for either
PIC simulator. These simulations create an initial velocity or position
distribution that illustrate common plasma physics phenomena:

1. Landau damping - initializes a sinusoidal electron density. This illustrates
the case where an electrostatic wave in a plasma is damped.
2. Two-stream instability - two streams of electrons traveling in opposite
directions. This illustrates the two stream instability, where the electrostatic
energy suddenly increases.
3. Single stream - single stream of particles traveling through a plasma with
a Maxwellian velocity distribution. This illustrates neutral beam injection,
where a beam of particles is shot into a plasma to increasee its temperature.

All simulations have multiple parameters such as the simulation resolution and
run time that can be customized. For a description of how to customize a 
simulation, see the section [Customizable parameters](#customizable-parameters).

## [Installation](#installation)

The installation instructions require the user to have some familiarity with
using the terminal. Note that this information is duplicated from the main
README page for Platypus. To install PyCli, first clone the repository using the
terminal:

```
git clone https://github.com/timsliu/platypus.git
```

Next, navigate to the python source file and install the ```PyPlatypus```
python library:

```
cd platypus/py_platypus/src/
pip3 install -e .
```

Then, you will need to add the root of the platypus directory to your
```~/.bashrc``` file as an environment variable:

```
export PLATYPUS_HOME=`/path/to/platypus`
```

Finally, source your ```~/.bashrc``` for the new environment variable
to take effect:

```
source ~/.bashrc
```

## [Running your first simulation](#running-your-first-simulation)
Simulations should be placed in the ```PLATYPUS_HOME/py_platypus/simulations```
folder. Here you will also find several example simulations. This section 
describes each part of the ```two_stream_1d.py``` simulation.

To run the simulation, simply run the command:

```bash
cd platypus/py_platypus/simulations
python3 two_stream_1d.py
```


The program will generate a new folder in the ```py_platypus/out``` directory
containing output data, graphs, and a json file of the simulation parameters.
Below is a description of what each line of the program does.

At the start of the simulation are several import statements:

```python
import os

import py_platypus as plat
from py_platypus.utils.params import Parameters as Parameters
from py_platypus.vis.plotter import Plotter as Plotter
```

These line import the ```os``` library followed by the ```py_platypus```
library. The lines importing ```Parameters``` and ```Plotter``` allow these
classes to be used without using their full names (```plat.params.Parameters```).


```python
    # override default simulation parameters
    param_dict = {
        "runtime": 40,    # time to run simulation for
        "timestep": 0.04, # size of each timestep
        "save_every": 4,  # save outputs at every fourth step
    }
```
PyPlatypus uses a ```Parameters``` object to set up the simulation. There
are default values for all of the parameters, which can be overriden with
a user defined dictionary. Here, the dictionary is overriding how long the
simulation runs for, the size of each timestep, and how often to save the
simulation outputs. For a full list of the parameters, see 
[Customizable parameters](#customizable-parameters). 

```python
    # setup and run a 1D two stream simulation using the parameters
    # specified in param_dict. This method takes care of initialization
    # and runs the simulation
    plat.run_sim.two_stream("two_stream_1d", 1, param_dict=param_dict)
  
```

Next, the ```two_stream``` function from the ```run_sim``` module is called. 
This function creates an instance of the 1-D PIC simulator class and sets up a two-stream
instability. The arguments to this function are the name of the simulation
(which is used for naming the output directory), and the number of dimensions
for the simulation. These is also an optional argument ```param_dict``` which 
is a dictionary of parameters for customizing the simulation. If this argument
is not passed, then the default configuration for a two stream instability is
used. This single line instantiates and kicks off the simulation.

The full list of parameters used in a simulation is saved to a json file when a
simulation starts running. The following two lines create the path to the
parameters json, and loads them into an instance of the ```Parameters``` class.


```python
    # load the full parameters object from a json
    param_json = os.path.join(PLATYPUS_HOME, "py_platypus/out/two_stream_1d/params.json")
    params = Parameters(1, load_file=param_json)
```

PyPlatypus includes a visualization suite for generating graphs. The next lines
create an instance of the ```Plotter``` class. The ```add_animation``` method
tells the ```Plotter``` instance to generate animations for each quantity that
it is asked to plot.

```python
    # create instance of the Plotter class for generating charts
    plotter = Plotter("two_stream_1d", params)
    plotter.add_animation()   # outputs should be animations
```

Finally, there are four method calls that plot various quantities:

```python
    plotter.plot_velocity()   # plot the velocity distribution
    plotter.plot_phase()      # plot the phase chart (velocity and position)
    plotter.plot_density()    # plot electron density
    plotter.plot_energy()     # plot energy
```

The graphs are saved in the directory ```platypus/py_platypus/out/two_stream_1d/graphs```. 
The output of the simulation, including the density, electric
field, and velocity and particle distribution at each step, is saved as pickle
files. To make custom graphs, the user can directly load the data from the
pickle files.

## [Customizable parameters](#customizable-parameters)
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
    * ```"stream_frac": 1``` - (float 0 to 1) fraction of particles in the stream area that are in the stream
    * ```"stream_width": 0.5``` - (float) width of the stream in each direction; only for the 2D and 3D case.

## [Visualization](#visualization)

The ```Plotter``` class can be used to create several default visualizations of
simulation results. Users can also create their own visualizations by opening
the pickle output files saved in the ```PyPlatypus/out/<simulation_name>/data```
directory. For an example of using the ```Plotter``` class, please see the
section [Running your first simulation](#running-first-sim).

The methods for generating plots are:

* ```plot_electric_field()``` - electric field at several timesteps
* ```plot_energy()``` - line plot of electrostatic and kinetic energy
* ```plot_density()``` - electron density at several timesteps
* ```plot_phase()``` - electron velocity as a function of position
* ```plot_velocity()``` - histogram showing the velocity distribution at several
time steps
* ```plot_all()``` - calls all of the above plotting methods.

Below are several examples of plots generated by the ```Plotter``` class.

Phase plot for the 1D two stream instability. Notice the formation of the
"eye" structure as time progresses.

![Alt text](../images/two_stream_1d_phase.png?raw=true "Phase plot 1d.")

Electron density plot illustrating Landau damping.

![Alt text](../images/landau_2d_density.png?raw=true "Landau 2d density plot.")

## [Troubleshooting](#troubleshooting)

#### TypeError: expected str, bytes, or osPathLike object, not NoneType

This error occurs if you have not set the environment variable ```PLATYPUS_HOME```
in your ```~/.bashrc```. To resolve this error, be sure this environment variable
is set and source the ```~/.bashrc``` file:

```bash
source ~/.bashrc
```
