# Platypus
Platypus is an open source library of particle-in-cell (PIC) plasma simulators.
These simulators are designed to be both very easy to use through a simple
Python API and easy to understand, with well documented and cleanly written
code. Platypus is intended as an educational tool, and is not meant to be
a replacement for more sophisticated simulators.

Platypus currently includes the following simulators:
* 1-D electrostatic
* 2-D electrostatic
* 3-D electrostatic
* 2-D electromagnetic (in development)

The respository also includes scripts for creating graphs, videos, and other
visulatizations of the simulation results. There are also numerous example
simulations that users can copy and modify.

Platypus will ultimately include two sets of simulations. The first set of
simulations is PyPlatypus, a pure Python3 implementation that is 
designed to be easy to understand but not optimized for performance. All of the
simulations currently implemented and under development are PyPlatypus
simulations. An optimized,  high performance version of Platypus written in 
C++14 will be implemented at a later date. The ultimate goal is to have a
CUDA GPU implementation of Platypus. Currently, basically everything interesting
is located in the PyPlatypus directory.

New features will first be added to PyPlatypus and tested 
before they ``graduate`` and are implemented in C++. 

The following sections give a quick overview on how to install, 
run simulations, and create visualizations with PyPlatypus. Check out
the dedicated [README for PyPlatypus](https://github.com/timsliu/platypus/tree/main/PyPlatypus)
for a lot more details. Since only the Python3 PyPlatypus version is currently,
available, this README has some information duplicated from the dedicated
PyPlatypus page.

![Alt text](images/two_stream_1d_trunc.png?raw=true "Evolution of a two
stream instability in 1D")


## Table of contents
1. [Environment](#environment)
2. [Installation](#installation)
3. [Running your first simulation](#running-your-first-simulation)
4. [Repository contents](#repository-contents)
5. [Contributors and Contact](#contributors-and-contact)
6. [Note on Naming](#note-on-naming)
7. [License](#license)

## [Environment](#environment)

PyPlatypus has the following system requirements:

```
python3.6+
numpy
matplotlib
scipy
ffmpeg
```

Platypus was developed on macOS and is not guaranteed to run on Windows. It
should run fine on any Unix based OS, but this has yet to be tested.


## [Installation](#installation)

### PyPlatypus
The installation instructions require the user to have some familiarity with
using the terminal.To install PyCli, first clone the repository using the
terminal:

```
git clone https://github.com/timsliu/platypus.git
```

Next, navigate to the python source file and install the PyPlatypus
python library:

```
cd platypus/PyPlatypus/src/
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
PyPlatypus includes several example simulations that you can run without
writing any lines of code. This section gives a quick overview of how to run
a simulation, and be sure to checkout the dedicated [PyPlatypus README page](https://github.com/timsliu/platypus/tree/main/PyPlatypus#running-your-first-simulation) for more details.
To view the example simulations, navigate to the
```simulations``` sub directory of PyPlatypus:

```
cd platypus/PyPlatypus/simulations
ls
```

This will list the various simulations that are already included in the
repository. Each simulation is simply a Python file that can be executed like
any normal python file:

```
python3 two_stream_1d.py
```

When the simulation finishes, navigate to the ```out``` directory and find 
the folder with the same name as the simulation:

```
cd platypus/py_platypus/out/two_stream_1d
```

Here you will find a ```graphs``` directory, a ```data``` directory, and a file
called ```params.json```. Open the ```graphs``` directory to view various
visualizations gnerated by PyPlatypus.

## [Repository contents](#repository-contents)
This section describes the organization of the Platypus repository.

### bin
Output directory for binaries and object files.

### py-platypus
Py-platypus is a Python based collection of PIC simulators. They are a
more lightweight, low performance version of the Platypus PIC simulators.
Py-platypus is primarily meant for developing and testing the algorithms used
in the PIC simulators, and are also meant to be easier to understand than the
optimized Platypus simulators. PyPlatypus has the following organization:

```
PyPlatypus
|- out: outputs from running simulations
|- simulations: example simulations using PyPlatypus
|- src: actual PyPlatypus Python library and Python egg
  |- PyPlatypus:
    |- models: PIC simulators
    |- utils: various utility functions used by the simulators
    |- vis: visualization tools
|- test: test scripts for PyPlatypus
  |- unit_tests: unit tets for PyPlatypus
  |- demo_tests: simulations for demonstrating and testing specific features
```


### test
Folder containing the test harness and test scripts for Platypus.


### src
C++ source for Platypus.


### pysrc
Python scripts for parsing simulation arguments, checking that the passed
parameters are within bounds, and running the C++ Platypus binary.


### scripts
Python scripts for data analysis and creating visualizations.


## [Contributors and contact](#contributors-and-contact)
Platypus is an open source project created and maintained by Timothy Liu.
You can contact him at [timsliu.org/contact](https://timsliu.org/contact/).

## [Note on Naming](#note-on-naming)
The official name of the project is Platypus (capital P) and the name of the
Python implemention is PyPlatypus (capital P, camel case style). However,
the Python code in this repository follows the convention that variables are
snake_case and class names are CamelCase. The camel case "PyPlatypus" spelling
is found only in the documentation, while folders and code use the snake case
spelling "py_platypus".

## [License](#license)
Platypus is distributed under the MIT open source license. For the full
license contents, please see `LICENSE`.

