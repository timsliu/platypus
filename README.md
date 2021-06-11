# Platypus
Platypus is an open source library of particle in cell (PIC) plasma simulators.
Platypus will ultimately include 1-D, 2-D, and 3-D simulators that account
for both electrostatic and electromagnetic effects. The respository also
includes scripts for creating graphs and visualizations of the simulation
results. 

Platypus includes two sets of simulations. The optimized, high performance
version of Platypus is written in C++14 and located in the ```src``` directory.
The instructions below detail how to set up and run a simulation with 
Platypus. Thhe second set of simulations is py-platypus, a pure Python3
implementation that is lower performance. The py-platypus simulators are a 
proving ground for the algorithms used in the PIC simulation. These simulators
can be thought of as ```Platypus Light``` and are also designed to be easier
to understand. New features will first be added to py-platypus and tested 
before they ``graduate`` and are implemented in C++.

The following sections describe how to install, run simulations, and create
visualizations with Platypus.

## Table of contents
1. [Environment](#environment)
2. [Installation](#installation)
3. [Running your first simulation](#running-first-model)
4. [Folder organization](#folder-organization)
5. [Contributors and Contact](#contributors-and-contact)
6. [License](#license)

## [Environment](#environment)

Platypus has the following system requirements:

Python:
```


```

Platypus was developed on Unix based operating systems and is not guaranteed
to run on Windows.


## [Installation](#installation)

## [Running your first simulation](#running-first-simu)

## [Folder organization](#folder-organization)
This section describes the organization of the Platypus repository.

### bin
Output directory for binaries and object files.

### py-platypus
Py-platypus is a Python based collection of PIC simulators. They are a
more lightweight, low performance version of the Platypus PIC simulators.
Py-platypus is primarily meant for developing and testing the algorithms used
in the PIC simulators, and are also meant to be easier to understand than the
optimized Platypus simulators.


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

## [License](#license)
Platypus is distributed under the MIT open source license. For the full
license contents, please see `LICENSE`.

