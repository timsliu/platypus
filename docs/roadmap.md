# Roadmap

This document outlines the features, updates, and extensions planned for 
Platypus. The features are organized into major features and minor features.
Major features are significant expansions of the simulation capabilities, while
the minor features are smaller changes to the codebase. The planned expansions
for major features are further divided into phases, which roughly correspond
to the order they will be implemented.

## Major features

### Phase 1
* 2D py-platypus simulator - create a 2D PIC simulator class in py-platypus.
This simulator will be similar to the 1D simulator with immobile ions and only
electrostatic forces.
* Add examples of Landau damping, two stream instability, and single stream
instability for 2D PIC.
* Documentation explaining how to use py-platypus

### Phase 2
* 1D electrostatic Platypus PIC simulator
* Move the params class so it's usable by both py-platypus and platypus
* Set up system for specifying Platypus simulation in Python and running the
C++ Playtpus PIC
* Basic unit tests
* Translation layer to convert data saved from Platypus simulator to format
usable by the python visualization scripts
* Adapt plotting tools to work for both Platypus and py-platypus

### Phase 3
* Electromagnetic 2D py-platypus PIC
* Create a base 2D PIC class that's used by both the electrostatic and
electromagnetic simulator, or add electromagnetism to the existing PIC

### Phase 4
* 2D electrostatic Platypus PIC simulator

### Phase 5
* CUDA accelerated Platypus

## Minor features
* Make py-platypus graphing functions less annoying to use and configurable
using a vis-params class
* Add a way to plot many subplots that works with varying number of plots
* Create basic continuous integration (CI) tests to run when merging
* Virtual environment to manage installed Python packages
* Optionally mobile ions
* Check for parameter consistency/validity
