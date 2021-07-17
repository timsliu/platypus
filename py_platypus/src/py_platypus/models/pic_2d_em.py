# 2D electromagnetic PIC

import numpy as np

import py_platypus as pla
from py_platypus.models.pic_2d import PIC_2D as PIC_2D


class PIC_2D_EM(PIC_2D):
    def __init__(self, params):

        super().__init__(params)

        self.Bz = np.zeros(self.cells)    # magnetic field on cell centers
        self.dBz = np.zeros(self.cells)   # value to update magnetic field each  half step
        self.Jx = np.zeros([self.cells[0], self.nodes[1]])   # current density in x direction
        self.Jy = np.zeros([self.nodes[1], self.cells[0]])   # current density in y direction

        self.Ex_edges = np.zeros([self.cells[0], self.nodes[1]])
        self.Ey_edges = np.zeros([self.nodes[0], self.cellss[1]])

        # variables pointing to the current and previous electron positions
        self.electron_x_last = np.zeros(self.n_particles)
        self.electron_y_last = np.zeros(self.n_particles)

    def update_e(self):
        '''update electric field E using Faraday's law: curl(E) = -dB/dt'''

    def calc_B_update(self):
        '''calculates the B field update for a half step using Ampere's law:
        curl(B) = u(J + epsilon dE/dt); the actual B field is not updated'''
    
    def update_B_half(self):
        '''update magnetic field B using the values calculated by 
        calc_B_update'''

    def update_J(self):
        '''calculate the current density J using the current and last
        x position'''

    def init_B(self):
        '''set up initial conditions for the magnetic field B'''

        # set initial magnetic field to be zero for now
        self.Bz = np.zeros(self.cells)

        return

    def init_E(self):
        '''set up initial conditions for the E field'''
        
        # calculate the initial electric field from Poisson's law
        super().update_ni()
        super().update_ne()
        super().update_rho()
        super().update_phi()
        super().update_e()   # electric field one the nodes (corners)

        # convert electric field on the nodes to the edges of the Yee mesh

    def step(self):
        '''run the simulation for a single step, updating all parameters;
        methods for saving outputs must be called separately; overrides
        step method for PIC_2D class'''

        if self.step == 0:
            self.init_B()
            self.init_E()

        self.calc_B_update()
        self.update_B_half()
        self.update_v()
        self.update_x()
        self.update_J()
        self.update_B_half()
        self.update_e()
        self.step += 1

        return
