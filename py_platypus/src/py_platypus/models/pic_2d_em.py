# 2D electromagnetic PIC

import numpy as np

import py_platypus as pla
from py_platypus.models.pic_2d import PIC_2D as PIC_2D
from py_platypus.utils.charge_step import ChargeStep as ChargeStep
from py_platypus.utils.charge_step import ChargeStepDivider as ChargeStepDivider
from py_platypus.utils import math_utils

import matplotlib.pyplot as plt


class PIC_2D_EM(PIC_2D):
    def __init__(self, params):

        super(PIC_2D_EM, self).__init__(params)

        self.bz = np.zeros(self.nodes)    # magnetic field on cell corners
        self.delta_bz = np.zeros(self.nodes)   # value to update magnetic field each  half step
        self.jx = np.zeros([self.cells[0], self.nodes[1]])   # current density in positive x direction
        self.jy = np.zeros([self.nodes[1], self.cells[0]])   # current density in positive y direction

        # vertical cell boundaries holding Ex values
        self.ex_edges = np.zeros([self.cells[0], self.nodes[1]])
        # horizontal cell boundaries holding Ey values
        self.ey_edges = np.zeros([self.nodes[0], self.cells[1]])

        # variables pointing to the current and previous electron positions
        self.electron_x_last = np.zeros(self.n_particles)
        self.electron_y_last = np.zeros(self.n_particles)
       
        # electric and magnetic field interpolated to each particle
        self.e_particle = np.zeros((self.n_particles, self.dimensions))
        self.b_particle = np.zeros(self.n_particles)
        self.charge_divider = ChargeStepDivider(self.dx, self.cells)

    def update_e(self):
        '''calculates the B field update for a half step using Ampere's law:
        curl(B) = u(J + epsilon dE/dt)'''

        # iterate through electric field in ex direction
        for i in range(self.nodes[0]):
            for j in range(self.cells[1]):
                #self.ex_edges[i, j] += 
                pass

        # iterate through electric field in ey_direction

    def calc_B_update(self):
        '''calculates the B field update for a half step using Faraday's law: 
        curl(E) = -dB/dt'''
        delta_x = self.dx[1]
        delta_y = self.dx[0]

        # TODO figure out the constants
        # iterate through cells and calculate the B field update 
        for i in range(self.bz.shape[0]):
            for j in range(self.bz.shape[1]):
                ex0 = math_utils.wrap_idx_2d(self.ex_edges, i - 1, j)
                ex1 = math_utils.wrap_idx_2d(self.ex_edges, i, j) 
                ey0 = math_utils.wrap_idx_2d(self.ey_edges, i, j - 1) 
                ey1 = math_utils.wrap_idx_2d(self.ey_edges, i, j) 

                dey = ey1 - ey0 
                dex = ex1 - ex0 
                
                # curl in 2D is dex/dy - dey/dx
                self.delta_bz[i][j] = (dex/delta_y - dey/delta_x) * self.dt 
 
        # divide update in half to calculate B update at each half step
        self.delta_bz /= 2
        return

    def update_B_half(self):
        '''update magnetic field B using the values calculated by 
        calc_B_update'''

        self.bz += self.delta_bz

        return

    def update_j(self):
        '''calculate the current density J using the current and last
        x position'''

        # iterate through all particles and calculate which cell boundaries
        # are crossed 
        for i in range(self.n_particles):
            
            # initial and final particle position 
            x0, y0 = self.electron_x_last[i], self.electron_y_last[i]
            x1, y1 = self.electron_x[i], self.electron_y[i]

            substeps = self.charge_divider(x0, y0, x1, y1)

            # iterate over the substeps
            for step in substeps:
                execute_four_bound(step)

        return

    def execute_four_bound(self, cs):
        ''' 
        Calculates the charge density update from the motion of a charge
        inputs: ChargeStep object representing the motion of a charge 
        '''
        dx = self.dx[1]
        dy = self.dx[0]

        # get list of horizontal and vertical boundaries touched
        # method guarantees the lower indexed boundary is listed fist
        hori, vert = self.charge_divider.boundaries_touched(cs.x0, cs.y0)
       
        # coordinates of the grid point at the center of crossed boundaries
        local_origin = dx * vert[0][1], dy * hori[0][0]

        # starting particle coordinates relative to local origin
        local_x0 = cs.x0 - local_origin[0]
        local_y0 = cs.y0 - local_origin[1]

        # (Villasenor and Buneman 1991)
        # The signs for jy do not match b/c Villasenor and Buneman
        # defined the positive y direction opposite to what it used here
        jx1 = cs.dx * (0.5 - local_y0 - 0.5 * cs.dy)
        jx2 = cs.dx * (0.5 + local_y0 + 0.5 * cs.dy)
        jy1 = - cs.dy * (0.5 - local_x0 - 0.5 * cs.dx)
        jy2 = - cs.dy * (0.5 + local_x0 + 0.5 * cs.dx)

        # update the current densities
        self.jx[tuple(hori[0])] += jx2
        self.jx[tuple(hori[1])] += jx1
        self.jy[tuple(vert[0])] += jy1
        self.jy[tuple(vert[1])] += jy2

        return

    def init_B(self):
        '''set up initial conditions for the magnetic field B'''
        # set initial magnetic field to be zero for now
        self.bz = np.zeros(self.bz.shape)

        return

    def init_E(self):
        '''set up initial conditions for the E field'''
        
        # calculate the initial electric field from Poisson's law
        super().update_ni()
        super().update_ne()
        super().update_rho()
        super().update_phi()
        super().update_e()   # electric field on the nodes (corners)

        # convert electric field on the nodes to the edges of the Yee mesh
        # Ex fields
        for i in range(self.ex_edges.shape[0]):
            for j in range(self.ex_edges.shape[1]):
                self.ex_edges[i][j] = np.mean([
                    self.ex[i][j], 
                    self.ex[i + 1][j]])

        # Ey fields
        for i in range(self.ey_edges.shape[0]):
            for j in range(self.ey_edges.shape[1]):
                self.ey_edges[i][j] = np.mean([
                    self.ey[i][j],
                    self.ey[i][j + 1]])
        
        return

    def interpolate(self, x_n, y_n, corners, values):
        '''interpolates the value of a field property to the position of
        a particle
        inputs - x_n - x coordinate of the particle
                 y_n - y coordinate of the particle
                 corners - (4, array) coords of the corners [x0, x1, y0, y1]
                 values - (4, array) values of the field at the four corners
                          [upper left, upper right, lower left, lower right]

        returns - interpolated value at the particle'''

        x0, x1, y0, y1 = corners

        area_upper_left  = math_utils.points_to_area((x_n, y_n), (x0, y0))
        area_upper_right = math_utils.points_to_area((x_n, y_n), (x1, y0))
        area_lower_left  = math_utils.points_to_area((x_n, y_n), (x0, y1))
        area_lower_right = math_utils.points_to_area((x_n, y_n), (x1, y1))

        # total area of a cell
        total_area = self.dx[0] * self.dx[1]

        # calculate weight to be distributed to each quadrant
        # indexing is row, col (y, x)
        weight_00 = area_lower_right/total_area 
        weight_01 = area_lower_left/total_area
        weight_10 = area_upper_right/total_area 
        weight_11 = area_upper_left/total_area

        # list of weights and list of electric fields
        weights = np.array([weight_00, weight_01, weight_10, weight_11])

        return np.sum(list(map(lambda a, b: a * b, weights, values)), axis=0)

    def interpolate_e(self):
        '''interpolate the E field at each particle'''
      
        # iterate through the particles
        for i in range(self.n_particles): 
            x_n = self.electron_x[i]
            y_n = self.electron_y[i]

            # ==== interpolate Ex field ====
            # Ex field lies on vertical edges, effectively on x nodes and on
            # y cells
            node_l, node_r, x0, x1 = self.node_neighbors(x_n, "x")
            cell_u, cell_d, y0, y1 = self.cell_neighbors(y_n, "y")

            # look up Ex field at the four corners around particle
            # indexing is row, col (y, x)
            e_00 = math_utils.wrap_idx_2d(self.ex_edges, cell_u, node_l) 
            e_01 = math_utils.wrap_idx_2d(self.ex_edges, cell_u, node_r)
            e_10 = math_utils.wrap_idx_2d(self.ex_edges, cell_d, node_l)
            e_11 = math_utils.wrap_idx_2d(self.ex_edges, cell_d, node_r)
            
            self.e_particle[i][0] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [e_00, e_01, e_10, e_11])

            # ==== interpolate Ey field ====
            # Ey field lies on horizontal edges, effectively on y nodes and on
            # x cells
            cell_l, cell_r, x0, x1 = self.cell_neighbors(x_n, "x")
            node_u, node_d, y0, y1 = self.node_neighbors(y_n, "y")

            # look up Ey field at the four corners
            # indexing is row, col (y, x)
            e_00 = math_utils.wrap_idx_2d(self.ey_edges, node_u, cell_l) 
            e_01 = math_utils.wrap_idx_2d(self.ey_edges, node_u, cell_r)
            e_10 = math_utils.wrap_idx_2d(self.ey_edges, node_d, cell_l)
            e_11 = math_utils.wrap_idx_2d(self.ey_edges, node_d, cell_r)
            
            self.e_particle[i][1] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [e_00, e_01, e_10, e_11])
            
        return

    def interpolate_b(self):
        '''interpolate the B field at each particle'''
      
        for i in range(self.n_particles):

            x_n = self.electron_x[i]
            y_n = self.electron_y[i]

            # find indices of four nodes neighboring the particle
            node_u, node_d, y0, y1 = self.node_neighbors(y_n, "y")
            node_l, node_r, x0, x1 = self.node_neighbors(x_n, "x")
          
            # B field at each corner
            # indexing is row, col (y, x)
            b_00 = self.bz[node_u][node_l] 
            b_01 = self.bz[node_u][node_r] 
            b_10 = self.bz[node_d][node_l] 
            b_11 = self.bz[node_d][node_r] 

            # indices of neighboring cells
            self.b_particle[i] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [b_00, b_01, b_10, b_11])
            
        return

    def update_v(self):
        '''updates the velocity according to Lorentz's law 
        dv/dt = q(e + v x B) and the Boris method. Overrides the method
        in PIC_2D base class'''

        # TODO need to work in constants
        for i in range(self.n_particles):
            # 3D vector of electric field at the particle 
            e_i = np.concatenate([self.e_particle[i], [0]])
            v_i = np.array([self.electron_vx[i], self.electron_vy[i]]) 

            # Boris method for updating position
            q_prime = self.dt/2
            h = [0, 0, q_prime * self.b_particle[i]]
            s = 2 * h/(1 + np.linalg.norm(h) ** 2)
            u = v_i + q_prime * e_i
            u_prime = u + np.cross((u + np.cross(u, h)), s)
           
            # calculate updated velocity (3D vector) and update
            v = u_prime + q_prime * e_i
            self.electron_vx[i] = v[0]
            self.electron_vy[i] = v[1]

        return

    def step(self):
        '''run the simulation for a single step, updating all parameters;
        methods for saving outputs must be called separately; overrides
        step method for PIC_2D class'''

        # initialize B and E fields on the first step
        if self.step == 0:
            self.init_B()
            self.init_E()

        self.calc_B_update()
        self.update_B_half()
        self.interpolate_e()
        self.interpolate_b()
        self.update_v()
        self.update_x()
        self.update_J()
        self.update_B_half()
        self.update_e()
        self.step += 1

        return
