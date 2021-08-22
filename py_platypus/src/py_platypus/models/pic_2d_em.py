# 2D electromagnetic PIC

import numpy as np

import py_platypus as pla
from py_platypus.models.pic_2d import PIC_2D as PIC_2D

import matplotlib.pyplot as plt

class charge_step:
    """
    Class representing a charge moving
    """
    def __init__(self, x0, y0, dx, dy):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x0 + dx
        self.y1 = y0 + dy
        self.dx = dx
        self.dy = dy


class PIC_2D_EM(PIC_2D):
    def __init__(self, params):

        super().__init__(params)

        self.Bz = np.zeros(self.cells)    # magnetic field on cell centers
        self.delta_Bz = np.zeros(self.cells)   # value to update magnetic field each  half step
        self.jx = np.zeros([self.cells[0], self.nodes[1]])   # current density in x direction
        self.jy = np.zeros([self.nodes[1], self.cells[0]])   # current density in y direction

        self.Ex_edges = np.zeros([self.nodes[0], self.cells[1]])
        self.Ey_edges = np.zeros([self.cells[0], self.nodes[1]])

        # variables pointing to the current and previous electron positions
        self.electron_x_last = np.zeros(self.n_particles)
        self.electron_y_last = np.zeros(self.n_particles)
       
        # electric and magnetic field interpolated to each particle
        self.e_particle = np.zeros((self.n_particles, self.dimensions))
        self.b_particle = np.zeros(self.n_particles)

    def update_e(self):
        '''calculates the B field update for a half step using Ampere's law:
        curl(B) = u(J + epsilon dE/dt)'''


    

    def calc_B_update(self):
        '''calculates the B field update for a half step using Faraday's law: 
        curl(E) = -dB/dt'''

        # TODO figure out the constants
        # iterate through cells and calculate the B field at the center of 
        # cells
        for i in range(self.cells[0]):
            for j in range(self.cells[1]):
                dy = self.Ey_edges[i][j + 1] - self.Ey_edges[i][j]
                dx = self.Ex_edges[i + 1][j] - self.Ex_edges[i][j]
                delta_x = self.dx[1]
                delta_y = self.dx[0]
                
                # curl in 2D is dy/dx - dx/dy
                self.delta_Bz[i][j] = -(dy/delta_x - dx/delta_y) * self.dt 
 
        # divide update in half to calculate B update at each half step
        self.delta_Bz /= 2
        return

    def update_B_half(self):
        '''update magnetic field B using the values calculated by 
        calc_B_update'''

        self.Bz += self.delta_Bz

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
            # change in position 
            delta_x = x1 - x0
            delta_y = y1 - y0
            # grid size in x and y directions 
            dx = self.dx[1]
            dy = self.dy[0]

            # list of boundaries touched at the starting and ending positions
            hori_0, vert_0 = boundaries_touched(x0, y0)
            hori_1, vert_1 = boundaries_touched(x1, y1)

            # list of unique boundaries touched
            hori_all = list(set(hori_0, hori_1))
            vert_all = list(set(vert_0, vert_1))

            # count unique boundaries touched at start and stop positions
            boundaries_touched = len(hori_all) + len(vert_all)
            
            if boundaries_touched == 4:
                # get the list of four boundary crossing motions
                substeps = get_submotions_four(x0, y0, x1, y1)
            elif boundaries_touched == 7:
                # get the list of four boundary crossing motions
                substeps = get_submotions_seven(x0, y0, x1, y1)
            elif boundaries_touched == 8:
                # touching 8 unique boundaries at start and stop positions
                # corresponds to crossing 10 boundaries
                substeps = get_submotions_ten(x0, y0, x1, y1)
            else:
                print("Warning! {} boundaries touched, expected 4, 7, or 8")

            # iterate over the substeps
            for step in substeps:
                execute_four_bound(step)

    def execute_four_bound(self, step):
        '''
        calculate the charge density contribution of add
        '''


    def get_submotions_four(self, x0, y0, x1, y1):
        '''
        returns a list of charge steps for the four boundary case
        '''
        delta_x = x1 - x0
        delta_y = y1 - y0

        return [charge_step(x0, y0, delta_x, delta_y)]
    
    def get_submotions_seven(self, x0, y0, x1, y1):
        '''
        returns a list of charge steps for the seven boundary case
        '''
        # determine the direction leading to two steps

        delta_x = x1 - x0
        delta_y = y1 - y0

        # grid size in x and y directions 
        dx = self.dx[1]
        dy = self.dy[0]

        if abs(dx) > abs(dy):
            # four vertical boundaries, three horizontal boundaries
            # sub motions dependent on x
            delta_x0 = abs((x0 % dx) - dx/2) * np.sign(delta_x)
            delta_y0 = delta_x0 * (delta_y/delta_x) 
        
        else:
            # four horizontal boundaries, three vertical boundaries
            # sub motions depdendent on y
            delta_y0 = abs((y0 % dy) - dy/2) * np.sign(delta_y)
            delta_x0 = delta_y1 * (delta_x/delta_y) 

        delta_x1 = delta_x - delta_x0
        delta_y1 = delta_y - delta_y0

        c0 = charge_step(x0, y0, delta_x0, delta_y0)
        c1 = charge_step(x0 + delta_x0, y0 + delta_y0, delta_x1, delta_y1)

        return [c0, c1]
    
    def get_submotions_ten(self, x0, y0, x1, y1):
        '''
        returns a list of charge steps for the ten boundary case
        '''
        delta_x = x1 - x0
        delta_y = y1 - y0

        # grid size in x and y directions 
        dx = self.dx[1]
        dy = self.dy[0]

        # fraction of the move to execute to get to the midpoint of a cell
        # along either dimension
        x_frac = dx/2 - (x0 % dx)/delta_x
        y_frac = dy/2 - (y0 % dy)/delta_y

        # shorter to get to the x midpoint - movement to x midpoint is first
        if x_frac < y_frac:
            delta_x0 = delta_x * x_frac
            delta_y0 = delta_y * x_frac

            delta_x1 = delta_x * y_frac - delta_x0
            delta_y1 = delta_y * y_frac - delta_y0

        # shorter to get to the y midpoint - movement to y midpoint is first
        else:
            delta_x0 = delta_x * y_frac
            delta_y0 = delta_y * y_frac
            
            delta_x1 = delta_x * x_frac - delta_x0
            delta_y1 = delta_y * x_frac - delta_y0

        delta_x2 = delta_x - delta_x1 - delta_x0 
        delta_y2 = delta_y - delta_y1 - delta_y0 

        c0 = charge_step(x0, y0, delta_x0, delta_y0)
        c1 = charge_step(c0.x1, c0.y1, delta_x1, delta_y1)
        c2 = charge_step(c2.x1, c0.y1, delta_x2, delta_y2)

        return [c0, c1, c2]

    def execute_four_bound(self, charge_step):
        ''' 
        Calculates the charge density update from the motion of a charge
        inputs: charge_step object representing the motion of a charge 
        '''


        return

    def boundaries_touched(self, x0, y0):
        '''
        returns a list of boundaries touched
        inputs: x, y - coordinates of the particle
        returns: horizontal - list of indices of horizontal edges touched
                 vertical - list of indices of vertical edges touched
        '''
        dx = self.dx[1]
        dy = self.dy[0]
        col_vert_idx = np.ceil((x - dx/2)/dx) # column index of the vertical boundary
        row_hori_idx = np.ceil((y - dy/2)/dy) # row index of the horizontal boundary

        horizontal = [[row_hori_indx, col_vert_idx - 1], [row_hori_idx, col_vert_idx]]
        vertical   = [[row_hori_indx - 1, col_vert_idx], [row_hori_idx, col_vert_idx]]

        horizontal, vertical


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
        # Ex fields
        for i in range(self.nodes[0]):
            for j in range(self.cells[1]):
                self.Ex_edges[i][j] = np.mean([
                    self.ex[i][j], 
                    self.ex[(i - 1) % self.nodes[0]][j]])

        # Ey fields
        for i in range(self.cells[0]):
            for j in range(self.nodes[1]):
                self.Ey_edges[i][j] = np.mean([
                    self.ey[i][j], 
                    self.ey[i][(j - 1) % self.nodes[1]]])
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

        area_upper_left  = pla.utils.points_to_area((x_n, y_n), (x0, y0))
        area_upper_right = pla.utils.points_to_area((x_n, y_n), (x1, y0))
        area_lower_left  = pla.utils.points_to_area((x_n, y_n), (x0, y1))
        area_lower_right = pla.utils.points_to_area((x_n, y_n), (x1, y1))

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

            # left right indices of neighboring vertical edges
            edge_left  = int(np.floor(x_n/self.dx[1]))
            edge_right = int(np.ceil(x_n/self.dx[1]))
            
            # up down indices of neighboring horizontal edges
            edge_up  = int(np.floor(y_n/self.dx[0]))
            edge_down = int(np.ceil(y_n/self.dx[0]))
            
            
            # ==== interpolate Ex field ====
            cell_x = edge_left
            
            # column indices of the horizontal edges to use 
            edge_hleft, edge_hright, x0, x1 = self.cell_neighbors(x_n, 0, cell_x)
            
            y0 = edge_up * self.dx[0]
            y1 = edge_down * self.dx[0]

            # wrap the indices of the edges
            edge_hleft = edge_hleft % self.Ex_edges.shape[1] 
            edge_hright = edge_hright % self.Ex_edges.shape[1] 

            # look up Ex field at the four corners
            e_00 = self.Ex_edges[edge_up][edge_hleft] 
            e_01 = self.Ex_edges[edge_up][edge_hright]
            e_10 = self.Ex_edges[edge_down][edge_hleft]
            e_11 = self.Ex_edges[edge_down][edge_hright]
            
            self.e_particle[i][0] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [e_00, e_01, e_10, e_11])

            # ==== interpolate Ey field ====
            cell_y = edge_up
            edge_vup, edge_vdown, y0, y1 = self.cell_neighbors(y_n, 1, cell_y) 
            
            x0 = edge_left * self.dx[1]
            x1 = edge_right * self.dx[1] 
            
            # wrap the indices of the edges
            edge_vup = edge_vup % self.Ey_edges.shape[0]
            edge_vdown = edge_vdown % self.Ey_edges.shape[0]
           
            # look up Ey field at the four corners
            e_00 = self.Ey_edges[edge_vup][edge_left] 
            e_01 = self.Ey_edges[edge_vup][edge_right]
            e_10 = self.Ey_edges[edge_vdown][edge_left]
            e_11 = self.Ey_edges[edge_vdown][edge_right]
            
            self.e_particle[i][1] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [e_00, e_01, e_10, e_11])

        return

    def interpolate_b(self):
        '''interpolate the B field at each particle'''
       
        for i in range(self.n_particles):

            x_n = self.electron_x[i]
            y_n = self.electron_y[i]
            #print(x_n, y_n)

            # find the cell the particle is in
            cell_y = int(np.floor(y_n/self.dx[0]))
            cell_x = int(np.floor(x_n/self.dx[1]))
            # print(cell_x, cell_y)

            # find four cell centers neighboring the particle
            cell_u, cell_d, y0, y1 = self.cell_neighbors(y_n, 0, cell_y)
            cell_l, cell_r, x0, x1 = self.cell_neighbors(x_n, 1, cell_x)
          
            # wrap the cell indices to get actual cell index
            cell_l = cell_l % self.cells[1]
            cell_r = cell_r % self.cells[1]
            cell_u = cell_u % self.cells[0]
            cell_d = cell_d % self.cells[0]

            # B field at each corner
            # indexing is row, col (y, x)
            b_00 = self.Bz[cell_u][cell_l] 
            b_01 = self.Bz[cell_u][cell_r] 
            b_10 = self.Bz[cell_d][cell_l] 
            b_11 = self.Bz[cell_d][cell_r] 

            # indices of neighboring cells
            self.b_particle[i] = self.interpolate(
                x_n, y_n, [x0, x1, y0, y1], [b_00, b_01, b_10, b_11])
            
            #print(self.b_particle[i])
            break
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
