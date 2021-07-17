# pic_2d.py
# 2D electrostatic particle in cell (PIC) plasma simulation
#

import numpy as np
import copy
from scipy import fft, ifft
import matplotlib.pyplot as plt

import py_platypus as plat

MIN_J = 1e-8    # minimum value for index J when building k array

class PIC_2D:
    def __init__(self, params):
        
        # random seed
        np.random.seed(params["seed"])
        self.params = params
        self.dimensions = 2
        self.step = 0                             # current step

        # domain parameters 
        self.dx = params["dx"]                    # size of cells
        self.dt = params["timestep"]
        self.steps = params["steps"]              # time steps to run for
        self.cells = params["cells"]              # number of cells in x direction
        self.nodes = [x + 1 for x in params["cells"]]
        self.n_particles = params["n_particles"]  # total number of particles
        self.xmax = self.dx[1] * self.cells[1]
        self.ymax = self.dx[0] * self.cells[0]
        
        self.particle_weight = 1/(self.n_particles/(self.cells[0] * self.cells[1]))  # density/particles per cell
        
        # state information
        self.electron_x  = np.zeros(self.n_particles) # electron positions
        self.electron_y  = np.zeros(self.n_particles) # electron positions
        self.electron_vx = np.zeros(self.n_particles) # electron velocities
        self.electron_vy = np.zeros(self.n_particles) # electron velocities
        
        self.ion_x  = np.zeros(self.n_particles)      # ion positions
        self.ion_y  = np.zeros(self.n_particles)      # ion positions
        self.ion_vx = np.zeros(self.n_particles)      # ion velocities
        self.ion_vy = np.zeros(self.n_particles)      # ion velocities
        self.ion_e = np.zeros(self.n_particles)       # e-field at particles
  
        # electron and ion number density at each cell
        self.ne  = np.zeros(self.cells)  
        self.ni  = np.zeros(self.cells)  
        # charge density at each cell center
        self.rho = np.zeros(self.cells)  
        # potential at cell centers
        self.phi = np.zeros(self.cells)  
        self.batch = []                       # batch of particles to follow

        # field quantities on nodes 
        self.ex = np.zeros(self.nodes)  # electric field Ex at each node
        self.ey = np.zeros(self.nodes)  # electric field Ey at each node

        # list of dictionaries holding output values
        self.output = {"electrostatic_energy" :[], "kinetic_energy": [], "batch_ke": []}

    def init_x_random(self):
        '''randomly initialize the positions of the macroparticles'''
        self.electron_x = np.random.rand(self.n_particles) * self.xmax
        self.electron_y = np.random.rand(self.n_particles) * self.ymax
        
        self.ion_x = np.random.rand(self.n_particles) * self.xmax
        self.ion_y = np.random.rand(self.n_particles) * self.ymax
        
        return

    def init_v_maxwellian(self):
        '''initializes the velocity distribution function as a maxwellian'''
       
        # confirm that 2D maxwell is maxwell in each dimension
        for i in range(self.n_particles):
            r1x = max(1e-8, np.random.rand())
            r2x = np.random.rand()
            
            r1y = max(1e-8, np.random.rand())
            r2y = np.random.rand()
            
            self.electron_vx[i] = np.sqrt(-np.log(r1x)) * np.cos(2 * np.pi * r2x)
            self.ion_vx[i] = 0
            
            self.electron_vy[i] = np.sqrt(-np.log(r1y)) * np.cos(2 * np.pi * r2y)
            self.ion_vy[i] = 0
        
        return
    
    def init_v_two_stream(self):
        '''initializes the velocity distribution of electrons as two 
        counter propagating streams. This function should be called after
        the x and v have already been initialized
        inputs: vpos - normalized velocity of positive stream
                vneg - normalized velocity of negative stream
                fraction - fraction of particles in the width part of streams'''

        vpos         = self.params["two_stream"]["vpos"]
        vneg         = self.params["two_stream"]["vneg"]
        fraction     = self.params["two_stream"]["stream_frac"]
        stream_width = self.params["two_stream"]["stream_width"]
        
        stream_y_start = self.ymax/2 - stream_width/2   # start and stop of stream
        stream_y_stop = self.ymax/2 + stream_width/2

        # iterate through particles and set the velocities
        for i in range(self.n_particles):
            if self.electron_y[i] < stream_y_stop and self.electron_y[i] > stream_y_start:
                r = np.random.rand()

                # particle is in the positive stream
                if r < fraction / 2:
                    self.electron_vx[i] = vpos 
                    self.electron_vy[i] = 0

                # particle is in the negative stream
                if r > fraction / 2 and r < fraction:
                    self.electron_vx[i] = vneg 
                    self.electron_vy[i] = 0
        
        return

    def init_v_single_stream(self):
        '''randomly sets a certain fraction of electrons to an identical
        velocity, simulating a single stream
        inputs: fraction - percent of particles to set velocity
                v - normalized velocity
                stream_width - width of the stream'''
        
        fraction = self.params["single_stream"]["stream_frac"]
        stream_width = self.params["single_stream"]["stream_width"] 
        v = self.params["single_stream"]["stream_v"]
        
        stream_y_start = self.ymax/2 - stream_width/2   # start and stop of stream
        stream_y_stop = self.ymax/2 + stream_width/2
       
        for i in range(self.n_particles):
            if self.electron_y[i] < stream_y_stop and self.electron_y[i] > stream_y_start:
                r = np.random.rand()
                if r < fraction:
                    self.electron_vx[i] = v
                    self.electron_vy[i] = 0
        return

    def density_perturbation(self):
        '''create a sinusoidal density perturbation along the x-axis
        delta_n - perturbation amplitude
        k - k wave number of perturbation'''
        
        delta_n = self.params["landau"]["amplitude"]
        k = self.params["landau"]["mode"]
        
        for i in range(self.n_particles):
            delta_x = delta_n/k * np.sin(k * self.electron_x[i])
            self.electron_x[i] += delta_x
            while self.electron_x[i] < 0:
                self.electron_x[i] += self.xmax

            while self.electron_x[i] > self.xmax:
                self.electron_x[i] -= self.xmax

    def update_ni(self):
        '''update the ion density in each cell'''
        self.update_n("ion")
        return
    
    def update_ne(self):
        '''update the electron density in each cell'''
        self.update_n("electron") 
        return

    def cell_neighbors(self, x, dim):
        '''calculates the indices of the two adjacent cells that a particle
        weight should be distributed between and the center of the cells. Note
        that the value returned may include "phantom cells" past the boundary
        of the system

        inputs: x - position
                dim - index of dimension (0 for x, 1 for y, 2 for z)
                cell - index of cell the
        outputs: cell_0 - index of lower cell the weight should be distributed
                 cell_1 - index of upper cell
                 cell_0_center - center of cell_0 along specified dimension
                 cell_1_center - center of cell_1 along specified dimension'''
            
        dx = self.dx[plat.DIM_MAP[dim]]       # cell size along this dimension
        cells = self.cells[plat.DIM_MAP[dim]] # number of cells along this dimension
        cell_in = int(np.floor(x/dx))         # index of cell particle is in

        # particle is beyond the cell midpoint
        if x > cell_in * dx + 0.5 * dx:
            cell_0 = cell_in 
            cell_1 = cell_in + 1
        # particle is not past the cell midpoint
        else:
            cell_0 = cell_in - 1
            cell_1 = cell_in
            
        # center of neighboring cells with lower and greater indices
        cell_0_center = cell_0 * dx + 0.5 * dx
        cell_1_center = cell_1 * dx + 0.5 * dx
        
        return cell_0, cell_1, cell_0_center, cell_1_center

    def update_n(self, particle_type):
        '''update the particle density
        particle_type (str) - "ion" or "electron" '''
     
        # created shallow copy of the particle type we're interested in
        if particle_type == "electron":
            particle_x = self.electron_x
            particle_y = self.electron_y
            self.ne = np.zeros(self.cells)
            densities = self.ne
        elif particle_type == "ion":
            particle_x = self.ion_x
            particle_y = self.ion_y
            self.ni = np.zeros(self.cells)
            densities = self.ni
        else:
            raise ValueError("Unrecognized particle type: ".format(particle_type))

        for i in range(self.n_particles):
            x_n = particle_x[i]
            y_n = particle_y[i]
            
            # find the cell the particle is in
            cell_y = int(np.floor(y_n/self.dx[0]))
            cell_x = int(np.floor(x_n/self.dx[1]))
            
            # find indices of cells to the left and right that the weight
            # will be distributed between
            cell_u, cell_d, y0, y1 = self.cell_neighbors(y_n, "y")
            cell_l, cell_r, x0, x1 = self.cell_neighbors(x_n, "x")

            # calculate area of each rectangle
            area_upper_left  = plat.math_utils.points_to_area((x_n, y_n), (x0, y0))
            area_upper_right = plat.math_utils.points_to_area((x_n, y_n), (x1, y0))
            area_lower_left  = plat.math_utils.points_to_area((x_n, y_n), (x0, y1))
            area_lower_right = plat.math_utils.points_to_area((x_n, y_n), (x1, y1))

            # total area of a cell
            total_area = self.dx[0] * self.dx[1]

            # calculate weight to be distributed to each quadrant
            weight_upper_left  = area_lower_right/total_area * self.particle_weight 
            weight_upper_right = area_lower_left/total_area * self.particle_weight
            weight_lower_left  = area_upper_right/total_area * self.particle_weight 
            weight_lower_right = area_upper_left/total_area * self.particle_weight
            # get actual cell index, accounting for wraparound
            
            cell_l = cell_l % self.cells[1]
            cell_r = cell_r % self.cells[1]
            cell_u = cell_u % self.cells[0]
            cell_d = cell_d % self.cells[0]
            
            # update the densities
            densities[cell_u][cell_l] += weight_upper_left
            densities[cell_u][cell_r] += weight_upper_right
            densities[cell_d][cell_l] += weight_lower_left
            densities[cell_d][cell_r] += weight_lower_right
       
        return

    def update_rho(self):
        '''update the charge density'''
        raw_rho = self.ni - self.ne            # charge density
        self.rho = raw_rho - np.mean(raw_rho)  # normalize charge density
        return

    def update_phi(self):
        '''update the electric potential at each cell center'''
        
        R = np.fft.fft2(-self.rho)                  # fft of rho deviation 

        # build intermediate k array
        k = np.zeros(self.cells)
        for i in range(self.cells[0]):
            for j in range(self.cells[1]):
                coeff   = (np.sin(np.pi * i/self.cells[0])/self.dx[0]) ** 2 +\
                          (np.sin(np.pi * j/self.cells[1])/self.dx[1]) ** 2

                k[i][j] = -4 * max(MIN_J, coeff)
        
        Y = R/k                      # divide Fourier transform of rho by coef
        Y_hat = np.fft.ifft2(Y)      # take inverse Fourier transform
        potential = np.real(Y_hat)   # potential is the real part
        avg_potential = np.mean(potential)
        self.phi = (potential - avg_potential)
        
        
        return

    def update_e(self):
        '''update electric field at each node'''
        
        rows = self.cells[0]
        cols = self.cells[1]

        for i in range(self.nodes[0]):
            for j in range(self.nodes[1]):
                
                ul_potential = self.phi[(i - 1) % rows][(j - 1) % cols]
                ur_potential = self.phi[(i - 1) % rows][j % cols]
                dl_potential = self.phi[i % rows][(j - 1) % cols]
                dr_potential = self.phi[i % rows][j % cols]

                # iterpolate potential adjacent to the node in each direction
                up_potential    = np.mean([ul_potential, ur_potential]) 
                down_potential  = np.mean([dl_potential, dr_potential])
                left_potential  = np.mean([ul_potential, dl_potential])
                right_potential = np.mean([ur_potential, dr_potential])

                # E = -(phi_i - phi_i-1)/dx
                self.ex[i][j] = -(right_potential - left_potential)/self.dx[1] 
                self.ey[i][j] = -(down_potential - up_potential)/self.dx[0] 
       
        return
    
    def update_v(self):
        '''update velocity of particles based on electric fields'''
        for i in range(self.n_particles):
            x_n = self.electron_x[i]
            y_n = self.electron_y[i]

            # indices of neighboring nodes
            node_left  = int(np.floor(x_n/self.dx[1]))
            node_right = int(np.ceil(x_n/self.dx[1]))
            
            node_up   = int(np.floor(y_n/self.dx[0]))
            node_down = int(np.ceil(y_n/self.dx[0]))

            # coordinates of surrounding nodes
            x0 = node_left * self.dx[1]
            x1 = node_right * self.dx[1]
            y0 = node_up * self.dx[0]
            y1 = node_down * self.dx[0]

            # area of each rectangle
            area_upper_left  = plat.math_utils.points_to_area((x_n, y_n), (x0, y0))
            area_upper_right = plat.math_utils.points_to_area((x_n, y_n), (x1, y0))
            area_lower_left  = plat.math_utils.points_to_area((x_n, y_n), (x0, y1))
            area_lower_right = plat.math_utils.points_to_area((x_n, y_n), (x1, y1))

            # total area of a cell
            total_area = self.dx[0] * self.dx[1]

            # calculate weight to be distributed to each quadrant
            weight_00 = area_lower_right/total_area 
            weight_01 = area_lower_left/total_area
            weight_10 = area_upper_right/total_area 
            weight_11 = area_upper_left/total_area

            # electric field at each corner
            e_00 = [self.ex[node_up][node_left], self.ey[node_up][node_left]]
            e_01 = [self.ex[node_up][node_right], self.ey[node_up][node_right]]
            e_10 = [self.ex[node_down][node_left], self.ey[node_down][node_left]]
            e_11 = [self.ex[node_down][node_right], self.ey[node_down][node_right]]
            
            # list of weights and list of electric fields
            weights = np.array([weight_00, weight_01, weight_10, weight_11])
            e_nodes = np.array([e_00, e_01, e_10, e_11])

            e_particle = np.sum(list(map(lambda a, b: a * b, weights, e_nodes)), axis=0)
            self.electron_vx[i] -= e_particle[0] * self.dt
            self.electron_vy[i] -= e_particle[1] * self.dt

        return

    def update_x(self):
        '''update position of particles based on v_(n + 0.5)'''
        for i in range(self.n_particles):

            self.electron_x[i] += self.electron_vx[i] * self.dt
            self.electron_y[i] += self.electron_vy[i] * self.dt

            # particle past boundary condition; circular boundary 
            while self.electron_x[i] < 0:
                self.electron_x[i] += self.xmax

            while self.electron_x[i] > self.xmax:
                self.electron_x[i] -= self.xmax
            
            while self.electron_y[i] < 0:
                self.electron_y[i] += self.ymax

            while self.electron_y[i] > self.ymax:
                self.electron_y[i] -= self.ymax

        return

    def calc_bulk_u(self):
        '''calculate and save the bulk velocity'''
        # TODO
        return

    def calc_electrostatic_energy(self):
        '''calculate and save the electrostatic energy'''
        electrostatic_energy = 0 
        for i in range(self.cells[0]):
            for j in range(self.cells[1]):
                ex_cell = np.mean([
                    self.ex[i][j],
                    self.ex[i][j + 1],
                    self.ex[i + 1][j],
                    self.ex[i + 1][j + 1]
                ])
                
                ey_cell = np.mean([
                    self.ey[i][j],
                    self.ey[i][j + 1],
                    self.ey[i + 1][j],
                    self.ey[i + 1][j + 1]
                ])
               
                # U = 0.5 * epsilon * area * E^2
                electrostatic_energy += 0.5 * self.dx[0] * self.dx[1]\
                    * (ex_cell ** 2 + ey_cell ** 2)
        
        # save the value
        self.output["electrostatic_energy"].append(electrostatic_energy)
    
    def calc_kinetic_energy(self):
        '''calculate and save the kinetic energy'''
        ke_energy = 0.5 * self.particle_weight * \
            (sum(self.electron_vx * self.electron_vx) + \
             sum(self.electron_vy * self.electron_vy))
        
        ke_energy *= np.prod(self.dx)  # multiply by ratio of potential energy
                                       # to kinetic energy so total energy is
                                       # constant
        
        self.output["kinetic_energy"].append(ke_energy) 

        return
    
    def calc_batch_kinetic_energy(self):
        '''calculate and save the kinetic energy of the particles in the
        batch being studied'''
        
        ke_energy = 0.0 
        
        for i in self.batch: 
            ke_energy += 0.5 * self.particle_weight * \
                (self.electron_vx[i] * self.electron_vx[i] +\
                 self.electron_vy[i] * self.electron_vy[i])
        
        ke_energy *= np.prod(self.dx)  # multiply by ratio of potential energy
                                       # to kinetic energy so total energy is
                                       # constant

        self.output["batch_ke"].append(ke_energy) 
        return

    def step(self):
        '''run the simulation for a single step, updating all parameters;
           methods for saving outputs must be called separately'''
        self.update_ni()        # calculate e and i number densities
        self.update_ne()  
        self.update_rho()       # update charge density
        self.update_phi()       # calculate cell potential
        self.update_e()         # calculate electric field at nodes
        self.update_v()         # calculate velocity of each particle
        self.update_x()         # update positions
        self.step += 1


