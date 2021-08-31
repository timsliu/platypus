# pic_3d.py
# 3D electrostatic particle in cell (PIC) plasma simulation
#

import numpy as np
import copy
import matplotlib.pyplot as plt

import py_platypus as plat

MIN_J = 1e-8    # minimum value for index J when building k array

class PIC_3D:
    def __init__(self, params):
        
        # random seed
        np.random.seed(params["seed"])
        self.params = params
        self.dimensions = 3

        # domain parameters 
        self.dx = params["dx"]                    # size of cells
        self.dt = params["timestep"]
        self.steps = params["steps"]              # time steps to run for
        self.cells = params["cells"]              # number of cells in each direction
        self.nodes = [x + 1 for x in params["cells"]]
        self.n_particles = params["n_particles"]  # total number of particles
        self.xmax = self.dx[1] * self.cells[1]    # self.dx is in matrix order (i, j, k)
        self.ymax = self.dx[0] * self.cells[0]
        self.zmax = self.dx[2] * self.cells[2]
        
        self.particle_weight = 1/(self.n_particles/(np.prod(self.cells)))  # density/particles per cell
        
        # state information
        self.electron_x  = np.zeros(self.n_particles) # electron positions
        self.electron_y  = np.zeros(self.n_particles) # electron positions
        self.electron_z  = np.zeros(self.n_particles) # electron positions
        
        self.electron_vx = np.zeros(self.n_particles) # electron velocities
        self.electron_vy = np.zeros(self.n_particles) # electron velocities
        self.electron_vz = np.zeros(self.n_particles) # electron velocities
        
        self.ion_x  = np.zeros(self.n_particles)      # ion positions
        self.ion_y  = np.zeros(self.n_particles)      # ion positions
        self.ion_z  = np.zeros(self.n_particles)      # ion positions
        self.ion_vx = np.zeros(self.n_particles)      # ion velocities
        self.ion_vy = np.zeros(self.n_particles)      # ion velocities
        self.ion_vz = np.zeros(self.n_particles)      # ion velocities
  
        # electron and ion number density at each cell
        # TODO this doesn't work if the number of cells in each axis is different
        self.ne  = np.zeros(self.cells)       # x, y, z 
        self.ni  = np.zeros(self.cells)       # x, y, z 
        # charge density at each cell center
        self.rho = np.zeros(self.cells)  
        # potential at cell centers
        self.phi = np.zeros(self.cells)  
        self.batch = []                       # batch of particles to follow

        # field quantities on nodes 
        self.ex = np.zeros(self.nodes)  # electric field Ex at each node
        self.ey = np.zeros(self.nodes)  # electric field Ey at each node
        self.ez = np.zeros(self.nodes)  # electric field Ez at each node

        # list of dictionaries holding output values
        self.output = {"electrostatic_energy" :[], "kinetic_energy": [], "batch_ke": []}

    def init_x_random(self):
        '''randomly initialize the positions of the macroparticles'''
        self.electron_x = np.random.rand(self.n_particles) * self.xmax
        self.electron_y = np.random.rand(self.n_particles) * self.ymax
        self.electron_z = np.random.rand(self.n_particles) * self.zmax
        
        self.ion_x = np.random.rand(self.n_particles) * self.xmax
        self.ion_y = np.random.rand(self.n_particles) * self.ymax
        self.ion_z = np.random.rand(self.n_particles) * self.zmax
        
        return

    def init_v_maxwellian(self):
        '''initializes the velocity distribution function as a maxwellian'''
       
        for i in range(self.n_particles):
            r1x = max(1e-8, np.random.rand())
            r2x = np.random.rand()
            
            r1y = max(1e-8, np.random.rand())
            r2y = np.random.rand()
            
            r1z = max(1e-8, np.random.rand())
            r2z = np.random.rand()
            
            self.electron_vx[i] = np.sqrt(-np.log(r1x)) * np.cos(2 * np.pi * r2x)
            self.ion_vx[i] = 0
            
            self.electron_vy[i] = np.sqrt(-np.log(r1y)) * np.cos(2 * np.pi * r2y)
            self.ion_vy[i] = 0
            
            self.electron_vz[i] = np.sqrt(-np.log(r1z)) * np.cos(2 * np.pi * r2z)
            self.ion_vz[i] = 0
        
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
        
        stream_z_start = self.zmax/2 - stream_width/2   # start and stop of stream
        stream_z_stop = self.zmax/2 + stream_width/2

        # iterate through particles and set the velocities
        for i in range(self.n_particles):
            # particle is in the stream
            if self.electron_y[i] < stream_y_stop and \
                self.electron_y[i] > stream_y_start and \
                self.electron_z[i] < stream_z_stop and \
                self.electron_z[i] > stream_z_start:
                
                r = np.random.rand()

                # particle is in the positive stream
                if r < fraction / 2:
                    self.electron_vx[i] = vpos 
                    self.electron_vy[i] = 0
                    self.electron_vz[i] = 0

                # particle is in the negative stream
                if r > fraction / 2 and r < fraction:
                    self.electron_vx[i] = vneg 
                    self.electron_vy[i] = 0
                    self.electron_vz[i] = 0
        
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
        
        stream_z_start = self.zmax/2 - stream_width/2   # start and stop of stream
        stream_z_stop = self.zmax/2 + stream_width/2
       
        for i in range(self.n_particles):
            if self.electron_y[i] < stream_y_stop and \
                self.electron_y[i] > stream_y_start and \
                self.electron_z[i] < stream_z_stop and \
                self.electron_z[i] > stream_z_start:
                r = np.random.rand()
                if r < fraction:
                    self.electron_vx[i] = v
                    self.electron_vy[i] = 0
                    self.electron_vz[i] = 0
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

    def cell_neighbors(self, x, dim, cell):
        '''calculates the indices of the two adjacent cells that a particle
        weight should be distributed between and the center of the cells. Note
        that the value returned may include "phantom cells" past the boundary
        of the system

        inputs: x - position
                dim - index of dimension (0 for x, 1 for y, 2 for z)
                cell - index of cell to return the neighbors of
        outputs: cell_0 - index of lower cell the weight should be distributed
                 cell_1 - index of upper cell
                 cell_0_center - center of cell_0 along specified dimension
                 cell_1_center - center of cell_1 along specified dimension'''
        
        dx = self.dx[dim]           # cell size along this dimension
        cells = self.cells[dim]     # number of cells along this dimension

        # particle is to the right of cell midpoint
        if x > cell * dx + 0.5 * dx:
            cell_0 = cell
            cell_1 = cell + 1
        # particle is to the left of cell midpoint
        else:
            cell_0 = cell - 1
            cell_1 = cell
            
        # center of left and right cells
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
            particle_z = self.electron_z

            self.ne = np.zeros(self.cells)
            densities = self.ne
        
        elif particle_type == "ion":
            particle_x = self.ion_x
            particle_y = self.ion_y
            particle_z = self.ion_z
            
            self.ni = np.zeros(self.cells)
            densities = self.ni
        else:
            raise ValueError("Unrecognized particle type: ".format(particle_type))

        for i in range(self.n_particles):
            x_n = particle_x[i]
            y_n = particle_y[i]
            z_n = particle_z[i]
            
            # find the cell the particle is in
            cell_y = int(np.floor(y_n/self.dx[0]))
            cell_x = int(np.floor(x_n/self.dx[1]))
            cell_z = int(np.floor(z_n/self.dx[2]))
            
            # find indices and cell centers of cells surrounding the particle 
            # that the weight will be distributed between
            cell_y_0, cell_y_1, y0, y1 = self.cell_neighbors(y_n, 0, cell_y)
            cell_x_0, cell_x_1, x0, x1 = self.cell_neighbors(x_n, 1, cell_x)
            cell_z_0, cell_z_1, z0, z1 = self.cell_neighbors(z_n, 2, cell_z)

            # calculate area of each rectangular prism
            vol_x0y0z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y0, z0))
            vol_x0y0z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y0, z1))
            vol_x0y1z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y1, z0))
            vol_x0y1z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y1, z1))
            vol_x1y0z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y0, z0))
            vol_x1y0z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y0, z1))
            vol_x1y1z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y1, z0))
            vol_x1y1z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y1, z1))

            # total area of a cell
            total_volume = np.prod(self.dx)

            # calculate weight to be distributed to each quadrant
            weight_x0y0z0 = vol_x1y1z1/total_volume * self.particle_weight 
            weight_x0y0z1 = vol_x1y1z0/total_volume * self.particle_weight
            weight_x0y1z0 = vol_x1y0z1/total_volume * self.particle_weight 
            weight_x0y1z1 = vol_x1y0z0/total_volume * self.particle_weight
            weight_x1y0z0 = vol_x0y1z1/total_volume * self.particle_weight 
            weight_x1y0z1 = vol_x0y1z0/total_volume * self.particle_weight
            weight_x1y1z0 = vol_x0y0z1/total_volume * self.particle_weight 
            weight_x1y1z1 = vol_x0y0z0/total_volume * self.particle_weight


            # get actual cell index, accounting for wraparound
            cell_y_0 = cell_y_0 % self.cells[0]  
            cell_y_1 = cell_y_1 % self.cells[0] 
            cell_x_0 = cell_x_0 % self.cells[1]
            cell_x_1 = cell_x_1 % self.cells[1]
            cell_z_0 = cell_z_0 % self.cells[2]
            cell_z_1 = cell_z_1 % self.cells[2]

            # update the densities
            densities[cell_y_0][cell_x_0][cell_z_0] += weight_x0y0z0
            densities[cell_y_0][cell_x_0][cell_z_1] += weight_x0y0z1
            densities[cell_y_0][cell_x_1][cell_z_0] += weight_x0y1z0
            densities[cell_y_0][cell_x_1][cell_z_1] += weight_x0y1z1
            densities[cell_y_1][cell_x_0][cell_z_0] += weight_x1y0z0
            densities[cell_y_1][cell_x_0][cell_z_1] += weight_x1y0z1
            densities[cell_y_1][cell_x_1][cell_z_0] += weight_x1y1z0
            densities[cell_y_1][cell_x_1][cell_z_1] += weight_x1y1z1
       
        return

    def update_rho(self):
        '''update the charge density'''
        raw_rho = self.ni - self.ne            # charge density
        self.rho = raw_rho - np.mean(raw_rho)  # normalize charge density
        return

    def update_phi(self):
        '''update the electric potential at each cell center'''
        
        R = np.fft.fftn(-self.rho)                  # fft of rho deviation 

        # build intermediate k array
        kappa = np.zeros(self.cells)
        for i in range(self.cells[0]):
            for j in range(self.cells[1]):
                for k in range(self.cells[2]):
                    coeff   = (np.sin(np.pi * i/self.cells[0])/self.dx[0]) ** 2 +\
                              (np.sin(np.pi * j/self.cells[1])/self.dx[1]) ** 2 +\
                              (np.sin(np.pi * k/self.cells[2])/self.dx[2]) ** 2

                    kappa[i][j][k] = -4 * max(MIN_J, coeff)
        
        Y = R/kappa                  # divide Fourier transform of rho by coef
        Y_hat = np.fft.ifftn(Y)      # take inverse Fourier transform
        potential = np.real(Y_hat)   # potential is the real part
        avg_potential = np.mean(potential)
        self.phi = (potential - avg_potential)
        
        return

    def update_e(self):
        '''update electric field at each node'''
        
        rows = self.cells[0]
        cols = self.cells[1]
        lays = self.cells[2]    # number of layers

        for i in range(self.nodes[0]):
            for j in range(self.nodes[1]):
                for k in range(self.nodes[2]):

                    # look up potential in the 8 cells surrounding the node 
                    x0y0z0_potential = self.phi[(i - 1) % rows][(j - 1) % cols][(k - 1) % lays] 
                    x0y0z1_potential = self.phi[(i - 1) % rows][(j - 1) % cols][k % lays]
                    x1y0z0_potential = self.phi[(i - 1) % rows][j % cols][(k - 1) % lays]
                    x1y0z1_potential = self.phi[(i - 1) % rows][j % cols][k % lays]
                    x0y1z0_potential = self.phi[i % rows][(j - 1) % cols][(k - 1) % lays]
                    x0y1z1_potential = self.phi[i % rows][(j - 1) % cols][k % lays]
                    x1y1z0_potential = self.phi[i % rows][j % cols][(k - 1) % lays]
                    x1y1z1_potential = self.phi[i % rows][j % cols][k % lays]

                    # iterpolate potential adjacent to the node in each direction
                    x0_potential = np.mean([x0y0z0_potential, x0y0z1_potential, x0y1z0_potential, x0y1z1_potential]) 
                    x1_potential = np.mean([x1y0z0_potential, x1y0z1_potential, x1y1z0_potential, x1y1z1_potential])

                    y0_potential = np.mean([x0y0z0_potential, x0y0z1_potential, x1y0z0_potential, x1y0z1_potential])
                    y1_potential = np.mean([x0y1z0_potential, x0y1z1_potential, x1y1z0_potential, x1y1z1_potential])

                    z0_potential = np.mean([x0y0z0_potential, x0y1z0_potential, x1y0z0_potential, x1y1z0_potential])
                    z1_potential = np.mean([x0y0z1_potential, x0y1z1_potential, x1y0z1_potential, x1y1z1_potential])
                    
                    # E = -(phi_i - phi_i-1)/dx
                    self.ex[i][j][k] = -(x1_potential - x0_potential)/self.dx[1] 
                    self.ey[i][j][k] = -(y1_potential - y0_potential)/self.dx[0] 
                    self.ez[i][j][k] = -(z1_potential - z0_potential)/self.dx[2] 
       
        return
    
    def update_v(self):
        '''update velocity of particles based on electric fields'''
        
        # iterate through all particles
        for i in range(self.n_particles):
            x_n = self.electron_x[i]
            y_n = self.electron_y[i]
            z_n = self.electron_z[i]

            # indices of neighboring nodes
            node_x0 = int(np.floor(x_n/self.dx[1]))
            node_x1 = int(np.ceil(x_n/self.dx[1]))
            
            node_y0 = int(np.floor(y_n/self.dx[0]))
            node_y1 = int(np.ceil(y_n/self.dx[0]))
            
            node_z0 = int(np.floor(z_n/self.dx[2]))
            node_z1 = int(np.ceil(z_n/self.dx[2]))

            # coordinates of surrounding nodes
            x0 = node_x0 * self.dx[1]
            x1 = node_x1 * self.dx[1]
            y0 = node_y0 * self.dx[0]
            y1 = node_y1 * self.dx[0]
            z0 = node_z0 * self.dx[2]
            z1 = node_z1 * self.dx[2]

            # calculate area of each rectangular prism
            vol_x0y0z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y0, z0))
            vol_x0y0z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y0, z1))
            vol_x0y1z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y1, z0))
            vol_x0y1z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x0, y1, z1))
            vol_x1y0z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y0, z0))
            vol_x1y0z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y0, z1))
            vol_x1y1z0  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y1, z0))
            vol_x1y1z1  = plat.math_utils.points_to_volume((x_n, y_n, z_n), (x1, y1, z1))

            # total area of a cell
            total_volume = np.prod(self.dx)

            # calculate weight to be distributed to each quadrant
            weight_x0y0z0 = vol_x1y1z1/total_volume
            weight_x0y0z1 = vol_x1y1z0/total_volume
            weight_x0y1z0 = vol_x1y0z1/total_volume
            weight_x0y1z1 = vol_x1y0z0/total_volume
            weight_x1y0z0 = vol_x0y1z1/total_volume
            weight_x1y0z1 = vol_x0y1z0/total_volume
            weight_x1y1z0 = vol_x0y0z1/total_volume
            weight_x1y1z1 = vol_x0y0z0/total_volume

            # electric field [Ex, Ey, Ez] at each node surrounding the point
            e_x0y0z0 = [self.ex[node_y0][node_x0][node_z0], self.ey[node_y0][node_x0][node_z0], self.ez[node_y0][node_x0][node_z0]]
            e_x0y0z1 = [self.ex[node_y0][node_x0][node_z1], self.ey[node_y0][node_x0][node_z1], self.ez[node_y0][node_x0][node_z1]]
            e_x0y1z0 = [self.ex[node_y1][node_x0][node_z0], self.ey[node_y1][node_x0][node_z0], self.ez[node_y1][node_x0][node_z0]]
            e_x0y1z1 = [self.ex[node_y1][node_x0][node_z1], self.ey[node_y1][node_x0][node_z1], self.ez[node_y1][node_x0][node_z1]]
            e_x1y0z0 = [self.ex[node_y0][node_x1][node_z0], self.ey[node_y0][node_x1][node_z0], self.ez[node_y0][node_x1][node_z0]]
            e_x1y0z1 = [self.ex[node_y0][node_x1][node_z1], self.ey[node_y0][node_x1][node_z1], self.ez[node_y0][node_x1][node_z1]]
            e_x1y1z0 = [self.ex[node_y1][node_x1][node_z0], self.ey[node_y1][node_x1][node_z0], self.ez[node_y1][node_x1][node_z0]]
            e_x1y1z1 = [self.ex[node_y1][node_x1][node_z1], self.ey[node_y1][node_x1][node_z1], self.ez[node_y1][node_x1][node_z1]]

            # list of weights and list of electric fields
            weights = np.array([
                weight_x0y0z0, weight_x0y0z1, 
                weight_x0y1z0, weight_x0y1z1, 
                weight_x1y0z0, weight_x1y0z1, 
                weight_x1y1z0, weight_x1y1z1])

            e_nodes = np.array([
                e_x0y0z0, e_x0y0z1, e_x0y1z0, e_x0y1z1, 
                e_x1y0z0, e_x1y0z1, e_x1y1z0, e_x1y1z1])

            # averaged electric field at the particle
            e_particle = np.sum(list(map(lambda a, b: a * b, weights, e_nodes)), axis=0)
            self.electron_vx[i] -= e_particle[0] * self.dt
            self.electron_vy[i] -= e_particle[1] * self.dt
            self.electron_vz[i] -= e_particle[2] * self.dt

        return

    def update_x(self):
        '''update position of particles based on v_(n + 0.5)'''
        for i in range(self.n_particles):

            self.electron_x[i] += self.electron_vx[i] * self.dt
            self.electron_y[i] += self.electron_vy[i] * self.dt
            self.electron_z[i] += self.electron_vz[i] * self.dt

            # particle past boundary condition; circular boundary 
            while self.electron_x[i] < 0:
                self.electron_x[i] += self.xmax

            while self.electron_x[i] > self.xmax:
                self.electron_x[i] -= self.xmax
            
            while self.electron_y[i] < 0:
                self.electron_y[i] += self.ymax

            while self.electron_y[i] > self.ymax:
                self.electron_y[i] -= self.ymax
            
            while self.electron_z[i] < 0:
                self.electron_z[i] += self.zmax

            while self.electron_z[i] > self.zmax:
                self.electron_z[i] -= self.zmax

        return

    def calc_electrostatic_energy(self):
        '''calculate and save the electrostatic energy'''
        electrostatic_energy = 0
        
        # iterator over all 3 digit binary values representing neighbor indices
        array_3d = np.ones([2 for x in range(self.dimensions)])
        neighbor_it = np.nditer(array_3d, flags=['multi_index'])
        neighbor_indices = [neighbor_it.multi_index for i in neighbor_it]

        # iterator over indices of all cells
        cell_it = np.nditer(np.ones(self.cells), flags=['multi_index'])

        # iterate over all cells
        for x in cell_it:

            # list of the electric field of the 8 nodes bordering the cell
            ex_neighbors = []
            ey_neighbors = []
            ez_neighbors = []
           
            # iterate through relative indices of the 8 nodes around the cell
            for neighbor_index in neighbor_indices:
                # get absolute index
                j_neighbor, i_neighbor, k_neighbor = np.array(neighbor_index)\
                                                     + np.array(cell_it.multi_index)

                ex_neighbors.append(self.ex[j_neighbor][i_neighbor][k_neighbor])
                ey_neighbors.append(self.ey[j_neighbor][i_neighbor][k_neighbor])
                ez_neighbors.append(self.ez[j_neighbor][i_neighbor][k_neighbor])
       
            ex_cell = np.mean(ex_neighbors)
            ey_cell = np.mean(ey_neighbors)
            ez_cell = np.mean(ez_neighbors)

            # U = 0.5 * epsilon * volume * E^2
            electrostatic_energy += 0.5 * np.prod(self.dx)\
                * (ex_cell ** 2 + ey_cell ** 2 + ez_cell ** 2)
        
        # save the value
        self.output["electrostatic_energy"].append(electrostatic_energy)

        return
    
    def calc_kinetic_energy(self):
        '''calculate and save the kinetic energy'''
        ke_energy = 0.5 * self.particle_weight * \
            (sum(self.electron_vx * self.electron_vx) + \
             sum(self.electron_vy * self.electron_vy) + \
             sum(self.electron_vz * self.electron_vz))
        
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
                 self.electron_vy[i] * self.electron_vy[i] +\
                 self.electron_vz[i] * self.electron_vz[i])
        
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


