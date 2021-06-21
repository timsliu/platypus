# pic_1d.py
# 1D particle in cell plasma simulation
#

import numpy as np
import copy
from scipy import fft, ifft
import matplotlib.pyplot as plt

MIN_J = 1e-8    # minimum value for index J when building k array

class PIC_1D:
    def __init__(self, params):
        # TODO verify it's a valid params set
        
        # random seed
        #np.random.seed(params["seed"])

        # domain parameters 
        self.dx = params["dx"]
        self.dt = params["timestep"]
        self.steps = params["steps"]              # time steps to run for
        self.cells = params["cells"]              # number of cells
        self.nodes = [x + 1 for x in params["cells"]]
        self.n_particles = params["n_particles"]  # total number of particles
        self.xmax = self.dx[0] * self.cells[0]
       
        self.particle_weight = 1/(self.n_particles/np.prod(self.cells))  # density/particles per cell
        
        # state information
        self.electron_x = np.zeros(self.n_particles) # electron positions
        self.electron_v = np.zeros(self.n_particles) # electron velocities
        self.electron_e = np.zeros(self.n_particles) # e-field at particles
        
        self.ion_x = np.zeros(self.n_particles)      # ion positions
        self.ion_v = np.zeros(self.n_particles)      # ion velocities
        self.ion_e = np.zeros(self.n_particles)      # e-field at particles
        
        self.ne = np.zeros(self.cells)        # electron number density at each cell
        self.ni = np.zeros(self.cells)        # electron number density at each cell
        self.rho = np.zeros(self.cells)       # charge density at each cell center
        self.phi = np.zeros(self.cells)       # potential at cell centers
        self.batch = []                       # batch of particles to follow

        # field quantities on nodes 
        self.e = np.zeros(self.nodes)          # electric field at each node

        # list of dictionaries holding output values
        self.output = {"electrostatic_energy" :[], "kinetic_energy": [], "batch_ke": []}

    def init_x_random(self):
        '''randomly initialize the positions of the macroparticles'''
        self.electron_x = np.random.rand(self.n_particles) * self.xmax
        self.ion_x = np.random.rand(self.n_particles) * self.xmax
        
        return

    def init_x_uniform(self):
        '''uniformly initialize the positions of the macroparticles'''
        self.electron_x = np.linspace(0, self.xmax, num=self.n_particles, endpoint=False) 
        self.ion_x = np.linspace(0, self.xmax, num=self.n_particles, endpoint=False) 
        return

    def init_v_maxwellian(self):
        '''initializes the velocity distribution function as a maxwellian'''
        for i in range(self.n_particles):
            r1 = max(1e-8, np.random.rand())
            r2 = np.random.rand()
            self.electron_v[i] = np.sqrt(-np.log(r1)) * np.cos(2 * np.pi * r2)
            self.ion_v[i] = 0
        
        return
    
    def init_v_two_beams(self, vpos, vneg):
        '''initializes the velocity distribution of electrons as two 
        counter propagating beams
        inputs: vpos - normalized velocity of positive beam
                vneg - normalized velocity of negative beam'''

        # randomly select which half is positive
        pos_particles = np.random.choice(
            range(self.n_particles), 
            size=int(self.n_particles/2), 
            replace = False)

        # iterate through particles and set the velocities
        for i in range(self.n_particles):
            if i in pos_particles:
                self.electron_v[i] = vpos
            else:
                self.electron_v[i] = vneg
            self.ion_v[i] = 0
        
        return

    def init_v_single_stream(self, fraction, v):
        '''randomly sets a certain fraction of electrons to an identical
        velocity, simulating a single stream
        inputs: fraction - percent of particles to set velocity
                v - normalized velocity'''
        
        # randomly select which half is positive
        stream_particles = np.random.choice(
            range(self.n_particles), 
            size=int(self.n_particles * fraction), 
            replace = False)

        self.batch = stream_particles
        # iterate through particles and set the velocities
        for i in range(self.n_particles):
            if i in stream_particles:
                self.electron_v[i] = v
        
        return

    def density_perturbation(self, delta_n, k):
        '''create a sinusoidal density perturbation
        delta_n - perturbation amplitude
        k - k wave number of perturbation'''

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

    def update_n(self, particle_type):
        '''update the particle density
        particle_type (str) - "ion" or "electron" '''
        
        # copy the particle array we're interested in
        if particle_type == "electron":
            particle_x = np.copy(self.electron_x)
        elif particle_type == "ion":
            particle_x = np.copy(self.ion_x)
        else:
            raise ValueError("Unrecognized particle type: ".format(particle_type))

        # clear the array of densities
        densities = np.zeros(self.cells[0])
        for x_n in particle_x:
            
            # cell the particle is in
            cell = int(np.floor(x_n/self.dx[0]))
            # find indices of cells to the left and right that the weight
            # will be distributed between

            # particle is to the right of cell midpoint
            if x_n > cell * self.dx[0] + 0.5 * self.dx[0]:
                cell_left = cell
                cell_right = cell + 1
            # particle is to the left of cell midpoint
            else:
                cell_left = cell - 1
                cell_right = cell
            
            # center of left and right cells
            cell_left_x = cell_left * self.dx[0] + 0.5 * self.dx[0]
            cell_right_x = cell_right * self.dx[0] + 0.5 * self.dx[0]
           
            # weight to be distributed to left and right cells
            weight_left = (cell_right_x - x_n)/self.dx[0] * self.particle_weight 
            weight_right = (x_n - cell_left_x)/self.dx[0] * self.particle_weight 
        
            # get actual cell index, accounting for wraparound
            cell_left = cell_left % self.cells[0]
            cell_right = cell_right % self.cells[0]

            densities[cell_left] += weight_left
            densities[cell_right] += weight_right
       
        # copy the cell densities to appropriate array
        if particle_type == "electron":
            self.ne = copy.deepcopy(densities)
        if particle_type == "ion":
            self.ni = copy.deepcopy(densities)

        return

    def update_rho(self):
        '''update the charge density'''
        raw_rho = self.ni - self.ne            # charge density
        self.rho = raw_rho - np.mean(raw_rho)  # normalize charge density
        return

    def update_phi(self):
        '''update the electric potential at each cell center'''
        #self.rho = np.sin(np.linspace(0, 2 * np.pi, self.cells))
        R = fft(-self.rho)                  # fft of rho deviation 

        # build intermediate k array
        k = np.zeros(self.cells[0])
        for j in range(self.cells[0]):
            k[j] = np.pi/self.dx[0] * max(j, MIN_J)/(self.cells[0]/2)
            if j >= self.cells[0]/2:
                k[j] -= 2 * np.pi/self.dx[0]
        
        # intermediate kappa array
        kappa = np.sin(k * self.dx[0]/2)/(self.dx[0]/2)
        # intermediate Y array
        Y = - R/(kappa * kappa)
        Y_hat = ifft(Y)
        potential = np.real(Y_hat)   # potential is the real part
        avg_potential = np.mean(potential)
        self.phi = (potential - avg_potential)

        return

    def update_e(self):
        '''update electric field at each node'''
        for i in range(self.nodes[0]):
            if i == 0:
                # use the left potential boundary condition 
                left_potential = self.phi[-1]
            else:
                left_potential = self.phi[i-1]

            if i == (self.nodes[0] - 1):
                # use the right potential boundary condition 
                right_potential = self.phi[0]
            else:
                right_potential = self.phi[i]

            # E = -(phi_i - phi_i-1)/dx
            self.e[i] = -(right_potential - left_potential)/self.dx[0] 
        
        return
    
    def update_v(self):
        '''update velocity of particles based on electric fields'''
        for i in range(self.n_particles):
            x_n = self.electron_x[i]

            # indices of left and right nodes
            node_left = int(np.floor(x_n/self.dx[0]))
            node_right = int(np.ceil(x_n/self.dx[0]))


            # electric field at the left and right nodes
            e_left = self.e[node_left]
            e_right = self.e[node_right]

            # position of left and right nodes
            x_left = node_left * self.dx[0]
            x_right = node_right * self.dx[0]
           
            # calculate electric field at particle and update velocity
            e_particle = (x_right - x_n)/self.dx[0] * e_left + (x_n - x_left)/self.dx[0] * e_right
            self.electron_v[i] -= e_particle * self.dt

    def update_x(self):
        '''update position of particles based on v_(n + 0.5)'''
        for i in range(self.n_particles):

            self.electron_x[i] += self.electron_v[i] * self.dt

            # particle past boundary condition; circular boundary 
            while self.electron_x[i] < 0:
                self.electron_x[i] += self.xmax

            while self.electron_x[i] > self.xmax:
                self.electron_x[i] -= self.xmax

        return

    def calc_bulk_u(self):
        '''calculate and save the bulk velocity'''
        # TODO
        return

    def calc_electrostatic_energy(self):
        '''calculate and save the electrostatic energy'''
        electrostatic_energy = 0 
        for i in range(self.cells[0]):
            e_cell = np.mean([self.e[i], self.e[i + 1]]) # average E field in cell
            electrostatic_energy += 0.5 * self.dx[0] * (e_cell ** 2)
        
        # save the value
        self.output["electrostatic_energy"].append(electrostatic_energy)
    
    def calc_kinetic_energy(self):
        '''calculate and save the kinetic energy'''
        ke_energy = 0.5 * self.particle_weight * sum(self.electron_v * self.electron_v)
        self.output["kinetic_energy"].append(ke_energy) 

        return
    
    def calc_batch_kinetic_energy(self):
        '''calculate and save the kinetic energy'''
        ke_energy = 0.0 
        for i in self.batch: 
            ke_energy += 0.5 * self.particle_weight * self.electron_v[i] * self.electron_v[i]

        self.output["batch_ke"].append(ke_energy) 
        return

    def step(self):
        '''run the simulation for a single step, updating all parameters;
           methods for saving outputs must be called separately'''
        self.update_ni()        # calculate e and i number densities
        #print("Ni: ", self.ni)
        self.update_ne()  
        #print("Ne: ", self.ne)
        self.update_rho()       # update charge density
        #print("Charge density:", self.rho)
        #print("Mean charge density:", np.mean(self.rho))
        self.update_phi()       # calculate cell potential
        #print("Potential: ", self.phi)
        self.update_e()         # calculate electric field at nodes
        #print("Electric field: ", self.e)
        self.update_v()         # calculate velocity of each particle
        self.update_x()         # update positions


    def spectate(self):
        '''print velocity, position, electric field of a particle'''

        print("x: {:.3f}, v: {:.3f}, e_left: {:.3f}, e_right: {:.3f}".format(
            float(self.electron_x[10]), 
            float(self.electron_v[10]),
            float(self.e[int(np.floor(self.electron_x[10]/self.dx[0]))]),
            float(self.e[int(np.ceil(self.electron_x[10]/self.dx[0]))])))
    

