"""
Unit tests for the 2D electromagnetic PIC that test functions that calculate
deterministic values.
"""
import unittest
import py_platypus as plat
import numpy as np
import matplotlib.pyplot as plt
import unit_helpers



class Pic2dEmValueTester(unittest.TestCase):
    """
    Unit tests for electromagnetic 2D PIC that test calculated values.
    """
    
    def setUp(self):
        """
        set up function called before each test
        """

        self.pic = unit_helpers.setup("2d_em") 
        self.num_random_test = 10

    def test_interpolate(self):
        '''
        Test interpolating field properties at particles
        '''
    
        print("==== Testing Interpolate ====")
        corners = [0, 1, 0, 1]     # [x0, x1, y0, y1] coordinates of corners
        x_ns = [0.5, 0.75]         # list of x_n particle positions
        y_ns = [0.5, 0.75]         # list of y_n particle positions
        values = [1, 2, 3, 4]      # field values at the four corners
    
        # expected values
        expected = [0.25 * 1 + 0.25 * 2 + 0.25 * 3 + 0.25 * 4, 
                    1/16 * 1 + 3/16 * 2 + 3/16 * 3 + 9/16 * 4]
        
        # scale expected values by 1/ area of cells 
        expected = 1 / np.prod(self.pic.dx) * np.array(expected)
   

        for i in range(len(expected)):
        # iterate through test cases
            with self.subTest(i=i): 
                # calculate interpolated values 
                int_value = self.pic.interpolate(x_ns[i], y_ns[i], corners, values) 
                self.assertEqual(expected[i], int_value) 
    
        return
    
    def test_interpolate_ex(self):
        '''
        test interpolating the electric field at each particle from the field
        values
        '''
   
        print("==== Testing Interpolate Ex field ====")
        # randomly set up electric field and interpolate field at each particle 
        self.pic.ex_edges = np.random.rand(self.pic.ex_edges.shape[0], self.pic.ex_edges.shape[1])
        self.pic.init_x_random() 
        self.pic.interpolate_e()
   
        # calculate the expected Ex field at the first N particles
        expected = np.zeros(self.num_random_test)
        for i in range(self.num_random_test):

            # randomly place N particles in a shifted quarter grid cell
            offset = np.random.rand()
            dy, dx = self.pic.dx
            
            x_rand = np.random.rand()
            y_rand = np.random.rand()
            
            self.pic.electron_x[i] = x_rand * dx 
            self.pic.electron_y[i] = (y_rand + 0.5) * dy

            # calculate the expected electric field
            # TODO CHECK THIS
            expected[i] = (self.pic.ex_edges[0][0] * (1 - x_rand) * (1 - y_rand) +\
                           self.pic.ex_edges[0][1] * x_rand * (1 - y_rand) +\
                           self.pic.ex_edges[1][0] * (1 - x_rand) * y_rand +\
                           self.pic.ex_edges[1][1] * x_rand * y_rand) 
        self.pic.interpolate_e()
           
        # mean electric field at edges and mean electric field at particles
        # should be similar
        mean_ex_edges = np.mean(self.pic.ex_edges)
        mean_ex_parts = np.mean(self.pic.e_particle[:, 0])))
   
        self.assertTrue(np.isclose(mean_ex_edges, mean_ex_parts))

        for i in range(self.num_random_test):
            with self.subTest(i=i):
                self.assertTrue(
                    np.isclose(expected[i], self.pic.e_particle[i][0]), 
                    "Expected: {} Actual: {}".format(expected[i], self.pic.e_particle[i][0])
                )
    
        return
    
    def test_interpolate_ey(self):
        '''test interpolating the electric field at each particle from the field
        values'''
    
        return
        print("==== Testing Interpolate Ey field ====")
        self.pic.ey_edges = np.indices(self.pic.ey_edges.shape)[0]
        self.pic.init_x_random() 
        self.pic.interpolate_e()
    
        print("Average Ey field: {} Average Ey particle: {}".format(
            np.mean(self.pic.ey_edges), np.mean(self.pic.e_particle[:, 1])))
        
        expected = np.zeros(self.num_random_test)
        for i in range(self.num_random_test):
            offset = np.random.rand() / 2
            dy, dx = self.pic.dx
            self.pic.electron_x[i] = (1 + offset) * dx 
            self.pic.electron_y[i] = (1 - offset) * dy 
    
            short = 0.5 - offset
            expected[i] = (self.pic.ex_edges[0][0] * offset * short +\
                           self.pic.ex_edges[0][1] * (1 - short) * offset +\
                           self.pic.ex_edges[1][0] * (1 - offset) * short +\
                           self.pic.ex_edges[1][1] * (1 - offset) * (1 - short)) 
        self.pic.interpolate_e()
       
        for i in range(self.num_random_test):
            print("Expected: {} Actual: {}".format(expected[i], self.pic.e_particle[i][0]))
    
        return
    
    def test_interpolate_b(self):
        '''test interpolate the magnetic field at each particle from the field
        values'''
        return
    
        print("==== Testing Interpolate B field ====")
        self.pic.bz = np.indices(self.pic.bz.shape)[0]
        self.pic.init_x_random() 
        self.pic.interpolate_b() 
        
        print("Average B field: {} Average B particle: {}".format(
            np.mean(self.pic.bz), np.mean(self.pic.b_particle)))
       
        values = [[4, 1], [8, 2]]
    
        self.pic.bz[0][0] = values[0][0]
        self.pic.bz[0][1] = values[0][1]
        self.pic.bz[1][0] = values[1][0]
        self.pic.bz[1][1] = values[1][1]
    
        for i in range(5):
            offset = np.random.rand() / 2
            dy, dx = self.pic.dx
            self.pic.electron_x[0] = (1 + offset) * dx 
            self.pic.electron_y[0] = (1 - offset) * dy 
    
            short = 0.5 - offset
            expected = (self.pic.bz[0][0] * (1 - short) * short +\
                        self.pic.bz[0][1] * (1 - short) * (1 - short) +\
                        self.pic.bz[1][0] * short * short +\
                        self.pic.bz[1][1] * (1 - short) * short)
            self.pic.interpolate_b()
        
            print("Expected: {} Actual: {}".format(expected, self.pic.b_particle[0]))
        return

if __name__ == "__main__":
    unittest.main()
