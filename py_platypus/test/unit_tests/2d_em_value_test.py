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
   
        print("\n==== Testing Interpolate Ex field ====")
        # randomly set up electric field and interpolate field at each particle 
        self.pic.ex_edges = np.random.rand(self.pic.ex_edges.shape[0], self.pic.ex_edges.shape[1])
        self.pic.init_x_random() 
        self.pic.interpolate_e()
   
        # calculate the expected Ex field at the first N particles
        expected = np.zeros(self.num_random_test)
        for i in range(self.num_random_test):

            # randomly place N particles in a grid cell shifted half down
            dy, dx = self.pic.dx
            
            x_rand = np.random.rand()
            y_rand = np.random.rand()
            
            self.pic.electron_x[i] = x_rand * dx 
            self.pic.electron_y[i] = (y_rand + 0.5) * dy

            # calculate the expected electric field
            expected[i] = (self.pic.ex_edges[0][0] * (1 - x_rand) * (1 - y_rand) +\
                           self.pic.ex_edges[0][1] * x_rand * (1 - y_rand) +\
                           self.pic.ex_edges[1][0] * (1 - x_rand) * y_rand +\
                           self.pic.ex_edges[1][1] * x_rand * y_rand) 
        self.pic.interpolate_e()
           
        # mean electric field at edges and mean electric field at particles
        # should be similar
        mean_ex_edges = np.mean(self.pic.ex_edges)
        mean_ex_parts = np.mean(self.pic.e_particle[:, 0])
  
        self.assertTrue(
            np.isclose(mean_ex_edges, mean_ex_parts, rtol=1e-2, atol=1e-2), 
            "Mean ex edges: {} Mean ex particles: {}".format(mean_ex_edges, mean_ex_parts)
            )

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
    
        print("\n==== Testing Interpolate Ey field ====")
        self.pic.ey_edges = np.random.rand(self.pic.ey_edges.shape[0], self.pic.ey_edges.shape[1])
        self.pic.init_x_random() 
        self.pic.interpolate_e()
        
        # calculate the expected Ey field at the first N particles
        expected = np.zeros(self.num_random_test)
        for i in range(self.num_random_test):

            # randomly place N particles in a grid cell shifted half right
            dy, dx = self.pic.dx
            
            x_rand = np.random.rand()
            y_rand = np.random.rand()
            
            self.pic.electron_x[i] = (x_rand + 0.5) * dx 
            self.pic.electron_y[i] = y_rand * dy

            # calculate the expected electric field
            expected[i] = (self.pic.ey_edges[0][0] * (1 - x_rand) * (1 - y_rand) +\
                           self.pic.ey_edges[0][1] * x_rand * (1 - y_rand) +\
                           self.pic.ey_edges[1][0] * (1 - x_rand) * y_rand +\
                           self.pic.ey_edges[1][1] * x_rand * y_rand) 
        self.pic.interpolate_e()
           
        # mean electric field at edges and mean electric field at particles
        # should be similar
        mean_ey_edges = np.mean(self.pic.ey_edges)
        mean_ey_parts = np.mean(self.pic.e_particle[:, 1])

        # TODO these are not as close as the ex field for no apparent reason
        self.assertTrue(
            np.isclose(mean_ey_edges, mean_ey_parts, rtol=1e-2, atol=1e-2), 
            "Mean ey edges: {} Mean ey particles: {}".format(mean_ey_edges, mean_ey_parts)
            )

        for i in range(self.num_random_test):
            with self.subTest(i=i):
                self.assertTrue(
                    np.isclose(expected[i], self.pic.e_particle[i][1]), 
                    "Expected: {} Actual: {}".format(expected[i], self.pic.e_particle[i][1])
                )
    
        return
    
    def test_interpolate_b(self):
        '''
        test interpolating the magnetic field at several test particles
        from the field values
        '''
    
        print("\n==== Testing Interpolate B field ====")
        self.pic.bz = np.random.rand(self.pic.bz.shape[0], self.pic.bz.shape[1])
        self.pic.init_x_random()
      
        # calculate the expected Bz field at the first N particles
        expected = np.zeros(self.num_random_test)
        for i in range(self.num_random_test):

            # randomly place N particles in a grid cell shifted half right
            dy, dx = self.pic.dx
            
            x_rand = np.random.rand()
            y_rand = np.random.rand()
            
            self.pic.electron_x[i] = x_rand * dx 
            self.pic.electron_y[i] = y_rand * dy

            # calculate the expected electric field
            expected[i] = (self.pic.bz[0][0] * (1 - x_rand) * (1 - y_rand) +\
                           self.pic.bz[0][1] * x_rand * (1 - y_rand) +\
                           self.pic.bz[1][0] * (1 - x_rand) * y_rand +\
                           self.pic.bz[1][1] * x_rand * y_rand) 
        
        self.pic.interpolate_b() 
        mean_bz_grid = np.mean(self.pic.bz) 
        mean_bz_parts = np.mean(self.pic.b_particle)
           

        self.assertTrue(
            np.isclose(mean_bz_grid, mean_bz_parts, rtol=1e-3, atol=1e-3), 
            "Mean bz edges: {} Mean bz particles: {}".format(mean_bz_grid, mean_bz_parts)
            )

        
        for i in range(self.num_random_test):
            np.isclose(expected[i], self.pic.bz[i]), 
            with self.subTest(i=i):
                self.assertTrue(
                    np.isclose(expected[i], self.pic.b_particle[i]), 
                    "Expected: {} Actual: {}".format(expected[i], self.pic.b_particle[i])
                    )
          
        return

if __name__ == "__main__":
    unittest.main()
