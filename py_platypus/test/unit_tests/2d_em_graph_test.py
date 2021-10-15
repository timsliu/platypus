"""
Unit tests for the 2D electromagnetic PIC that may not be deterministic
and tests a graph output. Different quantities are computed and graphed,
and then the test asks for user input regarding whether the graphs are
correct.
"""


import unittest
import py_platypus as plat
import numpy as np
import matplotlib.pyplot as plt
import unit_helpers


class Pic2dEmGraphTester(unittest.TestCase):
    """
    Unit tests for electromagnetic 2D PIC that generates graphs.
    """

    def setUp(self):
        """
        set up function called before each test
        """

        self.pic = unit_helpers.setup("2d_em") 

    def test_random_x(self):
        """
        Tests initial random distribution of particles.
        """
        print("==== Testing random x ====")
        self.pic.init_x_random()
        plt.figure(1)
        
        plt.scatter(self.pic.electron_x, self.pic.electron_y, s=0.05)
        expect_str = "\nRandom distribution of points"

        self.assertEqual(unit_helpers.check_graph(expect_str), True, "random x image does not match expected")

        return

    def test_init_e(self):
        '''
        Tests calculating the inital electric field, both at the cell centers
        at initialization (from a Poisson solver) and interpolated to the
        edges.
        '''
        print("==== Test init E ====") 
        # create an initial density perturbation and calculate the electric field
        self.pic.init_x_random()
        self.pic.init_v_maxwellian()
        self.pic.density_perturbation()
        self.pic.init_E()
   

        # Ex electric field should be sinusoidal
        plt.figure(1)
        plt.imshow(self.pic.ex)
        plt.title("Ex field") 
        plt.colorbar()
    
        # Ey electric field should be random
        plt.figure(2) 
        plt.imshow(self.pic.ey)
        plt.title("Ey field") 
        plt.colorbar()
        
        # Ex electric field should be sinusoidal
        plt.figure(3)
        plt.imshow(self.pic.ex_edges)
        plt.title("Ex field edges") 
        plt.colorbar()
    
        # Ey electric field should be random
        plt.figure(4) 
        plt.imshow(self.pic.ey_edges)
        plt.title("Ey field edges") 
        plt.colorbar()
        
        expect_str = "\nEx and Ex field edges sinusoidal and similar; Ey and Ey fields edges random and similar"

        self.assertEqual(unit_helpers.check_graph(expect_str), True, "init E field image does not match expected")
    
        return
    
    def test_update_b(self):
        '''
        Tests calculating the magnetic field update based on the electric
        field.
        '''
   
        print("==== Test update B ====")
        # set up x vectors
        rows, cols = self.pic.ex_edges.shape
        for i in range(rows):
            for j in range(cols):
                x = i - rows/2
                y = j - cols/2
                
                self.pic.ex_edges[i][j] = y 
       
        # set up y vectors
        rows, cols = self.pic.ey_edges.shape
        for i in range(rows):
            for j in range(cols):
                x = i - rows/2
                y = j - cols/2
                
                self.pic.ey_edges[i][j] = -x 
    
        self.pic.calc_B_update()
        plt.figure(1)
        plt.imshow(self.pic.delta_bz)
        plt.title("B update for uniform curl")
        expect_str = "\nUniform graph"

        self.assertEqual(unit_helpers.check_graph(expect_str), True, "init B field image does not match expected")

if __name__ == "__main__":
    unittest.main()
