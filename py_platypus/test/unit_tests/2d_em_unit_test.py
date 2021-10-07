
import unittest
import py_platypus as plat
import numpy as np
import matplotlib.pyplot as plt
import unit_test_helpers



class Pic2dEmValueTester(unittest.TestCase):
    """
    Unit tests for electromagnetic 2D PIC that test calculated values.
    """
    
    def setUp(self):
        """
        set up function called before each test
        """

        self.pic = setup_2d_pic() 



def test_interpolate(pic):
    '''test helper function for interpolating field properties at particles'''

    print("==== Testing Interpolate ====")
    corners = [0, 1, 0, 1]     # [x0, x1, y0, y1] coordinates of corners
    x_ns = [0.5, 0.75]         # list of x_n particle positions
    y_ns = [0.5, 0.75]         # list of y_n particle positions
    values = [1, 2, 3, 4]      # field values at the four corners

    # expected values
    expected = [0.25 * 1 + 0.25 * 2 + 0.25 * 3 + 0.25 * 4, 
                1/16 * 1 + 3/16 * 2 + 3/16 * 3 + 9/16 * 4]
    
    # scale expected values by 1/ area of cells 
    expected = 1 / np.prod(pic.dx) * np.array(expected)

    # iterate through test cases
    for i in range(len(x_ns)):
        int_value = pic.interpolate(x_ns[i], y_ns[i], corners, values) 
        print("Expected: {} Actual: {}".format(expected[i], int_value))

    return

def test_interpolate_ex(pic):
    '''test interpolating the electric field at each particle from the field
    values'''

    print("==== Testing Interpolate Ex field ====")
    pic.ex_edges = np.random.rand(pic.ex_edges.shape[0], pic.ex_edges.shape[1])
    pic.init_x_random() 
    pic.interpolate_e()

    print("Average Ex field: {} Average Ex particle: {}".format(
        np.mean(pic.ex_edges), np.mean(pic.e_particle[:, 0])))
  
    expected = np.zeros(NUM_TESTS)
    for i in range(NUM_TESTS):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[i] = (1 + offset) * dx 
        pic.electron_y[i] = (1 - offset) * dy 

        short = 0.5 - offset
        expected[i] = (pic.ex_edges[0][0] * offset * short +\
                       pic.ex_edges[0][1] * (1 - short) * offset +\
                       pic.ex_edges[1][0] * (1 - offset) * short +\
                       pic.ex_edges[1][1] * (1 - offset) * (1 - short)) 
    pic.interpolate_e()
   
    for i in range(NUM_TESTS):
        print("Expected: {} Actual: {}".format(expected[i], pic.e_particle[i][0]))

    return

def test_interpolate_ey(pic):
    '''test interpolating the electric field at each particle from the field
    values'''

    print("==== Testing Interpolate Ey field ====")
    pic.ey_edges = np.indices(pic.ey_edges.shape)[0]
    pic.init_x_random() 
    pic.interpolate_e()

    print("Average Ey field: {} Average Ey particle: {}".format(
        np.mean(pic.ey_edges), np.mean(pic.e_particle[:, 1])))
    
    expected = np.zeros(NUM_TESTS)
    for i in range(NUM_TESTS):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[i] = (1 + offset) * dx 
        pic.electron_y[i] = (1 - offset) * dy 

        short = 0.5 - offset
        expected[i] = (pic.ex_edges[0][0] * offset * short +\
                       pic.ex_edges[0][1] * (1 - short) * offset +\
                       pic.ex_edges[1][0] * (1 - offset) * short +\
                       pic.ex_edges[1][1] * (1 - offset) * (1 - short)) 
    pic.interpolate_e()
   
    for i in range(NUM_TESTS):
        print("Expected: {} Actual: {}".format(expected[i], pic.e_particle[i][0]))

    return

def test_interpolate_b(pic):
    '''test interpolate the magnetic field at each particle from the field
    values'''

    print("==== Testing Interpolate B field ====")
    pic.bz = np.indices(pic.bz.shape)[0]
    pic.init_x_random() 
    pic.interpolate_b() 
    
    print("Average B field: {} Average B particle: {}".format(
        np.mean(pic.bz), np.mean(pic.b_particle)))
   
    values = [[4, 1], [8, 2]]

    pic.bz[0][0] = values[0][0]
    pic.bz[0][1] = values[0][1]
    pic.bz[1][0] = values[1][0]
    pic.bz[1][1] = values[1][1]

    for i in range(5):
        offset = np.random.rand() / 2
        dy, dx = pic.dx
        pic.electron_x[0] = (1 + offset) * dx 
        pic.electron_y[0] = (1 - offset) * dy 

        short = 0.5 - offset
        expected = (pic.bz[0][0] * (1 - short) * short +\
                    pic.bz[0][1] * (1 - short) * (1 - short) +\
                    pic.bz[1][0] * short * short +\
                    pic.bz[1][1] * (1 - short) * short)
        pic.interpolate_b()
    
        print("Expected: {} Actual: {}".format(expected, pic.b_particle[0]))
    return

if __name__ == "__main__":
    unittest.main()
