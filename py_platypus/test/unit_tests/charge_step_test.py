'''
Tests for the ChargeStepDivider helper class
'''
import numpy as np

import unittest
import py_platypus as plat
import py_platypus.utils.params
from py_platypus.utils.charge_step import ChargeStepDivider as ChargeStepDivider
from py_platypus.utils.charge_step import ChargeStep as ChargeStep

class ChargeStepTester(unittest.TestCase):

    def setUp(self):
        dx = [1.0, 1.0]         # size of the cells [y, x]
        cells = [32, 32]        # number of cells rows, cols
        self.charge_divider = ChargeStepDivider(dx, cells) 

    def test_boundaries_touched(self):
        """
        test method for determining which boundaries are touched by a particle
        """
        dx = self.charge_divider.dx 
        dy = self.charge_divider.dy 
       
        # positions to test
        # (TODO) add test cases for ghost boundaries 
        positions = [
            [0.75 * dx, 0.75 * dy], 
            [1.25 * dx, 0.75 * dy], 
            [0.75 * dx, 1.25 * dy], 
            [1.25 * dx, 1.25 * dy],
            [3.10 * dx, 1.75 * dy],
        ] 

        # list of horizontal edges intercepted
        h_expected = [
            [[1, 0], [1, 1]],
            [[1, 0], [1, 1]],
            [[1, 0], [1, 1]],
            [[1, 0], [1, 1]],
            [[2, 2], [2, 3]],
        ]
        
        # list of vertical edges intercepted
        v_expected = [
            [[0, 1], [1, 1]],
            [[0, 1], [1, 1]],
            [[0, 1], [1, 1]],
            [[0, 1], [1, 1]],
            [[1, 3], [2, 3]],
        ]

        # iterate through the test cases
        for i, pos in enumerate(positions):
            with self.subTest(i=i):
                horizontal, vertical = self.charge_divider.boundaries_touched(pos[0], pos[1]) 
                self.assertEqual(horizontal, h_expected[i], "Horizontal mismatch")
                self.assertEqual(vertical, v_expected[i], "Vertical mismatch")

    def test_four_boundaries(self):
        """
        test the case where four boundaries are crossed
        """
        dx = self.charge_divider.dx 
        dy = self.charge_divider.dy 
       
        # list of charge steps that only cross 4 boundaries and should
        # not be broken up
        charge_steps = [
            ChargeStep(0.75 * dx, 0.75 * dy, 0.5 * dx, 0.5 * dy),
            ChargeStep(0.51 * dx, 0.51 * dy, 0.78 * dx, 0.78 * dy),
            ChargeStep(1.6 * dx, 3 * dy, 0.8 * dx, 0.2 * dy),
            ChargeStep(2.25 * dx, 2.25 * dy, -0.5 * dx, -0.5 * dy),
            ChargeStep(2.25 * dx, 2.25 * dy, -0.5 * dx, 0.2 * dy),
        ]
        
        # iterate through the test cases
        for i, c in enumerate(charge_steps):
            with self.subTest(i=i):
                substeps = self.charge_divider.get_charge_steps(c.x0, c.y0, c.x1, c.y1) 
                self.assertEqual(len(substeps), 1, "Only one substep expected")
                self.assertEqual(substeps, [c], "Substeps do not match expectation")

    def test_seven_boundaries(self):
        """
        test the case where seven boundaries are crossed
        """
        dx = self.charge_divider.dx 
        dy = self.charge_divider.dy 
      
        # list of charge steps crossing 7 boundaries
        # note that some of these test cases violate the Courant condition
        # but are still valid tests
        charge_steps = [
            ChargeStep(0.75 * dx, 0.75 * dy, 1.5 * dx, 0.0 * dy),  # right
            ChargeStep(0.8 * dx, 0.75 * dy, 1.5 * dx, 0.15 * dy),  # right
            ChargeStep(0.8 * dx, 0.75 * dy, 0.15 * dx, 1.5 * dy),  # down
            ChargeStep(4.1 * dx, 0.75 * dy, -0.7 * dx, 0.35 * dy), # left
            ChargeStep(3.5 * dx, 5.1 * dy, -0.2 * dx, -0.8 * dy),  # up
        ]
       
        # lists of expected substeps that the charge_steps are divided
        # into
        split_steps = [
            [
                ChargeStep(0.75 * dx, 0.75 * dy, 0.75 * dx, 0.0 * dy),
                ChargeStep(1.5 * dx, 0.75 * dy, 0.75 * dx, 0.0 * dy),
            ],
            [
                ChargeStep(0.8 * dx, 0.75 * dy, 0.7 * dx, 0.07 * dy),
                ChargeStep(1.5 * dx, 0.82 * dy, 0.8 * dx, 0.08 * dy),
            ],
            [
                ChargeStep(0.8 * dx, 0.75 * dy, 0.075 * dx, 0.75 * dy),
                ChargeStep(0.875 * dx, 1.5 * dy, 0.075 * dx, 0.75 * dy),
            ],
            [
                ChargeStep(4.1 * dx, 0.75 * dy, -0.6 * dx, 0.3 * dy),
                ChargeStep(3.5 * dx, 1.05 * dy, -0.1 * dx, 0.05 * dy),
            ],
            [
                ChargeStep(3.5 * dx, 5.1 * dy, -0.15 * dx, -0.6 * dy),
                ChargeStep(3.35 * dx, 4.5 * dy, -0.05 * dx, -0.2 * dy),
            ],
        ]
        
        # iterate through the test cases
        for i, c in enumerate(charge_steps):
            with self.subTest(i=i):
                substeps = self.charge_divider.get_charge_steps(c.x0, c.y0, c.x1, c.y1) 
                self.assertEqual(len(substeps), 2, "Exactly two substeps expected")
                self.assertEqual(substeps, split_steps[i], "Substeps do not match expectation")
    
    def test_ten_boundaries(self):
        """
        test the case where ten boundaries are crossed
        """
        dx = self.charge_divider.dx 
        dy = self.charge_divider.dy 
      
        # list of charge steps crossing 10 boundaries
        # note that some of these test cases violate the Courant condition
        # but are still valid tests
        charge_steps = [
            ChargeStep(0.9 * dx, 0.9 * dy, 0.7 * dx, 0.8 * dy),   # lower right
            ChargeStep(5.1 * dx, 4.1 * dy, -0.8 * dx, -0.7 * dy), # upper left
            ChargeStep(5.1 * dx, 4.8 * dy, -0.8 * dx, 0.8 * dy),  # lower left
            ChargeStep(2.9 * dx, 3.2 * dy, 0.8 * dx, -0.8 * dy),  # upper right
        ]
       
        # lists of expected substeps that the charge_steps are divided
        # into
        split_steps = [
            [
                ChargeStep(0.9 * dx, 0.9 * dy, 0.525 * dx, 0.6 * dy),
                ChargeStep(1.425 * dx, 1.5 * dy, 0.075 * dx, 0.075 * 8/7 * dy),
                ChargeStep(1.5 * dx, (1.5 + 0.075 * 8/7) * dy, 0.1 * dx, (0.2 - 0.075 * 8/7) * dy),
            ],
            [
                ChargeStep(5.1 * dx, 4.1 * dy, -0.6 * dx, -0.525 * dy),
                ChargeStep(4.5 * dx, 3.575 * dy, -0.075 * 8/7 * dx, -0.075 * dy),
                ChargeStep((4.5 - 0.075 * 8/7) * dx, 3.5 * dy, (-0.2 + 0.075 * 8/7)* dx, -0.1 * dy),
            ],
            [
                ChargeStep(5.1 * dx, 4.8 * dy, -0.6 * dx, 0.6 * dy),
                ChargeStep(4.5 * dx, 5.4 * dy, -0.1 * dx, 0.1 * dy),
                ChargeStep(4.4 * dx, 5.5 * dy, -0.1 * dx, 0.1 * dy),
            ],
            [
                ChargeStep(2.9 * dx, 3.2 * dy, 0.6 * dx, -0.6 * dy),
                ChargeStep(3.5 * dx, 2.6 * dy, 0.1 * dx, -0.1 * dy),
                ChargeStep(3.6 * dx, 2.5 * dy, 0.1 * dx, -0.1 * dy),
            ]
        ]
        
        # iterate through the test cases
        for i, c in enumerate(charge_steps):
            with self.subTest(i=i):
                substeps = self.charge_divider.get_charge_steps(c.x0, c.y0, c.x1, c.y1) 
                self.assertEqual(len(substeps), 3, "Exactly three substeps expected")
                self.assertEqual(substeps, split_steps[i], "Substeps do not match expectation")

    def test_intermediate_correct(self):
        """
        Test to make sure that the start of the second charge step
        is between the starting and ending position
        """
        
        dx = [2 * np.pi/32, 2 * np.pi/32]         # size of the cells [y, x]
        cells = [32, 32]        # number of cells rows, cols
        charge_divider = ChargeStepDivider(dx, cells) 
        
        dx = charge_divider.dx 
        dy = charge_divider.dy

        x0 = 2.711
        x1 = 2.714
        y0 = 2.651
        y1 = 2.649

        substeps = charge_divider.get_charge_steps(x0, y0, x1, y1)
        self.assertEqual(len(substeps), 2, "Exactly 2 substeps expected")
        self.assertTrue(substeps[0].x1 > min([x0, x1]), "x1 too low")
        self.assertTrue(substeps[0].x1 < max([x0, x1]), "x1 too high")

if __name__ == "__main__":
    unittest.main()
