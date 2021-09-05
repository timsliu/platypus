'''
Tests for the ChargeStepDivider helper class
'''

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
        # TODO
        dx = self.charge_divider.dx 
        dy = self.charge_divider.dy 
       
        # list of charge steps that only cross 4 boundaries and should
        # not be broken up
        charge_steps = [
            ChargeStep(0.75 * dx, 0.75 * dy, 1.5 * dx, 0.0 * dy),
        ]
        
        split_steps = [
            [
                ChargeStep(0.75 * dx, 0.75 * dy, 0.75 * dx, 0.0 * dy),
                ChargeStep(1.5 * dx, 0.75 * dy, 0.75 * dx, 0.0 * dy),
            ]
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
        # TODO
        pass



if __name__ == "__main__":
    unittest.main()
