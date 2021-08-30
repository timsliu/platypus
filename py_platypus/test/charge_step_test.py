'''
Tests for the ChargeStepDivider helper class
'''

import unittest
import py_platypus.utils.params
import py_platypus.utils.current_helpers as current_helpers


class ChargeStepTester(unittest.TestCase):

    def setUp(self):
        dx = [0.5, 0.5]         # size of the cells [y, x]
        cells = [32, 32]        # number of cells rows, cols
        self.charge_divider = current_helpers.ChargeStepDivider(dx, cells) 

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
        # TODO 
        pass

    def test_seven_boundaries(self):
        """
        test the case where seven boundaries are crossed
        """
        # TODO
        pass
    
    def test_ten_boundaries(self):
        """
        test the case where ten boundaries are crossed
        """
        # TODO
        pass



if __name__ == "__main__":
    unittest.main()
