import main as main
import numpy as np
import unittest

class TestAll(unittest.TestCase):

    def testSwap(self):
        s = np.array([0,1,2,3,4,5,6,7,8,9])
        s = main.swap(s,0,9)
        self.assertTrue(np.array_equal(s, np.array([9,1,2,3,4,5,6,7,8,0])))    

        s = np.array([
            [0, 0],
            [1, 1],
            [2, 2],
            [3, 3]
        ])

        s = main.swap(s, 0, 3)

        ans = np.array([
            [3, 3],            
            [1, 1],
            [2, 2],
            [0, 0]
        ])

        self.assertTrue(np.array_equal, ans)


    def testObj(self):
        self.assertEqual(main.objective(10,8), 50)
        self.assertEqual(main.objective(-10, 8), 50)
        self.assertEqual(main.objective(10,-8), 50)
        self.assertEqual(main.objective(-10,-8), 50)

        self.assertEqual(main.objective(8,10), 58)
        self.assertEqual(main.objective(-8,10), 58)
        self.assertEqual(main.objective(8,-10), 58)
        self.assertEqual(main.objective(-8,-10), 58)

    def testComputeDxDy(self):
        pass

if __name__ == "__main__":
    unittest.main()
