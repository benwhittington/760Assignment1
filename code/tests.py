import main as main
import numpy as np
import unittest

class TestAll(unittest.TestCase):

    def testSwap(self):
        s = np.array([0,1,2,3,4,5,6,7,8,9])
        s = main.swap(s,0,9)
        self.assertTrue(np.array_equal(s, np.array([9,1,2,3,4,5,6,7,8,0])))        
    

if __name__=="__main__":
    unittest.main()
