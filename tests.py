import unittest
import copy
from Galois import *

"""Testing for the GaloisPy library"""

# TODO
GF2 = GF(2)
GF3 = GF(3)
GF4 = GF(4)
GF5 = GF(5)
GF7 = GF(7)
GF11 = GF(11)
a = "a"
b = "b"
class TestArithMethods(unittest.TestCase):
    def setUp(self):
        GF2 = GF(2)
        GF3 = GF(3)
        GF4 = GF(4)
        GF5 = GF(5)
        GF7 = GF(7)
        GF11 = GF(11)
        
    def test_add(self):
        self.assertEqual(GF11.add_scalar(0, 0), 0)
        self.assertEqual(GF2.add(1, 1), 0)
        self.assertEqual(GF2.add(0, 1), 1)
        self.assertEqual(GF3.add(1, 2), 0)
        self.assertEqual(GF3.add(5, 30), 2)
        self.assertEqual(GF3.add(-1, 0), 2)
        self.assertEqual(GF7.add(3, 4), 0)
        self.assertEqual(GF7.add(49, 6), 6)
        self.assertEqual(GF4.add(0, 0), 0)
        self.assertEqual(GF4.add(a, a), 0)
        self.assertEqual(GF4.add(1, 1), 0)
        self.assertEqual(GF4.add(b, b), 0)
        self.assertEqual(GF4.add(a, 1), b)
        self.assertEqual(GF4.add(1, a), b)
        self.assertEqual(GF4.add(0, a), a)
        
    def test_multiply(self):
        self.assertEqual(GF2.mult_scalar(0, 0), 0)
        self.assertEqual(GF4.mult_scalar(0, 0), 0)
        self.assertEqual(GF4.mult_scalar(0, a), 0)
        self.assertEqual(GF4.mult_scalar(b, 0), 0)
        self.assertEqual(GF4.mult_scalar(1, 0), 0)
        self.assertEqual(GF7.mult_scalar(1, 7), 0)
        self.assertEqual(GF7.mult_scalar(2, 5), 3)
        self.assertEqual(GF7.mult_scalar(346, 55), 4)
        self.assertEqual(GF4.mult_scalar(1, a), a)
        self.assertEqual(GF4.mult_scalar(a, a), b)
        self.assertEqual(GF4.mult_scalar(b, b), a)
        self.assertEqual(GF4.mult_scalar(b, a), 1)
        with self.assertRaises(ValueError):
            GF4.mult_scalar(a, 10)
        
    def test_mult_inv(self):
        self.assertEqual(GF2.mult_inverse(1), 1)
        self.assertEqual(GF4.mult_scalar(a, GF4.mult_inverse(a)), 1)
        self.assertEqual(GF4.mult_scalar(b, GF4.mult_inverse(b)), 1)
        self.assertEqual(GF4.mult_inverse(1), 1)
        self.assertEqual(GF7.mult_scalar(6, GF7.mult_inverse(6)), 1)
        self.assertEqual(GF7.mult_scalar(5, GF7.mult_inverse(5)), 1)
        self.assertEqual(GF7.mult_scalar(2, GF7.mult_inverse(2)), 1)
        self.assertEqual(GF7.mult_inverse(1), 1)
        self.assertEqual(GF11.mult_scalar(9, GF11.mult_inverse(9)), 1)
        self.assertEqual(GF11.mult_scalar(20, GF11.mult_inverse(9)), 1)
        self.assertEqual(GF11.mult_scalar(9, GF11.mult_inverse(20)), 1)
        with self.assertRaises(ZeroDivisionError):
            GF5.mult_inverse(0)
    
    def test_add_inv(self):
        self.assertEqual(GF2.add_inverse(1), 1)
        self.assertEqual(GF4.add_inverse(1), 1)
        self.assertEqual(GF4.add_inverse(a), a)
        self.assertEqual(GF4.add_inverse(b), b)
        self.assertEqual(GF7.add(5, GF7.add_inverse(5)), 0)
        self.assertEqual(GF7.add([1, 4, 2, 6, 9, 120, -1],
                                 GF7.add_inverse([1, 4, 2, 6, 9, 120, -1])),
                         [0, 0, 0, 0, 0, 0, 0])
    
    def test_vectors(self):
        # Fermat's Little Theorem
        self.assertEqual(GF5.exp_scalar(0, 4), 0)
        self.assertEqual(GF5.exp_scalar(2, 4), 1)
        self.assertEqual(GF5.exp_scalar(3, 4), 1)
        self.assertEqual(GF5.exp_scalar(4, 4), 1)
        self.assertEqual(GF11.exp_scalar(8, 10), 1)
        
        self.assertEqual(GF4.exp_scalar(a, 2), b)
        self.assertEqual(GF4.exp_scalar(a, 3), 1)
        self.assertEqual(GF4.exp_scalar(a, 4), a)
        self.assertEqual(GF4.exp_scalar(b, 3), 1)

class TestRREF(unittest.TestCase):
    def setUp(self):
        GF2 = GF(2)
        GF3 = GF(3)
        GF4 = GF(4)
        GF5 = GF(5)
        GF7 = GF(7)
        GF11 = GF(11)
        
    def test_binary_field(self):
        GF2.verbose = False
        # Inputs
        M_0 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        M_1 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        M_2 = [[0]]
        M_3 = [[1]]
        # Outputs
        M_0o = copy.deepcopy(M_0)
        M_1o = copy.deepcopy(M_1)
        M_2o = copy.deepcopy(M_2)
        M_3o = copy.deepcopy(M_3)
        M_4o = [[1, 0, 1, 0], [0, 1, 0, 0], [1, 0, 1, 0]]
        M_5o = [[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 0, 1]]
        GF2.rref(M_0o)
        GF2.rref(M_1o)
        GF2.rref(M_2o)
        GF2.rref(M_3o)
        GF2.rref(M_4o)
        
        self.assertEqual(M_0o, M_0)
        self.assertEqual(M_1o, M_1)
        self.assertEqual(M_2o, M_2)
        self.assertEqual(M_3o, M_3)
        self.assertEqual(M_4o, [[1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 0]])
        
    def test_prime_field(self):
        M_1o = [[1, 1, 2, 1, 2],
                [1, 0, 1, 1, 0],
                [1, 2, 0, 1, 1],
                [1, 1, 2, 0, 2],
                [2, 2, 1, 2, 1]]
        M_2o = [[9, -2],
                [0, 11]]
        M_3o = [[0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 2, 0]]
        GF3.rref(M_1o)
        GF7.rref(M_2o)
        GF3.rref(M_3o)
        
        self.assertEqual(M_1o, [[1, 0, 1, 0, 0],
                                [0, 1, 1, 0, 2],
                                [0, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0]])
        self.assertEqual(M_2o, [[1, 0], [0, 1]])
        self.assertEqual(M_3o, [[0, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0]])
    
    def test_GF4(self):
        GF4.verbose = True
        M_1o = [[a, b],
                [b, a]]
        M_2o = [[0, 0, b, 0],
                [0, 0, 0, 0],
                [a, 0, b, 1],
                [1, 0, a, b]]
        GF4.rref(M_1o)
        GF4.rref(M_2o)
        
        self.assertEqual(M_1o, [[1, 0], [0, 1]])
        self.assertEqual(M_2o, [[1, 0, 0, b],
                                [0, 0, 1, 0],
                                [0, 0, 0, 0],
                                [0, 0, 0, 0]])
    
    def test_zero(self):
        pass


class TestCodingMethods(unittest.TestCase):
    def setUp(self):
        GF2 = GF(2)
        GF3 = GF(3)
        GF4 = GF(4)
        GF5 = GF(5)
        GF7 = GF(7)
        GF11 = GF(11)
        
    def test_lin_dep(self):
        pass
    
    def test_rref(self):
        pass

if __name__ == '__main__':
    unittest.main()
