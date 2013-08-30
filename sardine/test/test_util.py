from unittest import TestCase
from numpy import zeros, allclose
from numpy.linalg import norm
from ..util import compute_distance_vector, compute_distance_matrix_from_vector

class TestDistanceMatrix(TestCase):
    def setUp(self):
        X = zeros([4,3])
        X[0,0] = -2.0
        X[1,0] = -1.0
        X[2,0] = 0.0
        X[3,0] = 4.0
        D = zeros([4,4])
        d_01 = norm( X[0,:] - X[1,:] )
        d_02 = norm( X[0,:] - X[2,:] )
        d_03 = norm( X[0,:] - X[3,:] )
        d_12 = norm( X[1,:] - X[2,:] )
        d_13 = norm( X[1,:] - X[3,:] )
        d_23 = norm( X[2,:] - X[3,:] )
        D[0,1] = d_01
        D[1,0] = d_01
        D[0,2] = d_02
        D[2,0] = d_02
        D[0,3] = d_03
        D[3,0] = d_03
        D[1,2] = d_12
        D[2,1] = d_12
        D[1,3] = d_13
        D[3,1] = d_13
        D[2,3] = d_23
        D[3,2] = d_23
        self.expected_D = D
        self.X = X

    def test_distance_matrix_is_correct(self):
        D_vec = compute_distance_vector(self.X)
        D = compute_distance_matrix_from_vector(D_vec)
        self.assertTrue(allclose(D, self.expected_D),
                        "\n%s\n%s" % (D, self.expected_D))
