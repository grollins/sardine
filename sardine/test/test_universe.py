from unittest import TestCase
from numpy import array, sqrt, dot, diagflat, allclose
from ..universe import UniverseFactory

class TestUniverse(TestCase):
    def setUp(self):
        uf = UniverseFactory()
        uf.add_atom(x=-1., y=0, z=0, mass=1.0, charge=0.0, radius=1.0, serial_num=1)
        uf.add_atom(x=0., y=0, z=0, mass=1.0, charge=0.0, radius=1.0, serial_num=2)
        self.universe = uf.create_universe()

    def test_coord_matrix_size_is_correct(self):
        X = self.universe.get_coords()
        self.assertEqual(X.ndim, 2)
        self.assertEqual(X.shape[0], len(self.universe))
        self.assertEqual(X.shape[1], 3)
