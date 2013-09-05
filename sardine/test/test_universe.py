from unittest import TestCase
from numpy import array, sqrt, dot, diagflat, allclose
from ..universe import UniverseFactory
from ..atom import AtomFactory

class TestUniverse(TestCase):
    def setUp(self):
        af = AtomFactory()
        atom1 = af.create_atom(x=-1., y=0, z=0, mass=1.0, charge=0.0,
                               radius=1.0, serial_num=1)
        atom2 = af.create_atom(x=0., y=0, z=0, mass=1.0, charge=0.0,
                               radius=1.0, serial_num=2)
        uf = UniverseFactory()
        uf.add_atom(atom1)
        uf.add_atom(atom2)
        self.universe = uf.create_universe()

    def test_coord_matrix_size_is_correct(self):
        X = self.universe.get_coords()
        self.assertEqual(X.ndim, 2)
        self.assertEqual(X.shape[0], len(self.universe))
        self.assertEqual(X.shape[1], 3)
