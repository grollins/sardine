from unittest import TestCase
from numpy import zeros, allclose
from ..energy import BondEnergyFactory

class TestBondEnergy(TestCase):
    def setUp(self):
        force_const = 0.1
        r_0 = 1.5
        r = 1.0
        num_atoms = 2
        self.expected_E = 0.5 * force_const * (r-r_0)**2
        bef = BondEnergyFactory()
        bef.add_term(1, 2, force_const, r_0)
        self.bond_energy_func = bef.create_func(num_atoms)
        D = zeros([num_atoms,num_atoms])
        D[0,1] = r
        D[1,0] = r
        self.D = D

    def test_computes_correct_energy(self):
        E = self.bond_energy_func(self.D)
        self.assertTrue(allclose(E, self.expected_E),
                        "%.2f\t%.2f" % (E, self.expected_E))
