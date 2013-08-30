from unittest import TestCase
from numpy import array, sqrt, dot, diagflat, allclose
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, EnergyFunctionFactory
from ..nma import compute_hessian, compute_force_constant_matrix

m1 = 1.0
m2 = 1.0
m3 = 1.0
k = 0.1
r_0 = 1.0

uf = UniverseFactory()
uf.add_atom(x=-1., y=0, z=0, mass=m1, charge=0.0, radius=1.0, serial_num=1)
uf.add_atom(x=0., y=0, z=0, mass=m2, charge=0.0, radius=1.0, serial_num=2)
uf.add_atom(x=1., y=0, z=0, mass=m3, charge=0.0, radius=1.0, serial_num=3)
universe = uf.create_universe()

bef = BondEnergyFactory()
bef.add_term(1, 2, k, r_0)
bef.add_term(2, 3, k, r_0)
bond_energy_func = bef.create_func(num_atoms=3)
eff = EnergyFunctionFactory()
eff.add_energy_term('bonds', bond_energy_func)
energy_func = eff.create_energy_func(['bonds'])

expected_H = array(([k, -k, 0],
                    [-k, 2*k, -k],
                    [0, -k, k]))

M_diag = array((m1, m2, m3))
inv_sqrt_diag = 1./sqrt(M_diag)
expected_M = diagflat(inv_sqrt_diag)

expected_F = dot( dot(expected_M, expected_H), expected_M )


class TestTriatomicNMA(TestCase):

    def test_computes_mass_matrix(self):
        M = universe.get_inv_sqrt_mass_matrix()
        self.assertTrue( allclose(M, expected_M),
                         msg="\n%s\n%s" % (M, expected_M) )

    def test_computes_hessian(self):
        X = universe.get_coords()
        H = compute_hessian(energy_func, X)
        self.assertTrue( allclose(H, expected_H),
                         msg="\n%s\n%s" % (H, expected_H) )

    def test_computes_force_constant_matrix(self):
        F = compute_force_constant_matrix(expected_H, expected_M)
        self.assertTrue( allclose(F, expected_F),
                         msg="\n%s\n%s" % (F, expected_F) )
