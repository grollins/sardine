from unittest import TestCase
from numpy import zeros, allclose, cos, sin, pi
from ..energy import BondEnergyFactory, AngleEnergyFactory, VDWEnergyFactory
from ..bonded_terms import BondFactory, AngleFactory
from ..util import deg2rad

class TestBondEnergy(TestCase):
    def setUp(self):
        force_const = 0.1
        r_0 = 1.5
        r = 1.0
        num_atoms = 2
        self.expected_E = 0.5 * force_const * (r-r_0)**2
        bf = BondFactory()
        bond = bf.create_bond(1, 2, force_const, r_0)
        bond_energy = BondEnergyFactory()
        bond_energy.add_bond(bond)
        self.bond_energy_func = bond_energy.create_energy_func(num_atoms)
        X = zeros([num_atoms,3])
        D = zeros([num_atoms,num_atoms])
        X[0,0] = 0.0 # x1
        X[1,0] = r   # x2
        D[0,1] = r
        D[1,0] = r
        self.X = X
        self.D = D

    def test_computes_correct_energy(self):
        E = self.bond_energy_func(self.X, self.D)
        self.assertTrue(allclose(E, self.expected_E),
                        "%.2f\t%.2f" % (E, self.expected_E))


class TestAngleEnergy(TestCase):
    def setUp(self):
        theta_0 = deg2rad(109.5)
        force_const = 64.71
        num_atoms = 3
        theta = deg2rad(112.)
        self.expected_E = 0.5 * force_const * (theta - theta_0)**2
        af = AngleFactory()
        angle = af.create_angle(1, 2, 3, force_const, theta_0)
        angle_energy = AngleEnergyFactory()
        angle_energy.add_angle(angle)
        self.angle_energy_func = angle_energy.create_energy_func()
        X = zeros([num_atoms,num_atoms])
        X[0,0] = -1.0 # x1
        X[2,0] = 1.0 * cos(pi - theta) # x3
        X[2,1] = 1.0 * sin(pi - theta) # y3
        self.X = X

    def test_computes_correct_energy(self):
        E = self.angle_energy_func(self.X, None)
        self.assertTrue( allclose(E, self.expected_E),
                         "%.2f\t%.2f" % (E, self.expected_E) )


class TestVDWEnergy(TestCase):
    def setUp(self):
        well_distance = 2.6 # angstroms
        well_depth = 0.1 # kcal/mol
        self.expected_E = -well_depth

        vdw_energy = VDWEnergyFactory()
        vdw_energy.set_well_distance(well_distance)
        vdw_energy.set_well_depth(well_depth)
        self.vdw_energy_func = vdw_energy.create_energy_func()

        num_atoms = 2
        X = zeros([num_atoms,num_atoms])
        X[0,0] = 0.0 # x1
        X[1,0] = well_distance # x2
        D = zeros([num_atoms,num_atoms])
        D[0,1] = well_distance
        D[1,0] = well_distance
        self.X = X
        self.D = D

    def test_computes_correct_energy(self):
        E = self.vdw_energy_func(self.X, self.D)
        self.assertTrue( allclose(E, self.expected_E),
                         "%.2f\t%.2f" % (E, self.expected_E) )
  