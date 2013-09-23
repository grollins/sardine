from unittest import TestCase
from numpy import zeros, zeros_like, allclose, cos, sin, pi, dot
from numpy.linalg import norm
from ..energy import BondEnergyFactory, AngleEnergyFactory, VDWEnergyFactory
from ..bonded_terms import BondFactory, AngleFactory
from ..util import deg2rad, compute_angle

class TestBondGradient(TestCase):
    def setUp(self):
        num_atoms = 2
        r = 1.0
        X = zeros([num_atoms,3])
        D = zeros([num_atoms,num_atoms])
        X[0,0] = 0.0 # x1
        X[1,0] = r   # x2
        D[0,1] = r
        D[1,0] = r
        self.X = X
        self.D = D

        force_const = 0.1
        r_0 = 1.5
        total_F = -force_const * (r - r_0)
        dX = X[1,:] - X[0,:]
        F_atom2 = total_F * (dX / r)
        F_atom1 = -F_atom2

        G_atom1 = -F_atom1
        G_atom2 = -F_atom2
        G = zeros_like(X)
        G[0,:] = G_atom1
        G[1,:] = G_atom2
        self.expected_G = G

        bf = BondFactory()
        bond = bf.create_bond(1, 2, force_const, r_0)
        bond_energy = BondEnergyFactory()
        bond_energy.add_bond(bond)
        self.bond_gradient_func = bond_energy.create_gradient_func(num_atoms)

    def test_computes_correct_gradient(self):
        G = self.bond_gradient_func(self.X, self.D)
        self.assertTrue(allclose(G, self.expected_G),
                        "\n%s\n%s" % (G, self.expected_G))


class TestAngleGradient(TestCase):
    def setUp(self):
        num_atoms = 3
        theta = deg2rad(112.)
        X = zeros([num_atoms,3])
        X[0,0] = -1.0 # x1
        X[2,0] = 1.0 * cos(pi - theta) # x3
        X[2,1] = 1.0 * sin(pi - theta) # y3
        self.X = X

        theta_0 = deg2rad(109.5)
        force_const = 64.71
        total_G = force_const * (theta - theta_0)
        vec12 = X[0,:] - X[1,:]
        vec32 = X[2,:] - X[1,:]
        length_vec12 = norm(vec12)
        length_vec32 = norm(vec32)
        vec12 /= length_vec12
        vec32 /= length_vec32
        cos_theta = cos(theta)

        G1 = total_G * (vec32 - cos_theta * vec12) / length_vec12
        G2 = total_G * (vec12 - cos_theta * vec32) / length_vec32
        G_atom1 = -G1
        G_atom2 = G1 + G2
        G_atom3 = -G2

        G = zeros_like(X)
        G[0,:] = G_atom1
        G[1,:] = G_atom2
        G[2,:] = G_atom3
        self.expected_G = G

        af = AngleFactory()
        angle = af.create_angle(1, 2, 3, force_const, theta_0)
        angle_energy_factory = AngleEnergyFactory()
        angle_energy_factory.add_angle(angle)
        self.angle_gradient_func = angle_energy_factory.create_gradient_func()

    def test_computes_correct_gradient(self):
        G = self.angle_gradient_func(self.X, None)
        self.assertTrue( allclose(G, self.expected_G),
                         "\n%s\n%s" % (G, self.expected_G) )


class TestVDWEnergy(TestCase):
    def setUp(self):
        well_distance = 2.6 # angstroms
        well_depth = 0.1 # kcal/mol

        vdw_energy = VDWEnergyFactory()
        vdw_energy.set_well_distance(well_distance)
        vdw_energy.set_well_depth(well_depth)
        self.vdw_gradient_func = vdw_energy.create_gradient_func()

        num_atoms = 2
        X = zeros([num_atoms,num_atoms])
        X[0,0] = 0.0 # x1
        X[1,0] = well_distance # x2
        D = zeros([num_atoms,num_atoms])
        D[0,1] = well_distance
        D[1,0] = well_distance
        self.X = X
        self.D = D

        self.expected_G = zeros_like(X)

    def test_computes_correct_energy(self):
        G = self.vdw_gradient_func(self.X, self.D)
        self.assertTrue( allclose(G, self.expected_G),
                         "\n%s\n%s" % (G, self.expected_G) )
